import numpy as np
from rasterstats import zonal_stats
from scipy.signal import savgol_filter

from polar_iceshelves.iceshelf_postprocessing import constants
from polar_iceshelves.iceshelf_postprocessing import timeseries_helpers


def calculate_basal_melt_timeseries_with_errors(
        elevation_change_folder_path,
        elevation_change_errors_folder_path,
        elevation_change_errors_folder_path_poca,
        polygon,
        filter_type=None):
    # 1 calculate mean melt rate
    melt_rate = zonal_stats(polygon, constants.BASAL_MELT_TIF, copy_properties=True, stats=['mean'])[0]['mean']
    if melt_rate is None or melt_rate is np.nan:
        return None
    melt_rate_error = zonal_stats(polygon, constants.BASAL_MELT_ERROR_TIF, copy_properties=True,
                                  add_stats={'errors': timeseries_helpers.zonal_stats_errors_calc})[0]['errors']
    print('Melt rate is: ', melt_rate, "+-", melt_rate_error)

    # 2 dhdt timeseries
    print("Calculating dhdt timeseries...")
    nodata = None
    df_timeseries = timeseries_helpers.generate_timeseries_from_tifs(
        elevation_change_folder_path, polygon=polygon, nodata=nodata)
    df_timeseries['dhdt'] = df_timeseries['averages']
    if (df_timeseries["dhdt"].notna()).value_counts().loc[True] < 10:  # don't bother when less than 10
        return None
    # 2.1 dhdt_anomalies
    # gets yearly rate based on fit on elevation change timeseries
    dhdt_rate = timeseries_helpers.get_timeseries_rate(df_timeseries, column='dhdt')
    df_timeseries['dhdt_linear'] = (df_timeseries['dates_decimal'] - df_timeseries['dates_decimal'].iloc[0]) * dhdt_rate
    df_timeseries['dhdt_anomalies'] = df_timeseries['dhdt'] - df_timeseries['dhdt_linear']  # remove rate
    mean_for_adjustment = timeseries_helpers.get_mean_with_outliers_removed(
        np.array(df_timeseries['dhdt_anomalies']))  # adjust heigh of timeseries
    df_timeseries['dhdt_anomalies'] = df_timeseries['dhdt_anomalies'] - mean_for_adjustment
    # remove outliers
    ol_removed = timeseries_helpers.set_outliers_to_nan(df_timeseries[['dhdt_anomalies']])
    ol_interp = timeseries_helpers.replace_nans_with_interpolation(np.array(ol_removed['dhdt_anomalies']))
    df_timeseries['dhdt_anomalies'] = ol_interp

    # 2.2 dhdt errors
    dhdt_errors = timeseries_helpers.generate_errors_from_tifs(
        elevation_change_errors_folder_path, polygon, nodata=None)
    dhdt_errors_poca = timeseries_helpers.generate_errors_from_tifs(
        elevation_change_errors_folder_path_poca, polygon, nodata=None)
    poca_error = np.array(dhdt_errors_poca['errors'])
    poca_error[[isinstance(x, np.ma.core.MaskedConstant) for x in poca_error]] = 0.0
    swath_error = np.array(dhdt_errors['errors'])
    if len(poca_error) < len(swath_error):
        poca_error = np.concatenate([poca_error, np.zeros((len(swath_error) - len(poca_error)))], axis=0)
    elif len(poca_error) > len(swath_error):
        poca_error = poca_error[:len(swath_error)]
    swath_error[[isinstance(x, np.ma.core.MaskedConstant) for x in swath_error]] = 0.0
    swath_error = np.array(list(swath_error))  # hack to sort out issues with inconsistent data types in np array
    poca_error = np.array(list(poca_error))  # hack to sort out issues with inconsistent data types in np array
    dhdt_errors['dhdt_errors'] = (swath_error**2 + poca_error**2)**0.5
    # add penetration uncertainty
    dhdt_errors['dhdt_errors'] = (dhdt_errors['dhdt_errors']**2
                                  + constants.ELEVATION_CHANGE_PENETRATION_UNCERTAINTY**2)**0.5
    df_timeseries = df_timeseries.merge(dhdt_errors[['dates_decimal', 'dhdt_errors']], on='dates_decimal', how='left')
    # 3.1 SMB timeseries
    print("Calculating SMB timeseries...")
    df_timeseries_smb = timeseries_helpers.generate_timeseries_from_tifs(
        constants.SMB_FILEPATH, polygon=polygon, nodata=constants.SMB_NO_DATA)  # This timeseries is in kg/m2
    df_timeseries_smb['changes'] = df_timeseries_smb['averages']
    df_timeseries_smb['changes'] = df_timeseries_smb['changes'].cumsum() - df_timeseries_smb['changes'].iloc[0]

    df_timeseries['smb'] = timeseries_helpers.match_timeseries_results_to_other_dates(
        list(df_timeseries.dates_decimal), df_timeseries_smb, 'changes')

    # 3.2 SMB anomalies
    # gets yearly rate based on fit on smb timeseries
    smb_rate = timeseries_helpers.get_timeseries_rate(df_timeseries, column='smb')
    df_timeseries['smb_linear'] = (df_timeseries['dates_decimal'] - df_timeseries['dates_decimal'].iloc[0]) * smb_rate
    df_timeseries['smb_anomalies'] = df_timeseries['smb'] - df_timeseries['smb_linear']
    arr = np.array(df_timeseries['smb_anomalies'])
    # adjust height of timeseries
    mean_for_adjustment = timeseries_helpers.get_mean_with_outliers_removed(arr[~np.isnan(arr)])
    df_timeseries['smb_anomalies'] = df_timeseries['smb_anomalies'] - mean_for_adjustment
    df_timeseries['smb_anomalies'] = timeseries_helpers.replace_nan_padding_with_closest_value(
        np.array(df_timeseries['smb_anomalies']))

    # 3.3 SMB errors
    df_timeseries['smb_errors'] = abs(np.insert(np.diff(df_timeseries['smb_anomalies']) * 0.1, 0, 0))

    # 4.1 Firn timeseries
    print("Calculating Firn air content timeseries...")
    # THIS DATASET IS IN [m], it describes contents of firn air in m, meaning it's a cumulative time series
    df_timeseries_firn = timeseries_helpers.generate_timeseries_from_tifs(
        constants.FIRN_FILEPATH, polygon=polygon, nodata=constants.FIRN_NO_DATA)
    df_timeseries['firn'] = timeseries_helpers.match_timeseries_results_to_other_dates(
        list(df_timeseries.dates_decimal), df_timeseries_firn, 'changes')
    # 4.2 FIRN anomalies
    # gets yearly rate based on fit on firn timeseries
    firn_rate = timeseries_helpers.get_timeseries_rate(df_timeseries, column='firn')
    df_timeseries['firn_linear'] = (df_timeseries['dates_decimal'] - df_timeseries['dates_decimal'].iloc[0]) * firn_rate
    df_timeseries['firn_anomalies'] = df_timeseries['firn'] - df_timeseries['firn_linear']
    arr = np.array(df_timeseries['firn_anomalies'])
    # adjust height of timeseries
    mean_for_adjustment = timeseries_helpers.get_mean_with_outliers_removed(arr[~np.isnan(arr)])
    df_timeseries['firn_anomalies'] = df_timeseries['firn_anomalies'] - mean_for_adjustment
    df_timeseries['firn_anomalies'] = timeseries_helpers.replace_nan_padding_with_closest_value(
        np.array(df_timeseries['firn_anomalies']))

    # 4.3 Firn errors
    df_timeseries['firn_errors'] = abs(np.insert(np.diff(df_timeseries['firn_anomalies']) * 0.1, 0, 0))

    # BASAL MELT
    print("Calculating basal melt timeseries...")
    # dh_air = 0 ## with firn compaction excluded
    dh_air = df_timeseries['firn_anomalies']  # in meters
    dh_anom = df_timeseries['dhdt_anomalies']  # in meters
    ms_anom = df_timeseries['smb_anomalies']  # in kg/m2
    velocity_anom = 0  # at the moment we don't include velocities

    density_ratio = constants.DENSITY_OF_SEA_WATER_KG_PER_M3 / (
        constants.DENSITY_OF_SEA_WATER_KG_PER_M3 - constants.DENSITY_OF_ICE_KG_PER_M3)

    w_b_anomalies = (
        ms_anom / constants.DENSITY_OF_ICE_KG_PER_M3) - (density_ratio * (dh_anom - dh_air)) - velocity_anom
    df_timeseries['bm_anomalies'] = w_b_anomalies

    # remove first few measurements and cut off
    df_timeseries.drop(df_timeseries.index[:2], inplace=True)  # remove first few measurements
    df_timeseries = df_timeseries.reset_index(drop=True)
    df_timeseries.drop(df_timeseries.index[-2:], inplace=True)  # remove last few measurements
    df_timeseries = df_timeseries.reset_index(drop=True)

    df_timeseries['bm_linear'] = (df_timeseries['dates_decimal'] - df_timeseries['dates_decimal'].iloc[0]) * melt_rate
    df_timeseries['bm_cum'] = df_timeseries['bm_linear'] + df_timeseries['bm_anomalies']

    # Error propagation
    smb_err = df_timeseries['smb_errors']
    dhdt_err = df_timeseries['dhdt_errors']
    firn_err = df_timeseries['firn_errors']
    df_timeseries['bm_anomalies_errors'] = (
        (smb_err / constants.DENSITY_OF_ICE_KG_PER_M3)**2 + (
            density_ratio * (((dhdt_err**2) + (firn_err**2))**0.5))**2)**0.5
    df_timeseries = apply_filters(df_timeseries, 'bm_anomalies_errors', filter_type=filter_type)

    # final bm timeseries yearly and monthly melt rates
    df_timeseries['bm_anomaly_derivative'] = df_timeseries['bm_anomalies'].diff()
    df_timeseries['bm_monthly'] = melt_rate + 12 * df_timeseries['bm_anomaly_derivative']
    df_timeseries = apply_filters(df_timeseries, 'bm_monthly', filter_type=filter_type)

    df_timeseries['bm_monthly_errors'] = np.sqrt(
        (melt_rate_error / 12)**2
        + np.sqrt(
            (df_timeseries['bm_anomalies_errors'])**2
            + np.concatenate(
                ([0], np.array(df_timeseries['bm_anomalies_errors'])[:-1]))**2)**2)  # mean melt rate, t-1 and t errors

    return df_timeseries


def apply_filters(df, column, filter_type):
    arr = np.array(df[column])
    nan_mask = np.isnan(arr)
    arr[nan_mask] = np.interp(np.flatnonzero(nan_mask), np.flatnonzero(~nan_mask), arr[~nan_mask])  # fill no data
    df[column] = savgol_filter(arr, 5, 1)
    arr = np.array(df[column])
    nan_mask = np.isnan(arr)
    arr[nan_mask] = np.interp(np.flatnonzero(nan_mask), np.flatnonzero(~nan_mask), arr[~nan_mask])  # fill no data
    if filter_type == "event":
        df[column] = [*savgol_filter(arr, constants.FILTER_2, 1)[:50], *savgol_filter(arr, constants.FILTER_1, 1)[50:]]
    else:
        df[column] = savgol_filter(arr, constants.FILTER_1, 1)
    return df
