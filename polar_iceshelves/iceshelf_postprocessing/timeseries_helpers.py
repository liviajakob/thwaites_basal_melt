from datetime import datetime
import glob
import os

import numpy as np
import pandas as pd
import rasterio
from rasterstats import zonal_stats
from scipy.signal import savgol_filter
import statsmodels.api as sm
import xarray as xr

from polar_iceshelves.iceshelf_postprocessing.date_helpers import datetime2year
from polar_iceshelves.iceshelf_postprocessing.date_helpers import year2datetime
from polar_iceshelves.iceshelf_postprocessing import constants
from shapely.errors import ShapelyDeprecationWarning
import warnings
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)


def generate_timeseries_from_tifs(in_filepath, polygon, nodata=None, calculate_coverage=False):
    timeseries_tifs = sorted(glob.glob(in_filepath + "/*.tif"))
    dates = [datetime.strptime(os.path.splitext(os.path.split(f)[1])[0], '%Y-%m-%d') for f in timeseries_tifs]
    # calculate zonal stats
    averages = []
    coverages = []
    for f in timeseries_tifs:
        stats = zonal_stats(polygon, f, copy_properties=True, stats=['mean'], nodata=nodata)
        averages.append(stats[0]['mean'])
        if calculate_coverage:
            stats_covered = zonal_stats(polygon, f, copy_properties=True, stats=['count'], nodata=nodata)
            stats_all = zonal_stats(polygon, f, copy_properties=True, stats=['count'], nodata=1)
            coverages.append(stats_covered[0]["count"] / stats_all[0]["count"])

    df_timeseries = pd.DataFrame()
    df_timeseries['dates'] = dates
    df_timeseries['dates_datetime'] = dates
    df_timeseries['dates_decimal'] = [round(datetime2year(d), 6) for d in df_timeseries.dates_datetime]
    df_timeseries['averages'] = averages
    df_timeseries['changes'] = df_timeseries['averages'] - df_timeseries['averages'].iloc[0]
    if calculate_coverage:
        df_timeseries['coverages'] = coverages
    return df_timeseries


def generate_errors_from_tifs(in_filepath, polygon, nodata=None):
    timeseries_tifs = sorted(glob.glob(in_filepath + "/*.tif"))
    dates = [datetime.strptime(os.path.splitext(os.path.split(f)[1])[0], '%Y-%m-%d') for f in timeseries_tifs]
    # calculate zonal stats
    errors = []
    for f in timeseries_tifs:
        stats = zonal_stats(polygon, f, copy_properties=True, add_stats={
                            'errors': zonal_stats_errors_calc}, nodata=nodata)
        errors.append(stats[0]['errors'])
    df_timeseries = pd.DataFrame()
    df_timeseries['dates'] = dates
    df_timeseries['dates_datetime'] = dates
    df_timeseries['dates_decimal'] = [round(datetime2year(d), 6) for d in df_timeseries.dates_datetime]
    df_timeseries['errors'] = errors
    return df_timeseries


def generate_timeseries_from_netcdf(in_filepath, polygon, column='elev_c', date_column='time'):
    f = rasterio.open(in_filepath)
    affine = f.transform
    # load and read netCDF-file to dataset and get datarray for variable
    nc_fo = in_filepath
    nc_ds = xr.open_dataset(nc_fo)
    nc_var = nc_ds[column]
    years = nc_ds[date_column].values
    # go through all years
    averages = []
    for year in years:
        # get values of variable pear year
        nc_arr = nc_var.sel(time=year)
        nc_arr_vals = nc_arr.values
        print(affine)
        # go through all geometries and compute zonal statistics
        stats = zonal_stats(polygon, nc_arr_vals, affine=affine, stats=['mean'], nodata=f.nodata)
        averages.append(stats[0]['mean'])
    df_timeseries = pd.DataFrame()
    datetime_dates = ([pd.Timestamp(year2datetime(d)).round(freq='D') for d in years])
    df_timeseries['dates'] = datetime_dates
    df_timeseries['dates_datetime'] = datetime_dates
    df_timeseries['dates_decimal'] = years
    df_timeseries['averages'] = averages
    df_timeseries['changes'] = df_timeseries['averages'] - df_timeseries['averages'].iloc[0]
    return df_timeseries


def get_timeseries_rate(df_timeseries, column='changes'):
    '''
    :param df_timeseries: dataframe to calculate rate on (will have dates_datetime and timeseries columns),
    :param column: changes column in df_timeseries
    :return: timeseries rate per year (float)
    '''
    df_timeseries = df_timeseries[df_timeseries[column].notna()]
    vals = np.asarray([df_timeseries['dates_decimal']])
    vals = sm.add_constant(np.transpose(vals))
    ts = df_timeseries[column]
    # Create model and fit it (least squares)
    model = sm.RLM(ts, vals)
    results = model.fit()
    # ts_rate_yearly = results.params.x1*31556926.08
    ts_rate_yearly = results.params.x1
    return ts_rate_yearly


def match_timeseries_results_to_other_dates(dates: list, df_timeseries: pd.DataFrame, column: str, dayspan: int = 90):
    '''match a timeseries to a different set of dates and compute averages of a column over specified timespan
    :param dates: dates to be matched to
    :param df_timeseries: pandas Dataframe with a timeseries, must contain dates_decimal column
    :param column: columns from df_timeseries to be used for matching calc
    :param dayspan: span to consider for generating an average value, number of days
    :return:
    '''
    averages = []
    span_decimal = (dayspan / 2) / 365.25
    for d in dates:
        startdate = d - span_decimal
        enddate = d + span_decimal
        dfsub = df_timeseries[(df_timeseries['dates_decimal'] >= startdate)
                              & (df_timeseries['dates_decimal'] < enddate)]
        averages.append(dfsub[column].mean())
    return averages


def calculate_smooth_changes(df_timeseries, window_size=9, poly_order=1, column='changes'):
    arr = np.array(df_timeseries[column])
    nan_mask = np.isnan(arr)
    if not nan_mask.all():
        arr[nan_mask] = np.interp(np.flatnonzero(nan_mask), np.flatnonzero(~nan_mask), arr[~nan_mask])  # fill no data
        smoothed = savgol_filter(arr, window_size, poly_order)
        return smoothed
        # smoothed_start_value = smoothed[0]
        # return smoothed-smoothed_start_value
    else:
        return np.nan


def get_mean_with_outliers_removed(arr):
    arr = arr[~np.isnan(arr)]
    mean = np.mean(arr, axis=0)
    sd = np.std(arr, axis=0)
    final_list = [x for x in arr if (x > mean - 2 * sd)]
    final_list = [x for x in final_list if (x < mean + 2 * sd)]
    return np.mean(final_list)


def zonal_stats_errors_calc(arr, neff_division=constants.SPATIAL_AUTOCORRELATION_NEFF):
    ''' Helper for zonal_stats to calculate errors from a number of pixels
    '''
    count = np.ma.count(arr)
    if count > 0:
        # sum_of_square_of_arr = np.sum(np.square(np.array(arr, dtype=np.float128)))
        sum_of_square_of_arr = np.ma.sum(np.square(arr))
        sqrt = np.sqrt(sum_of_square_of_arr)
        errors = sqrt / (count / neff_division)
    else:
        errors = None
    return errors


def replace_nans_with_interpolation(data):
    ''' Linear interpolation
    '''
    mask = np.isnan(data)
    data[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), data[~mask])
    return data


def replace_nan_padding_with_closest_value(a):
    '''This replaces nan values _at the start and end of the array only_ with the closest value to them.
    '''
    ind = np.where(~np.isnan(a))[0]
    first, last = ind[0], ind[-1]
    a[:first] = a[first]
    a[last + 1:] = a[last]
    return a


def set_outliers_to_nan(df: pd.DataFrame):
    for _, col in df.items():
        # Quartile
        q1 = col.describe(percentiles=[0.05])["5%"]
        q3 = col.describe(percentiles=[0.95])["95%"]
        iqr = q3 - q1  # Interquartile range

        # Outlier reference point
        outlier_min = q1 - (iqr) * 1.5
        outlier_max = q3 + (iqr) * 1.5

        # Excludes values that are out of range
        col[col < outlier_min] = None
        col[col > outlier_max] = None
    return df
