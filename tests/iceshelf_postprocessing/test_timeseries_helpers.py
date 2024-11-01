import numpy as np
import pandas as pd
import pytest

from polar_iceshelves.iceshelf_postprocessing import timeseries_helpers


@pytest.fixture()
def example_timeseries_dataframe() -> pd.DataFrame:
    np.random.seed(2501)
    dates = np.arange(2020, 2022, 1 / 12)
    noise = np.random.rand(len(dates)) * 100
    signal = np.linspace(1000, 500, len(dates))
    # set a few items nan to test handling
    signal[3] = np.nan
    signal[7] = np.nan
    return pd.DataFrame({
        'dates_decimal': dates,
        'changes': noise + signal})


def test_get_timeseries_rate(example_timeseries_dataframe: pd.DataFrame):
    np.testing.assert_almost_equal(timeseries_helpers.get_timeseries_rate(
        example_timeseries_dataframe), -273.773, decimal=3)


def test_match_timeseries_results_to_other_dates(example_timeseries_dataframe: pd.DataFrame):
    np.testing.assert_almost_equal(timeseries_helpers.match_timeseries_results_to_other_dates(
        [2021.4, 2021.75], example_timeseries_dataframe, column='changes'), [672.245, 571.394], decimal=3)


def test_calculate_smooth_changes(example_timeseries_dataframe: pd.DataFrame):
    np.testing.assert_almost_equal(timeseries_helpers.calculate_smooth_changes(
        example_timeseries_dataframe)[[5, 10, 15]], [942.392, 840.353, 723.231], decimal=3)


def test_get_mean_with_outliers_removed(example_timeseries_dataframe: pd.DataFrame):
    np.testing.assert_almost_equal(timeseries_helpers.get_mean_with_outliers_removed(
        example_timeseries_dataframe['changes']), 789.569, decimal=3)


def test_replace_nans_with_interpolation(example_timeseries_dataframe: pd.DataFrame):
    nan_free_measurements = timeseries_helpers.replace_nans_with_interpolation(
        example_timeseries_dataframe['changes'])
    assert not any(np.isnan(nan_free_measurements))
    assert len(np.unique(nan_free_measurements)) == len(nan_free_measurements)


def test_replace_nan_padding_with_closest_value():
    test_data = np.array([np.nan, np.nan, 1, 2, np.nan, 4, 5, np.nan])
    de_padded_test_data = timeseries_helpers.replace_nan_padding_with_closest_value(test_data)
    np.testing.assert_array_equal(de_padded_test_data, [1, 1, 1, 2, np.nan, 4, 5, 5])


def test_set_outliers_to_nan(example_timeseries_dataframe: pd.DataFrame):
    # add an outlier
    example_timeseries_dataframe['changes'][13] = 1e6

    # it doesn't make sense to apply this to time
    result = timeseries_helpers.set_outliers_to_nan(example_timeseries_dataframe.drop(columns=['dates_decimal']))

    # check it's now NaN
    assert np.isnan(result['changes'][13])
