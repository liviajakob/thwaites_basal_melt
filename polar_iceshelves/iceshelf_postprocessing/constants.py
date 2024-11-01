import os

# density values
DENSITY_OF_SEA_WATER_KG_PER_M3 = 1028
DENSITY_OF_ICE_KG_PER_M3 = 917

# Uncertainty of elevation change for radar penetration
ELEVATION_CHANGE_PENETRATION_UNCERTAINTY = 0.3

# neff for assumption on spation autocorrelation
SPATIAL_AUTOCORRELATION_NEFF = 2

# File paths
BASE_PATH = '/data/ox1/working/polarice/outputs'
STATIC_MASKS_FP = os.path.join(
    BASE_PATH, 'shapefiles/ice-shelves-static/intersection-mask-updated-with-bens-new-masks-thwaites-new.gpkg')
STATIC_MASK_TIF = os.path.join(BASE_PATH, 'shapefiles/intersection-mask-raster/intersection-mask-updated-gl.tif')
BASAL_MELT_TIF = os.path.join(
    BASE_PATH, 'shapefiles/thwaites_datasets/basal_melt/filter_melt_0_9_shear_0_4_medians_5_thwaites_replaced.tif')
BASAL_MELT_ERROR_TIF = os.path.join(
    BASE_PATH, 'shapefiles/basal_melt_map_firn_calc/final/AIS_basal_melt_errors.tif')
SMB_FILEPATH = os.path.join(BASE_PATH, 'shapefiles/smb-models/smb_racmo_updated/interpolated')
SMB_NO_DATA = None
FIRN_FILEPATH = os.path.join(BASE_PATH, 'shapefiles/firn-air/firn_air_updated/interpolated')
FIRN_NO_DATA = None

# savgol filters
FILTER_1 = 15
FILTER_2 = 8
