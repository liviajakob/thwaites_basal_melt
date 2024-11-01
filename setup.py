"""Setup file for polar-iceshelves."""
import os

from setuptools import find_packages
from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

if os.path.exists("full_version.txt"):
    with open("full_version.txt", "r", encoding="utf-8") as fh:
        """
        Note that this file is generated by the CI chain based on the git tag
        (by .github/scripts/define_new_version_number.py)
        It should not be present in the repository by default.
        """
        version_number = fh.read()
else:
    version_number = 'v0.0.0'  # default value when under development

setup(
    name='thwaites_basal_melt',
    version=version_number,
    description='Package containing work completed for the Polar Iceshelves project.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Livia Jakob',
    author_email='livia@earthwave.co.uk',
    url='https://github.com/liviajakob/thwaites_basal_melt',
    python_requires=">=3.9",
    packages=find_packages(),
    # note requirements listed ininstall_requires should be the *minimum required*
    # in order to allow pip to resolve multiple installed packages properly.
    # requirements.txt should contain a specific known working version instead.
    install_requires=[
        'gdal',
        'geopandas',
        'matplotlib',
        'netcdf4',
        'numpy',
        'pandas',
        'rasterio',
        'rasterstats',
        'shapely',
        'statsmodels',
        'xarray'
    ],
)
