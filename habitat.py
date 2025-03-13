# Load packages
import os
from glob import glob
import pathlib
from getpass import getpass
from math import floor, ceil
import time
import zipfile

# Third-party imports
import contextily as ctx
import earthaccess
import geopandas as gpd
import holoviews as hv
import hvplot.pandas
import hvplot.xarray
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
import pygbif.occurrences as occ
import pygbif.species as species
import requests
import regionmask
import rioxarray as rxr
import rioxarray.merge as rxrm
from shapely.geometry import MultiPolygon, Polygon
import skfuzzy
import xrspatial
import xarray as xr
import xrspatial as xrs

# make reproducible file paths
data_dir = os.path.join(
    pathlib.Path.home(),
    'earth-analytics',
    'data',
    'habitat'
)

# make the directory
os.makedirs(data_dir, exist_ok=True)

# set gbif directory
gbif_dir = os.path.join(data_dir, 'gbif_cholla')

# Access gbif
reset_credentials = False

# Enter gbif username, password, and email
credentials = dict(
    GBIF_USER=(input, 'GBIF username:'),
    GBIF_PWD=(getpass, 'GBIF password: '),
    GBIF_EMAIL=(input, 'GBIF email: '),
)
for env_variable, (prompt_func, prompt_text) in credentials.items():

    # Delete credential from the environment if requested
    if reset_credentials and (env_variable in os.environ):
        os.environ.pop(env_variable)

    # Ask for credential and save to environment
    if not env_variable in os.environ:
        os.environ[env_variable] = prompt_func(prompt_text)

# assign species code
species_key = 7282673

# set a file pattern
gbif_pattern = os.path.join(gbif_dir,
                            '*.csv')

# Download it once
if not glob(gbif_pattern):
    # Submit my query to GBIF
    gbif_query = occ.download([
        f"speciesKey = {species_key}",
        "hasCoordinate = True",
    ])

    # Only download once
    if not 'GBIF_DOWNLOAD_KEY' in os.environ:
        os.environ['GBIF_DOWNLOAD_KEY'] = gbif_query[0]
        download_key = os.environ['GBIF_DOWNLOAD_KEY']

        # Wait for download to build
        wait = occ.download_meta(download_key)['status']
        while not wait == 'SUCCEEDED':
            wait = occ.download_meta(download_key)['status']
            time.sleep(5)

    # Download the data
    download_info = occ.download_get(
        os.environ['GBIF_DOWNLOAD_KEY'],
        path=data_dir
    )

    # Unzip it
    with zipfile.ZipFile(download_info['path']) as download_zip:
        download_zip.extractall(path=gbif_dir)


# Find csv file path
gbif_path = glob(gbif_pattern)[0]

# Open gbif data
gbif_df = pd.read_csv(
    gbif_path,
    delimiter='\t'
)

# Make it spatial
gbif_gdf = (
    gpd.GeoDataFrame(
        gbif_df,
        geometry=gpd.points_from_xy(
            gbif_df.decimalLongitude,
            gbif_df.decimalLatitude
        ),
        crs='EPSG:4326'
    )
)

# Reproject GeoDataFrame to the same CRS as the basemap (Web Mercator, EPSG:3857)
gbif_gdf_merc = gbif_gdf.to_crs(epsg=3857)

# Create a plot
ax = gbif_gdf_merc.plot(figsize=(10, 10), edgecolor='black', alpha=0.7)
ctx.add_basemap(ax, crs=gbif_gdf_merc.crs.to_string(),
                source=ctx.providers.Esri.WorldImagery)
plt.title("Jumping Cholla Occurances in GBIF")
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# Show the plot
plt.show()

# URL for the dataset (direct link) note this url is not working and I can't find the right one you have to go to the catalog and directly download it to your machine. 
# https://www.usgs.gov/programs/gap-analysis-project/science/pad-us-data-download

pa_url = "https://www.sciencebase.gov/catalog/item/652d4f80d34e44db0e2ee45c"

# Download the file
pa_zip_path = os.path.join(data_dir, "PADUS4_0_State_AZ_GDB.zip")
if not os.path.exists(pa_zip_path):
    response = requests.get(pa_url, stream=True)
    if response.status_code == 200:
        with open(pa_zip_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=128):
                f.write(chunk)
        print("Downloaded Zip")
    else:
        print(f"Failed to download file. Status code: {response.status_code}")

# Unzip file
pa_dir = os.path.join(data_dir, "sites_cholla")
if not os.path.exists(pa_dir):
    with zipfile.ZipFile(pa_zip_path, 'r') as zip_ref:
        zip_ref.extractall(pa_dir)
    print(f"Extracted data to {pa_dir} directory")

## open pa path
pa_path = os.path.join(pa_dir, 'PADUS4_0_StateAZ.gdb')

# open polygon
pa_shp = gpd.read_file(pa_path)

# convert crs
pa_shp = pa_shp.to_crs(epsg=4326)

# fix invalid geoms
pa_shp['geometry'] = pa_shp['geometry'].apply(lambda geom: geom.make_valid() if not isinstance(geom,
                                                                                              MultiPolygon) and not geom.is_valid else geom)

# drop remaining invalid geometries
pa_shp = pa_shp[pa_shp.geometry.is_valid]

# drop rows with missing geometries
pa_shp = pa_shp.dropna(subset=['geometry'])

# simplify columns
pa_shp = pa_shp[['Own_Name', 'Mang_Name',
                 'Unit_Nm', 'Loc_Nm',
                 'geometry']]

# intersect cholla  occurrence with Arizona PAs
cholla_az = gpd.overlay(gbif_gdf, pa_shp, how='intersection')

# Determine how many occurrences per site
value_counts = cholla_az['Loc_Nm'].value_counts()

# Subset and drop unncessary columns
orpi_gdf = pa_shp[pa_shp['Loc_Nm'] == 'ORPI']
orpi_gdf = orpi_gdf[['Loc_Nm', 'geometry']]
tmp_gdf = pa_shp[pa_shp['Loc_Nm'] == 'Tucson Mountain Park']
tmp_gdf = tmp_gdf[['Loc_Nm', 'geometry']]

# Plot the sites
orpi_gdf_merc = orpi_gdf.to_crs(epsg=3857)
ax = orpi_gdf_merc.plot(figsize=(10, 10), edgecolor='black', alpha=0.7)
ctx.add_basemap(ax, crs=orpi_gdf_merc.crs.to_string(),
                source=ctx.providers.Esri.WorldImagery)
plt.title("Organ Pipe Cactus National Monument Site")
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

tmp_gdf_merc = tmp_gdf.to_crs(epsg=3857)
ax = tmp_gdf_merc.plot(figsize=(10, 10), edgecolor='black', alpha=0.7)
ctx.add_basemap(ax, crs=tmp_gdf_merc.crs.to_string(),
                 source=ctx.providers.Esri.WorldImagery)
plt.title("Tuscon Mountain Park Site")
ax.set_ylabel('Latitude')
ax.set_xlabel('Longitude')

# Download POLARIS data

def polaris_download(gdf, soil_prop):
    # Template for the URL
    soil_url_template = (
        "http://hydrology.cee.duke.edu/POLARIS/PROPERTIES/v1.0"
        "/{soil_prop}"
        "/mean"
        "/5_15"
        "/lat{min_lat}{max_lat}_lon{min_lon}{max_lon}.tif"
    )
    
    # Get the bounds from the GeoDataFrame
    bounds_min_lon, bounds_min_lat, bounds_max_lon, bounds_max_lat = gdf.total_bounds
    
    # Creat rasters list
    soil_rasters = []

    # Loop through the longitude and latitude ranges based on the bounds
    for min_lon in range(floor(bounds_min_lon), ceil(bounds_max_lon)):
        for min_lat in range(floor(bounds_min_lat), ceil(bounds_max_lat)):
            # Format the URL with the given soil property and bounding box values
            soil_url = soil_url_template.format(
                soil_prop=soil_prop, 
                min_lat=min_lat, max_lat=min_lat + 1, 
                min_lon=min_lon, max_lon=min_lon + 1
            )
            
            # Open the raster using rioxarray
            soil_da = rxr.open_rasterio(soil_url, mask_and_scale=True).squeeze()
            soil_rasters.append(soil_da)

    # Merge rasters
    soil_raster = rxrm.merge_arrays(soil_rasters)    
    
    # Clip rasters
    soil_raster = soil_raster.rio.clip(gdf.geometry.values)
    
    return soil_raster

orpi_clay = polaris_download(orpi_gdf, "clay")
orpi_clay.name = 'clay'
orpi_ph = polaris_download(orpi_gdf, "ph")
orpi_ph.name = 'ph'
tmp_clay = polaris_download(tmp_gdf, "clay")
tmp_clay.name = 'clay'
tmp_ph = polaris_download(tmp_gdf, "ph")
tmp_ph.name = 'ph'

# Plot the soil data
# Create a 2x2 grid for the subplots
fig, axs = plt.subplots(2, 2, figsize=(15, 12))

# Plot pH levels (orpi_ph)
pc1 = orpi_ph.plot(ax=axs[0, 0], cmap='coolwarm')
axs[0, 0].set_title('pH Levels — ORPI ')
axs[0, 0].set_xlabel('Longitude')
axs[0, 0].set_ylabel('Latitude')

# Plot clay content (orpi_clay)
pc2 = orpi_clay.plot(ax=axs[0, 1], cmap='Blues')
axs[0, 1].set_title('Clay Content — ORPI')
axs[0, 1].set_xlabel('Longitude')
axs[0, 1].set_ylabel('Latitude')

# Plot pH levels (tmp_ph)
pc3 = tmp_ph.plot(ax=axs[1, 0], cmap='coolwarm')
axs[1, 0].set_title('pH Levels — TMP')
axs[1, 0].set_xlabel('Longitude')
axs[1, 0].set_ylabel('Latitude')

# Plot clay content (tmp_clay)
pc4 = tmp_clay.plot(ax=axs[1, 1], cmap='Blues')
axs[1, 1].set_title('Clay Content — TMP')
axs[1, 1].set_xlabel('Longitude')
axs[1, 1].set_ylabel('Latitude')

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plot
plt.show()

# Set up earthaccess conncection
earthaccess.login(strategy="interactive", persist=True)


def srtm_download(directory_name, gdf):

    files_list = []
    directory = os.path.join(data_dir, directory_name)
    if not os.path.exists(directory):
        # Search and downlaod earthaccess data
        results = earthaccess.search_data(
            doi="10.5067/MEaSUREs/SRTM/SRTMGL1.003",
            bounding_box=tuple(gdf.total_bounds),
            )
        files = earthaccess.download(results, directory)
        
        # Unzip then delete .zip files
        for file in files:
            with zipfile.ZipFile(file, 'r') as zip_ref:
                zip_ref.extractall(directory)

            # After extracting, add the unzipped files to the list
            for extracted_file in zip_ref.namelist():
                extracted_file_path = os.path.join(directory, extracted_file)
                files_list.append(extracted_file_path)

            if file.endswith('.zip'):
                os.remove(file)
    else:
        # If the data is already downloaded, populate the file list directly
        files_list = [
            os.path.join(directory, f) for f in os.listdir(directory)
            if f.endswith('.hgt')
            ]

    # Open rasters
    srtm_rasters = [
        rxr.open_rasterio(file, mask_and_scale=True).squeeze()
        for file in files_list
    ]

    # Merge the rasters
    srtm_raster = rxrm.merge_arrays(srtm_rasters)
    
    # Clip the rasters
    srtm_raster = srtm_raster.rio.clip(gdf.geometry.values)

    # Reproject to UTM
    srtm_raster_utm = srtm_raster.rio.reproject(32612)

    # Calculate slope
    srtm_slope = xrs.slope(srtm_raster_utm)

    return srtm_raster, srtm_raster_utm, srtm_slope

orpi_srtm, orpi_srtm_utm, orpi_slope = srtm_download("orpi_srtm", orpi_gdf)
tmp_srtm, tmp_srtm_utm, tmp_slope = srtm_download("tmp_srtm", tmp_gdf)

# Plot elevation
# Create a figure with two subplots side by side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))

# Plot the orpi_srtm data
pc1 = orpi_srtm.plot(ax=ax1, cmap='viridis')
ax1.set_title('ORPI SRTM Elevation')
ax1.set_xlabel('Longitude')
ax1.set_ylabel('Latitude')

# Plot the tmp_srtm data
pc2 = tmp_srtm.plot(ax=ax2, cmap='viridis')
ax2.set_title('TMP SRTM Elevation')
ax2.set_xlabel('Longitude')
ax2.set_ylabel('Latitude')

# Adjust layout to prevent overlap
plt.tight_layout()

# Display the plots
plt.show()

# Create a figure with two subplots side by side
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))

# Plot the orpi_slope data
pc1 = orpi_slope.plot(ax=ax1, cmap='viridis')
ax1.set_title('ORPI Slope')
ax1.set_xlabel('Meters')
ax1.set_ylabel('Meters')

# Plot the tmp_slope data
pc2 = tmp_slope.plot(ax=ax2, cmap='viridis')
ax2.set_title('TMP Slope')
ax2.set_xlabel('Meters')
ax2.set_ylabel('Meters')

# Adjust layout to prevent overlap
plt.tight_layout()

# Display the plots
plt.show()

def convert_longitude(longitude):
    """
    Convert longitude range from 0-360 degrees to -180 to 180 degrees.

    Longitude values in some datasets are given in the range 0-360 degrees. 
    This function converts those to the range -180 to 180 degrees, which is 
    often the standard for representing longitudes.

    Args:
        longitude (float): The longitude value to convert (in the range 0-360).

    Returns:
        float: The converted longitude in the range -180 to 180.
    """
    return (longitude - 360) if longitude > 180 else longitude


def download_maca(site_name, site_gdf, emissions_scenario, gcms, year_min, year_max):
    """
    Download and process projected climate data for a given site, emissions scenario,
    and Global Climate Models (GCMs) for a specified time range.

    This function constructs URLs to access climate projection data from the MACA 
    (Multivariate Adaptive Constructed Analogs) dataset, clips the data to the 
    specified site boundary, and returns a DataArray with climate variables.

    Args:
        site_name (str): The name of the site for which climate data is requested.
        site_gdf (GeoDataFrame): A GeoDataFrame that contains the geographical boundary of the site.
        emissions_scenario (str): The emissions scenario to use (e.g., 'RCP8.5').
        gcms (list of str): A list of GCM names (e.g., ['IPSL-CM5A-MR', 'CanESM2']).
        year_min (int): The start year of the period for which data is requested.
        year_max (int): The end year of the period for which data is requested.

    Returns:
        tuple: A tuple containing four DataArrays, one for each emissions scenario
               and model combination: (hd, hw, cd, cw).
               Each DataArray contains the climate data for the given site and time range.
    """
    results = []

    # Loop through each GCM (Global Climate Model)
    for gcm in gcms:
        maca_das = []

        # Loop through each 5-year period from year_min to year_max
        for start_year in range(year_min, year_max, 5):
            end_year = start_year + 4

            # Construct the URL to access the climate data from MACA
            maca_url = (
                f'http://thredds.northwestknowledge.net:8080/thredds/dodsC/'
                f'MACAV2/{gcm}/macav2metdata_tasmax_{gcm}_r1i1p1_'
                f'{emissions_scenario}_{start_year}_{end_year}_CONUS_monthly.nc'
            )

            # Open the dataset from the constructed URL using xarray
            maca_da = xr.open_dataset(maca_url, mask_and_scale=True).squeeze()

            # Clip the dataset to the site boundary
            bounds = site_gdf.to_crs(maca_da.rio.crs).total_bounds

            # Update coordinate range with converted longitudes
            maca_da = maca_da.assign_coords(
                lon=("lon", [convert_longitude(longitude_value) for longitude_value in maca_da.lon.values])
            )
            maca_da = maca_da.rio.set_spatial_dims(x_dim='lon', y_dim='lat')
            maca_da = maca_da.rio.clip_box(*bounds)

            # Extract only the relevant variable from the dataset
            maca_da = maca_da['air_temperature']
            maca_da = maca_da.rio.write_crs(4326)

            # Store the DataArray in the list
            maca_das.append(maca_da)

        # Concatenate all the DataArrays along the 'time' dimension and calculate the maximum over time
        maca_da = xr.concat(maca_das, dim='time').max('time')

        results.append(maca_da)

    # Assign each resulting dataset to a variable for readability
    hd = results[0]
    hw = results[1]
    cd = results[2]
    cw = results[3]

    return hd, hw, cd, cw


models = [
    'IPSL-CM5A-MR', 'CanESM2', 'inmcm4', 'MRI-CGCM3'
]

# Download historical and future climate data for two different sites ('orpi' and 'tmp')
orpi_hist_hd, orpi_hist_hw, orpi_hist_cd, orpi_hist_cw = download_maca('orpi', orpi_gdf, 'historical', models, 1950, 1979)
orpi_85_hd, orpi_85_hw, orpi_85_cd, orpi_85_cw = download_maca('orpi', orpi_gdf, 'rcp85', models, 2064, 2096)
tmp_hist_hd, tmp_hist_hw, tmp_hist_cd, tmp_hist_cw = download_maca('tmp', tmp_gdf, 'historical', models, 1950, 1979)
tmp_85_hd, tmp_85_hw, tmp_85_cd, tmp_85_cw = download_maca('tmp', tmp_gdf, 'rcp85', models, 2064, 2096)


def plot_4(p1, title1, p2, title2, p3, title3, p4, title4):
    """
    Plot four datasets in a 2x2 grid of subplots.

    This function takes four datasets (DataArrays) and plots them in a 2x2 grid 
    of subplots, with each dataset having its own title and axis labels.

    Args:
        p1, p2, p3, p4 (DataArray): The datasets to plot.
        title1, title2, title3, title4 (str): The titles for each subplot.

    Returns:
        None
    """
    # Create a 2x2 grid for the subplots
    fig, axs = plt.subplots(2, 2, figsize=(15, 12))

    # Plot 1
    p1.plot(ax=axs[0, 0], cmap='viridis')
    axs[0, 0].set_title(title1)
    axs[0, 0].set_xlabel('Longitude')
    axs[0, 0].set_ylabel('Latitude')

    # Plot 2
    p2.plot(ax=axs[0, 1], cmap='viridis')
    axs[0, 1].set_title(title2)
    axs[0, 1].set_xlabel('Longitude')
    axs[0, 1].set_ylabel('Latitude')

    # Plot 3
    p3.plot(ax=axs[1, 0], cmap='viridis')
    axs[1, 0].set_title(title3)
    axs[1, 0].set_xlabel('Longitude')
    axs[1, 0].set_ylabel('Latitude')

    # Plot 4
    p4.plot(ax=axs[1, 1], cmap='viridis')
    axs[1, 1].set_title(title4)
    axs[1, 1].set_xlabel('Longitude')
    axs[1, 1].set_ylabel('Latitude')

    # Adjust layout to prevent overlap of subplots
    plt.tight_layout()

    # Display the plot
    plt.show()

    return

def harmonize_rasters(rasters, reference_raster):
    """
    Harmonizes a list of rasters to match the coordinates of a reference dataset.

    Parameters
    ----------
    rasters : list of xarray.DataArray or xarray.Dataset
        A list of xarray datasets or data arrays to be harmonized.
        Each dataset will be reprojected to match the coordinates of xds_match.

    xds_match : xarray.DataArray or xarray.Dataset
        The reference dataset to harmonize all rasters to.

    Returns
    -------
    list of xarray.DataArray or xarray.Dataset
        A list of harmonized xarray datasets or data arrays.
    """
    harmonized_rasters = []
    for raster in rasters:
        
        xds_repr_match = raster.rio.reproject_match(reference_raster)
        xds_repr_match = xds_repr_match.assign_coords({
            "x": reference_raster.x,
            "y": reference_raster.y,
        })
        harmonized_rasters.append(xds_repr_match)

    return harmonized_rasters

orpi_rasters = [orpi_clay, orpi_ph, orpi_slope, orpi_hist_hd, orpi_hist_hw, orpi_hist_cd, orpi_hist_cw, orpi_85_hd, orpi_85_hw, orpi_85_cd, orpi_85_cw]
harm_orpi_rasters = harmonize_rasters(orpi_rasters, orpi_slope)
tmp_rasters = [tmp_clay, tmp_ph, tmp_slope, tmp_hist_hd, tmp_hist_hw, tmp_hist_cd, tmp_hist_cw, tmp_85_hd, tmp_85_hw, tmp_85_cd, tmp_85_cw]
harm_tmp_rasters = harmonize_rasters(tmp_rasters, tmp_slope)

def fuzzy_logic_rasters(rasters, variable_names):
    """
    Applies fuzzy logic to a list of raster data based on their variable names.

    Parameters:
    - rasters: list of 2D array-like (e.g., xarray.DataArray, numpy array)
      List of raster data to which fuzzy logic will be applied.
    - variable_names: list of str
      List of variable names (e.g., 'precip', 'sand', 'ph', etc.) to match with rasters.
    Returns:
    - fuzzy_rasters: list of 2D array-like
      List of transformed rasters after applying the fuzzy logic function.
    """
    
    # Define breakpoints for different variables
    breakpoints = {
        'clay': [-1, 0, 35, 55],
        'ph': [5, 6, 7, 7.5],
        'slope': [-1, 0, 30, 50], 
        'air_temperature': [268, 288 , 311, 316],
        # 'precipitation': [25, 100, 425, 500], 
    }
    
    fuzzy_rasters = []

    # Process each raster and match to variable name based on raster name
    for raster in rasters:
        # Try to match the raster name to a variable name
        variable_name = None
        for name in variable_names:
            if name.lower() in raster.name.lower():  # Check if variable name is in the raster name
                variable_name = name
                break
        
        if variable_name is None:
            raise ValueError(f"No matching variable name found for raster: {raster.name}")

        # Ensure the variable name exists in the breakpoints dictionary
        if variable_name not in breakpoints:
            raise ValueError(f"Invalid variable name: {variable_name}")
        
        # Get the breakpoints for the specified variable
        params = breakpoints[variable_name]
        
        # Get the shape of the raster data
        shape = raster.values.shape

        # Make a copy of the raster data
        raster_fuzz = raster.copy()

        # Apply the trapezoidal fuzzy membership function
        raster_fuzz.values = np.reshape(
            skfuzzy.trapmf(
                raster.values.flatten(),  # Flatten the raster data for processing
                params  # Trapezoidal membership function breakpoints
            ),
            shape  # Reshape the resulting membership values to match the original shape
        )
        
        # Append the fuzzified raster to the list
        fuzzy_rasters.append(raster_fuzz)
        
    return fuzzy_rasters

variable_names = ['clay', 'ph', 'slope', 'air_temperature']
# Apply fuzzy logic to rasters
orpi_fuzz = fuzzy_logic_rasters(harm_orpi_rasters, variable_names)
tmp_fuzz = fuzzy_logic_rasters(harm_tmp_rasters, variable_names)

# Combine rasters
orpi_hist_hd = orpi_fuzz[0]*orpi_fuzz[1]*orpi_fuzz[2]*orpi_fuzz[3]
orpi_hist_hw = orpi_fuzz[0]*orpi_fuzz[1]*orpi_fuzz[2]*orpi_fuzz[4]
orpi_hist_cd = orpi_fuzz[0]*orpi_fuzz[1]*orpi_fuzz[2]*orpi_fuzz[5]
orpi_hist_cw = orpi_fuzz[0]*orpi_fuzz[1]*orpi_fuzz[2]*orpi_fuzz[6]

orpi_85_hd = orpi_fuzz[0]*orpi_fuzz[1]*orpi_fuzz[2]*orpi_fuzz[7]
orpi_85_hw = orpi_fuzz[0]*orpi_fuzz[1]*orpi_fuzz[2]*orpi_fuzz[8]
orpi_85_cd = orpi_fuzz[0]*orpi_fuzz[1]*orpi_fuzz[2]*orpi_fuzz[9]
orpi_85_cw = orpi_fuzz[0]*orpi_fuzz[1]*orpi_fuzz[2]*orpi_fuzz[10]

tmp_hist_hd = tmp_fuzz[0] * tmp_fuzz[1] * tmp_fuzz[2] * tmp_fuzz[3]
tmp_hist_hw = tmp_fuzz[0] * tmp_fuzz[1] * tmp_fuzz[2] * tmp_fuzz[4]
tmp_hist_cd = tmp_fuzz[0] * tmp_fuzz[1] * tmp_fuzz[2] * tmp_fuzz[5]
tmp_hist_cw = tmp_fuzz[0] * tmp_fuzz[1] * tmp_fuzz[2] * tmp_fuzz[6]

tmp_85_hd = tmp_fuzz[0] * tmp_fuzz[1] * tmp_fuzz[2] * tmp_fuzz[7]
tmp_85_hw = tmp_fuzz[0] * tmp_fuzz[1] * tmp_fuzz[2] * tmp_fuzz[8] 
tmp_85_cd = tmp_fuzz[0] * tmp_fuzz[1] * tmp_fuzz[2] * tmp_fuzz[9]
tmp_85_cw = tmp_fuzz[0] * tmp_fuzz[1] * tmp_fuzz[2] * tmp_fuzz[10]

# Plot the habitat suitability 
plot_4(
     orpi_hist_hd, 'ORPI IPSL-CM5A-MR 1950-1979', orpi_hist_hw, 'ORPI CanESM2 1950-1979',
     orpi_hist_cd, 'ORPI inmcm4 1950-1979', orpi_hist_cw, 'ORPI MRI-CGCM3 1950-1979'
     )
plot_4(
     orpi_85_hd, 'ORPI IPSL-CM5A-MR 2064-2096', orpi_85_hw, 'ORPI CanESM2 2064-2096',
     orpi_85_cd, 'ORPI inmcm4 2064-2096', orpi_85_cw, 'ORPI MRI-CGCM3 2064-2096'
     )
plot_4(
    tmp_hist_hd, 'TMP IPSL-CM5A-MR 1950-1979', tmp_hist_hw, 'TMP CanESM2 1950-1979',
    tmp_hist_cd, 'TMP inmcm4 1950-1979', tmp_hist_cw, 'TMP MRI-CGCM3 1950-1979'
    )
plot_4(
    tmp_85_hd, 'TMP IPSL-CM5A-MR 2064-2096', tmp_85_hw, 'TMP CanESM2 2064-2096',
    tmp_85_cd, 'TMP inmcm4 2064-2096', tmp_85_cw, 'TMP MRI-CGCM3 2064-2096'
    )
