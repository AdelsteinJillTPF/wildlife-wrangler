import os
import pandas as pd
import numpy as np
import geopandas as gpd
import shapely
import sqlite3
import matplotlib.pyplot as plt
output_db_conn = sqlite3.connect("T:/Occurrence_Records/withRTest3.sqlite")
cursor = output_db_conn.cursor()
#############################

# Get the record coordinates as a data frame
record_coordinates = (pd.read_sql("""SELECT decimalLatitude, decimalLongitude, radius_meters
                                     FROM occurrence_records""",
                                  con=output_db_conn)
                      .astype({'decimalLongitude': 'float', 'decimalLatitude': 'float',
                               'radius_meters': 'float'}))

# Make the data frame spatial
gdf = gpd.GeoDataFrame(record_coordinates, geometry=gpd.points_from_xy(record_coordinates['decimalLongitude'],
                                                   record_coordinates['decimalLatitude']))
gdf.crs={'init' :'epsg:4326'}

# Create world map
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
base = world.plot()
gdf.plot(ax=base, marker='o', color='k', markersize=5)
plt.show()

# Plot within USA
usa_bbox = np.array([-124.725839,   24.498131,  -66.949895,   49.384358])
fig, ax = plt.subplots()
ax.set_xlim(([usa_bbox[0],  usa_bbox[2]]))
ax.set_ylim(([usa_bbox[1],  usa_bbox[3]]))
world.plot(ax=ax)
gdf.plot(ax=ax, marker='o', color='k', markersize=5)
plt.show()

# Zoomed in plot <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Get a states basemap
states = gpd.read_file(os.getcwd() + '/data/us_states.shp')

# Buffer the coordinates with radius_meters to get a new geometry column
# reproject states and record coordinates
states = states.to_crs(epsg=5070)
gdf_5070 = gdf.to_crs(epsg=5070)

# new footprint column for the polygon geometry
footprints=gdf_5070.copy()
footprints['footprint']=footprints.apply(lambda x: x.geometry.buffer(x.radius_meters), axis=1)
footprints.set_geometry(col='footprint', inplace=True, drop=True)

# Get bbox of coordinates
coordinate_bbox = footprints.geometry.total_bounds

# Map the buffered points/footprints
fig, ax = plt.subplots()
ax.set_xlim(([coordinate_bbox[0],  coordinate_bbox[2]]))
ax.set_ylim(([coordinate_bbox[1],  coordinate_bbox[3]]))
states.plot(ax=ax)
footprints.boundary.plot(ax=ax, color='k')
plt.show()
