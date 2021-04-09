import pandas as pd
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
                      .astype({'decimalLongitude': 'float', 'decimalLatitude': 'float'}))

# Make the data frame spatial
gdf = gpd.GeoDataFrame(record_coordinates, geometry=gpd.points_from_xy(records['decimalLongitude'],
                                                   records['decimalLatitude']))
gdf.crs="EPSG:4326"

# Create world map
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
base = world.plot()
gdf.plot(ax=base, marker='o', color='k', markersize=5)
plt.show()

# Get bounding box for coordinates for use later )))))))))))))))))))))))))))))))))))))))))))))
coordinate_bbox = gdf.geometry.total_bounds

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

# Buffer the coordinates with radius_m to get a new geometry column
# reproject
# new footprint column for the polygon geometry
gdf2["footprint"] = gdf2.apply(lambda x: x.geometry.buffer(x.radius_m), axis=1)

# Get bbox of coordinates
coordinate_bbox = gdf.geometry.total_bounds

# Map the buffered points/footprints
fig2, ax2 = plt.figure(figsize=(10,6))
gdf.set_geometry("footprint")
