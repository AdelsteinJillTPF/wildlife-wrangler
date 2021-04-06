import pandas as pd
import geopandas as gpd
import shapely
import sqlite3
import matplotlib.pyplot as plt

conn = sqlite3.connect("T:/Occurrence_Records/withRTest3.sqlite")
cursor = conn.cursor()

records = (pd.read_sql("""SELECT decimalLatitude, decimalLongitude FROM occurrence_records""", con=conn)
               .astype({'decimalLongitude': 'float', 'decimalLatitude': 'float'}))

print(records.head())

# Make data frame spatial
gdf = gpd.GeoDataFrame(records, geometry=gpd.points_from_xy(records['decimalLongitude'],
                                                   records['decimalLatitude']))
gdf.crs="EPSG:4326"

world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
base = world.plot()
gdf.plot(ax=base, marker='o', color='k', markersize=5)
plt.show()

# Plot within USA
usa_bbox = np.array([-124.725839,   24.498131,  -66.949895,   49.384358])
fig, ax = plt.subplots()
ax.set_xlim(([usa_bbox[0],  usa_bbox[2]]))
ax.set_ylim(([usa_bbox[1],  usa_bbox[3]]))
world.plot()
gdf.plot(ax=ax, marker='o', color='k', markersize=5)
plt.show()
