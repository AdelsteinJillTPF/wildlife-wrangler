def MapShapefilePolygons(map_these, title):
    """
    Displays shapefiles on a simple CONUS basemap.  Maps are plotted in the order
    provided so put the top map last in the list.  You can specify a column
    to map as well as custom colors for it.  This function may not be very robust
    to other applications.

    NOTE: The shapefiles have to be in WGS84 CRS.

    (dict, str) -> displays maps, returns matplotlib.pyplot figure

    Arguments:
    map_these -- list of dictionaries for shapefiles you want to display in
                CONUS. Each dictionary should have the following format, but
                some are unnecessary if 'column' doesn't = 'None'.  The critical
                ones are file, column, and drawbounds.  Column_colors is needed
                if column isn't 'None'.  Others are needed if it is 'None'.
                    {'file': '/path/to/your/shapfile',
                     'alias': 'my layer'
                     'column': None,
                     'column_colors': {0: 'k', 1: 'r'}
                    'linecolor': 'k',
                    'fillcolor': 'k',
                    'linewidth': 1,
                    'drawbounds': True
                    'marker': 's'}
    title -- title for the map.
    """
    # Packages needed for plotting
    import warnings
    warnings.filterwarnings('ignore')
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import PathPatch

    # Basemap
    fig = plt.figure(figsize=(15,12))
    ax = plt.subplot(1,1,1)
    map = Basemap(projection='aea', resolution='l', lon_0=-95.5, lat_0=39.0,
                  height=3200000, width=5000000)
    map.drawcoastlines(color='grey')
    map.drawstates(color='grey')
    map.drawcountries(color='grey')
    map.fillcontinents(color='#a2d0a2',lake_color='#a9cfdc')
    map.drawmapboundary(fill_color='#a9cfdc')

    for mapfile in map_these:
        if mapfile['column'] == None:
            # Add shapefiles to the map
            if mapfile['fillcolor'] == None:
                map.readshapefile(mapfile['file'], 'mapfile',
                                  drawbounds=mapfile['drawbounds'],
                                  linewidth=mapfile['linewidth'],
                                  color=mapfile['linecolor'])
                # Empty scatter plot for the legend
                plt.scatter([], [], c='', edgecolor=mapfile['linecolor'],
                            alpha=1, label=mapfile['alias'], s=100,
                            marker=mapfile['marker'])

            else:
                map.readshapefile(mapfile['file'], 'mapfile',
                          drawbounds=mapfile['drawbounds'])
                # Code for extra formatting -- filling in polygons setting border
                # color
                patches = []
                for info, shape in zip(map.mapfile_info, map.mapfile):
                    patches.append(Polygon(np.array(shape), True))
                ax.add_collection(PatchCollection(patches,
                                                  facecolor= mapfile['fillcolor'],
                                                  edgecolor=mapfile['linecolor'],
                                                  linewidths=mapfile['linewidth'],
                                                  zorder=2))
                # Empty scatter plot for the legend
                plt.scatter([], [], c=mapfile['fillcolor'],
                            edgecolors=mapfile['linecolor'],
                            alpha=1, label=mapfile['alias'], s=100,
                            marker=mapfile['marker'])

        else:
            map.readshapefile(mapfile['file'], 'mapfile', drawbounds=mapfile['drawbounds'])
            for info, shape in zip(map.mapfile_info, map.mapfile):
                for thang in mapfile['column_colors'].keys():
                    if info[mapfile['column']] == thang:
                        x, y = zip(*shape)
                        map.plot(x, y, marker=None, color=mapfile['column_colors'][thang])

            # Empty scatter plot for the legend
            for seal in mapfile['column_colors'].keys():
                plt.scatter([], [], c=mapfile['column_colors'][seal],
                            edgecolors=mapfile['column_colors'][seal],
                            alpha=1, label=mapfile['value_alias'][seal],
                            s=100, marker=mapfile['marker'])

    # Legend -- the method that works is ridiculous but necessary; you have
    #           to add empty scatter plots with the symbology you want for
    #           each shapefile legend entry and then call the legend.  See
    #           plt.scatter(...) lines above.
    plt.legend(scatterpoints=1, frameon=True, labelspacing=1, loc='lower left',
               framealpha=1, fontsize='x-large')

    # Title
    plt.title(title, fontsize=20, pad=-40, backgroundcolor='w')
    return

def getGBIFcode(name, rank='species'):
    """
    Returns the GBIF species code for a scientific name.

    Example: gbifcode = getGBIFcode(name = "Dendroica cerulea")
    """
    from pygbif import species
    key = species.name_backbone(name = name, rank='species')['usageKey']
    return key

def getRecordDetails(key):
    """
    Returns a dictionary holding all GBIF details about the record.

    Example: details = getRecordDetails(key = 1265907957)
    """
    from pygbif import occurrences
    details = occurrences.get(key = key)
    return details

def drop_duplicates_latlongdate(df):
    '''
    Function to find and remove duplicate occurrence records within the
    wildlife wrangler workflow.  When duplicates exist, the record with the
    higher decimal precision is kept, and if precision values are equal, then the
    record with the higher individualCount is retained. Accounts for existence
    of records with a mix of decimal precision in latitude and longitude
    values. The process is a little complex.   The first data frame is cleaned
    up by dropping duplicates based on which record has greater individual count.
    Before doing that, records with unequal decimal precision in the latitude
    and longitude fields and those fields are truncated to the shorter precision
    present. An input data frame likely contains records with equal decimal
    precision in latitude and longitude fields, but that is lower than the rest
    (i.e. latitude and longitude have 3 places right of the decimal whereas most
    records have 4).  Duplication may occur between lower and higher precision
    records at the lower precision.  Therefore, duplication must be assessed at
    each of the lower precision levels present. The strategy for that is to, at
    each precision level, split the main data frame in two: one with records
    having the precision level of the investigation and another with records
    greater than the precision level. The "greater than" data frame records'
    latitude and longitude values are then truncated to the precision level.
    Records are identified from the "equals precision" data frame that have
    their latitude, longitude, and date values represented in the "greater than"
    df, and such records ID’s are collected in a list of records to be removed
    from the input/main data frame.  This process is iterated over all precision
    levels present in the data.

    Parameters
    ----------
    df : Input pandas dataframe.

    Returns
    -------
    df2 : A dataframe equal to df but without duplicates.  Use to drop records
    from the occurrences table.

    '''
    from datetime import datetime
    import pandas as pd
    startduptime = datetime.now()

    # Record df length before removing duplicates
    initial_length = len(df)

    """
    ############ RECTIFY UNEQUAL LAT-LONG PRECISION
    First, trim decimal length in cases where decimal length differs between
    latitude and longitude values, result is equal latitude and longitude
    length.  Record the trimmed decimal precision in a temp column for use later
    as a record to "native" precision.
    """
    df['dup_latPlaces'] = [len(x.split(".")[1]) for x in df['latitude']]
    df['dup_lonPlaces'] = [len(x.split(".")[1]) for x in df['longitude']]
    df['dup_OGprec'] = df['dup_latPlaces']
    df22 = df[df['dup_latPlaces'] != df['dup_lonPlaces']]
    for i in df22.index:
        x = df22.loc[i]
        if x['dup_latPlaces'] < x['dup_lonPlaces']:
            trim_len = int(x['dup_latPlaces'])
        else:
            trim_len = int(x['dup_lonPlaces'])
        df.loc[i, 'latitude'] = x['latitude'][:trim_len + 3]
        df.loc[i, 'longitude'] = x['longitude'][:trim_len + 4]
        # Record the resulting precision for reference later
        df.loc[i, 'dup_OGprec'] = trim_len
    df.drop(['dup_latPlaces', 'dup_lonPlaces'], axis=1, inplace=True)

    """
    ########  INITIAL DROP OF DUPLICATES
    Initial drop of duplicates on 'latitude', 'longitude', 'occurrenceDate',
    keeping the first (highest individual count)
    Sort so that the highest individual count is first ## ADD OCCURRENCEDATE BACK IN
    """
    df.sort_values(by=['latitude', 'longitude', 'occurrenceDate',
                        'individualCount'],
                    ascending=False, inplace=True, kind='mergesort',
                    na_position='last')

    df.drop_duplicates(subset=['latitude', 'longitude', 'occurrenceDate'],
                       keep='first', inplace=True)

    """
    #########  FIND IMPRECISE DUPLICATES
    Get a list of "native" precisions that are present in the data to loop through.
    Next, iterate through this list collecting id's of records that need to be
    removed from the main df.
    """
    # Get list of unique precisions.  Order is important: descending.
    precisions = list(set(df['dup_OGprec']))
    precisions.sort(reverse=True)
    # The highest precisions listed at this point has already been done: drop it.
    precisions = precisions[1:]

    # List for collecting records that are duplicates
    duplis = []

    # The precision-specific duplicate testing happens repeatedly, so make it a
    # function.
    def drop_duplicates(precision, df):
        """
        Function to find undesirable duplicates at a particular decimal precision.

        Parameters
        ----------
        precision : The level of precision (places right of decimal) in latitude
        and longitude values for the assessment of duplicates.
        df : dataframe to assess and drop duplicates from.  This function works
              'inplace'.
        """
        # Create a df with records from the input df having decimal precision > the
        # precision level being assessed.
        dfLonger = df[df['dup_OGprec'] > precision].copy()
        # Truncate lat and long values
        dfLonger['latitude'] = [x[:precision + 3] for x in dfLonger['latitude']]
        dfLonger['longitude'] = [x[:precision + 4] for x in dfLonger['longitude']]

        # Create a df with records having the precision being
        # investigated
        dfShorter1 = df[df['dup_OGprec'] == precision]

        # Find records in dfShorter1 with latitude, longitude, date combo
        # existing in dfLonger and append to list of duplis
        dfduplis = pd.merge(dfShorter1, dfLonger, how='inner',
                            on=['latitude', 'longitude', 'occurrenceDate'])
        dups_ids = dfduplis['occ_id_x']
        for d in dups_ids:
            duplis.append(d)

    # Drop latitude longitude duplicates at lower decimal precisions
    for p in precisions:
        drop_duplicates(p, df)

    # Drop rows from the current main df that have been identified as duplicates.
    df2 = df[df['occ_id'].isin(duplis) == False].copy()

    # Drop excess columns
    df2.drop(['dup_OGprec'], inplace=True, axis=1)

    duptime = datetime.now() - startduptime
    print(str(initial_length - len(df2)) + " duplicate records dropped: {0}".format(duptime))
    return df2

def exportSHP(database, table, column, outFile):
    '''
    Exports a spatialite geometry column as a shapefile.

    Parameters:
    database -- the sqlite database to use.  Must have spatial data.
    table -- name of the table with geometry in it.
    column -- column name of the geometry to export as a shapefile.
    outFile -- Path (and name) of the file to be created.
    '''
    from datetime import datetime
    import sqlite3
    import platform
    import os
    # Environment variables need to be handled
    if platform.system() == 'Windows':
        os.environ['PATH'] = os.environ['PATH'] + ';' + 'C:/Spatialite'
        os.environ['SPATIALITE_SECURITY'] = 'relaxed'
    if platform.system() == 'Darwin':
        os.environ['SPATIALITE_SECURITY'] = 'relaxed'
    exporttime1 = datetime.now()
    conn = sqlite3.connect(database, isolation_level='DEFERRED')
    conn.enable_load_extension(True)
    conn.execute('SELECT load_extension("mod_spatialite")')
    cursor = conn.cursor()
    cursor.execute("""SELECT ExportSHP('{0}', '{1}', '{2}',
                    'utf-8');""".format(table, column, outFile))
    conn.commit()
    conn.close()
    print("Exported shapefile: " + str(datetime.now() - exporttime1))

def ccw_wkt_from_shp(shapefile, out_txt):
    """
    Creates wkt with coordinates oriented counter clockwise for a given shapefile.
    Shapefiles are oriented clockwise, which is incompatible with spatial queries
    in many database management systems.  Use this to generate wkt that you can
    copy and paste into queries.

    (str, str) = text written to shpefile

    Arguments:
    shapefile -- path to the shpefile to read.
    out_txt -- path to the text file to write the wkt to.
    """

    import fiona
    import shapely
    from shapely.geometry import shape, Polygon, LinearRing
    #from shapely.wkb import dumps, loads

    # Read in a shapefile of polygon of interest.  It must be in CRS 4326
    # First get a fiona collection
    c = fiona.open(shapefile, 'r')

    if c.crs['init'] == 'epsg:4326':
        # Next make it a shapely polygon object
        poly = shape(c[0]['geometry'])

        # Use LinearRing to determine if coordinates are listed clockwise
        coords = c[0]["geometry"]["coordinates"][0]
        lr = LinearRing(coords)
        if lr.is_ccw == False:
            # Reverse coordinates to make them counter clockwise
            print("Points were clockwise......reversing")
            #coords.reverse()
            # Make the polygon's outer ring counter clockwise
            poly2 = shapely.geometry.polygon.orient(poly, sign=1.0)
            # Get the well-known text version of the polygon
            wkt = poly2.wkt
        else:
            print("Points were already counter clockwise")
            # Get the well-known text version of the polygon
            wkt = poly.wkt

        # Write WKT to text file
        with open(out_txt, 'w+') as file:
            file.write(wkt)
            print("WKT written to {0}".format(out_txt))

        # close the collections
        c.close()
    else:
        print("You need to reproject the shapefile to EPSG:4326")
    return

def retrieve_gbif_occurrences(codeDir, taxon_id, paramdb, spdb,
                              gbif_req_id, gbif_filter_id, default_coordUncertainty,
                              outDir, summary_name, username, password, email):
    """
    Retrieves species occurrence records from the GBIF API.  Filters occurrence
    records, buffers the xy points, and saves them in a database.  Finally,
    exports some Shapefiles.

    Arguments:
    codeDir -- directory of this code repo.
    taxon_id -- project-specific identifier for the taxon concept.
    paramdb -- path to the parameter database.
    spdb -- occurrence record database to be created by this function.
    gbif_req_id -- GBIF request ID for the process.
    gbif_filter_id -- GBIF filter ID for the process.
    default_coordUncertainty -- distance in meters to use if no coordinate
        Uncertainty is specified for a record.
    outDir -- where to save maps that are exported by this process.
    summary_name -- a short name for some file names.
    sp_geometry -- True or False to use geometry saved with taxon concept when
        filtering records.  Request geometry is always used if provided.
        NOTE: This option is forthcoming and is currently hard-coded to "True".
    """
    import pandas as pd
    pd.set_option('display.width', 1000)
    import sqlite3
    from pygbif import occurrences
    import os
    os.chdir('/')
    import json
    import platform
    import shapely
    from shapely.wkt import dumps, loads
    from datetime import datetime
    import sys
    import shutil
    from dwca.read import DwCAReader
    import numpy as np

    sp_geometry=True
    # Environment variables need to be handled
    if platform.system() == 'Windows':
        os.environ['PATH'] = os.environ['PATH'] + ';' + 'C:/Spatialite'
        os.environ['SPATIALITE_SECURITY'] = 'relaxed'

    if platform.system() == 'Darwin':
        print("Support for Darwin operating systems is forthcoming")
        #os.putenv('SPATIALITE_SECURITY', 'relaxed')
        os.environ['SPATIALITE_SECURITY'] = 'relaxed'

    print("SPATIALITE_SECURITY set to " + os.environ['SPATIALITE_SECURITY'])


    #############################################################################
    #                              Taxon-concept
    #############################################################################
    os.chdir(codeDir)
    # Get species info from requests database
    conn2 = sqlite3.connect(paramdb, isolation_level='DEFERRED')
    cursor2 = conn2.cursor()
    sql_tax = """SELECT gbif_id, common_name, scientific_name,
                        detection_distance_meters, geometry
                 FROM taxa_concepts
                 WHERE taxon_id = '{0}';""".format(taxon_id)
    concept = cursor2.execute(sql_tax).fetchall()[0]
    gbif_id = concept[0]
    common_name = concept[1]
    scientific_name = concept[2]
    det_dist = concept[3]
    sp_geom =concept[4]


    ############################################################################
    #####################    Create Occurrence Database    #####################
    ############################################################################
    """
    Description: Create a database for storing occurrence and taxon concept
    data.  Needs to have spatial querying functionality.
    """
    makedb1 = datetime.now()
    spdb = spdb
    # Delete the database if it already exists
    if os.path.exists(spdb):
        os.remove(spdb)

    # Create or connect to the database
    conn = sqlite3.connect(spdb, isolation_level='DEFERRED')
    conn.enable_load_extension(True)
    conn.execute('SELECT load_extension("mod_spatialite")')
    cursor = conn.cursor()

    # Make database spatial
    conn.executescript('''SELECT InitSpatialMetaData(1);''')
    conn.commit()

    ################################################# Create tables
    sql_cdb = """
            /* Create a table for occurrence records, WITH GEOMETRY */
            CREATE TABLE IF NOT EXISTS occurrences (
                    occ_id INTEGER NOT NULL PRIMARY KEY,
                    taxon_id INTEGER NOT NULL,
                    basisOfRecord TEXT,
                    issues TEXT,
                    collectionCode TEXT,
                    institutionCode TEXT,
                    datasetName TEXT,
                    identificationQualifier TEXT,
                    source TEXT NOT NULL,
                    request_id TEXT NOT NULL,
                    filter_id TEXT NOT NULL,
                    latitude TEXT,
                    longitude TEXT,
                    coordinateUncertaintyInMeters INTEGER,
                    occurrenceDate TEXT,
                    retrievalDate TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
                    individualCount INTEGER DEFAULT 1,
                    dataGeneralizations TEXT,
                    remarks TEXT,
                    detection_distance INTEGER,
                    radius_meters INTEGER,
                    footprintWKT TEXT,
                    weight INTEGER DEFAULT 10,
                    weight_notes TEXT,
                    doi_search TEXT,
                        FOREIGN KEY (taxon_id) REFERENCES taxa(taxon_id)
                        ON UPDATE RESTRICT
                        ON DELETE NO ACTION);

            SELECT AddGeometryColumn('occurrences', 'geom_xy4326', 4326, 'POINT',
                                     'XY');
    """
    cursor.executescript(sql_cdb)
    makedb2 = datetime.now()
    print("Created occurrence db: " + str(makedb2 - makedb1))

    requesttime1 = datetime.now()

    ############################################################################
    ###########################   Get Filters   ################################
    ############################################################################
    """
    Retrieve filter parameters from the parameters database.
    """
    def get_filter(column, where_column, table):
        '''
        Get the values of a filter from the parameters database

        Arguments:
        column -- string name of column to select from.
        where_column -- string name of column to condition selection on.
        table -- string name of table to query
        '''
        if table == 'gbif_requests':
            sql = """ SELECT {0} FROM {2} WHERE request_id = '{1}'""".format(column,
                                                                         where_column,
                                                                         table)
        if table == 'gbif_filters':
            sql = """ SELECT {0} FROM {2} WHERE filter_id = '{1}'""".format(column,
                                                                         where_column,
                                                                         table)
        filter = cursor2.execute(sql).fetchone()[0]
        return filter

    ############################# RETRIEVE REQUEST PARAMETERS
    # Up-front filters are an opportunity to lighten the load from the start.
    latRange = get_filter('lat_range', gbif_req_id, 'gbif_requests')
    lonRange = get_filter('lon_range', gbif_req_id, 'gbif_requests')
    years = get_filter('years_range', gbif_req_id, 'gbif_requests')
    months = get_filter('months_range', gbif_req_id, 'gbif_requests')
    geoIssue = get_filter('geoissue', gbif_req_id, 'gbif_requests')
    if geoIssue == 'None':
        geoIssue = None
    coordinate = get_filter('coordinate', gbif_req_id, 'gbif_requests')
    country = get_filter('country', gbif_req_id, 'gbif_requests')
    if country == "None":
        country = None
    poly0 = get_filter('geometry', gbif_req_id, 'gbif_requests')

    ################################################ SORT OUT GEOMETRIES
    # A geometry could also be stated for the species, assess what to do
    # It could also be that user opted not to use species geometry.
    if sp_geometry == False:
        sp_geom = None
    if poly0 == None and sp_geom == None:
        poly = None
    elif poly0 != None and sp_geom == None:
        poly = poly0
    elif poly0 == None and sp_geom != None:
        poly = sp_geom
    elif poly0 != None and sp_geom != None:
        # Get/use the intersection of the two polygons
        filter_polygon = shapely.wkt.loads(poly0)
        sp_polygon = shapely.wkt.loads(sp_geom)
        poly_intersection = filter_polygon.intersection(sp_polygon)
        poly = shapely.wkt.dumps(poly_intersection)

    ###################  RETRIEVE POST-REQUEST FILTER PARAMETERS
    filt_coordUncertainty = get_filter('has_coordinate_uncertainty',
                                       gbif_filter_id, 'gbif_filters')
    filt_maxcoord = get_filter('max_coordinate_uncertainty', gbif_filter_id,
                               'gbif_filters')
    filt_collection = get_filter('collection_codes_omit', gbif_filter_id,
                                 'gbif_filters')
    if type(filt_collection) == str:
        filt_collection = list(filt_collection.split(', '))
    else:
        filt_collection = []
    filt_instit = get_filter('institutions_omit', gbif_filter_id, 'gbif_filters')
    if type(filt_instit) == str:
        filt_instit = list(filt_instit.split(', '))
    else:
        filt_instit = []
    filt_bases = get_filter('bases_omit', gbif_filter_id, 'gbif_filters')
    if type(filt_bases) == str:
        filt_bases = list(filt_bases.split(', '))
    else:
        filt_bases = []
    filt_issues = get_filter('issues_omit', gbif_filter_id, 'gbif_filters')
    if type(filt_issues) == str:
        filt_issues = list(filt_issues.split(', '))
    else:
        filt_issues = []
    filt_sampling = get_filter('sampling_protocols_omit', gbif_filter_id, 'gbif_filters')
    if type(filt_sampling) == str:
        filt_sampling = list(filt_sampling.split(', '))
    else:
        filt_sampling = []
    print("Got request params and sorted out geometry constraints: " + str(datetime.now() - requesttime1))
    requesttime2 = datetime.now()

    # List of informative df columns/dictionary keys to keep (used later)
    keeper_keys = ['basisOfRecord', 'individualCount', 'scientificName',
                   'decimalLongitude', 'decimalLatitude',
                   'coordinateUncertaintyInMeters',
                   'eventDate', 'issue', 'issues', 'gbifID', 'id',
                   'dataGeneralizations', 'eventRemarks', 'locality',
                   'locationRemarks', 'collectionCode',
                   'samplingProtocol', 'institutionCode', 'establishmentMeans',
                   'institutionID', 'footprintWKT', 'identificationQualifier',
                   'occurrenceRemarks', 'datasetName']
    keeper_keys.sort()

    ############################################################################
    #######################      Get GBIF Records      #########################
    ############################################################################
    """
    Retrieve GBIF records for a species and save appropriate attributes in the
    occurrence database.
    """
    ####################################### HOW MANY RECORDS EXIST TO PULL FROM?
    ############################################################################
    # First, find out how many records there are that meet criteria
    occ_search = occurrences.search(gbif_id,
                                    year=years,
                                    month=months,
                                    decimalLatitude=latRange,
                                    decimalLongitude=lonRange,
                                    hasGeospatialIssue=geoIssue,
                                    hasCoordinate=coordinate,
                                    country=country,
                                    geometry=poly)
    occ_count=occ_search['count']
    print(str(occ_count) + " records available")


    ############################################################################
    #                         < 100,000 RECORDS (small)
    ############################################################################
    if occ_count <100000:
        # Get occurrences in batches, saving into master list
        alloccs = []
        batches = range(0, occ_count, 300)
        for i in batches:
            occ_json = occurrences.search(gbif_id,
                                          limit=300,
                                          offset=i,
                                          year=years,
                                          month=months,
                                          decimalLatitude=latRange,
                                          decimalLongitude=lonRange,
                                          hasGeospatialIssue=geoIssue,
                                          hasCoordinate=coordinate,
                                          country=country,
                                          geometry=poly)
            occs = occ_json['results']
            alloccs = alloccs + occs

        print("Downloaded records: " + str(datetime.now() - requesttime2))

        ######################################  LOAD JSON RECORDS INTO DATAFRAME
        ########################################################################
        dfRaw = pd.DataFrame(columns=keeper_keys)
        insertDict = {}
        for x in keeper_keys:
            insertDict[x] = []
        for x in alloccs:
            present_keys = list(set(x.keys()) & set(keeper_keys))
            for y in present_keys:
                insertDict[y] = insertDict[y] + [str(x[y])]
            missing_keys = list(set(keeper_keys) - set(x.keys()))
            for z in missing_keys:
                insertDict[z] = insertDict[z] + ["UNKNOWN"]
        insertDF = pd.DataFrame(insertDict)
        df0 = dfRaw.append(insertDF, ignore_index=True, sort=False)
        df0copy = df0.copy() # a copy for gbif_fields_returned below

        ###########################################  RENAME & DELETE FIELDS ETC.
        ########################################################################
        df0.rename(mapper={"gbifID": "occ_id",
                           "decimalLatitude": "latitude",
                           "decimalLongitude": "longitude",
                           "eventDate": "occurrenceDate"}, inplace=True, axis='columns')
        df0.drop(["issue", "id"], inplace=True, axis=1)
        df0['coordinateUncertaintyInMeters'].replace(to_replace="UNKNOWN",
                                                     value=np.NaN, inplace=True)
        df0 = df0.astype({'coordinateUncertaintyInMeters': 'float',
                          'latitude': 'string', 'longitude': 'string'})
        df0['individualCount'].replace(to_replace="UNKNOWN", value=1,
                                       inplace=True)

        ############################  SUMMARY TABLE OF KEYS/FIELDS RETURNED (SMALL)
        ########################################################################
        # Count entries per attribute(column), reformat as new df with appropriate
        # columns.  Finally, insert into database.
        # NOTE: When pulling from df0copy, only a specified subset of keys are
        # assessed (keeper_keys).  For a more complete picture, alloccs must be
        # assessed.  That has historically been very slow.
        """ # Fastest, but least informative method for gbif_fields_returned
        newt = datetime.now()
        df0copy.where(df0copy != 'UNKNOWN', inplace=True)
        df_populated1 = pd.DataFrame(df0copy.count(axis=0).T.iloc[1:])
        #df_populated1['included(n)'] = df_populated1[0] # Can this be determined from alloccs?  Quickly?
        df_populated1['populated(n)'] = df_populated1[0]
        df_populated2 = df_populated1.filter(items=['included(n)', 'populated(n)'], axis='columns')
        df_populated2.index.name = 'attribute'
        df_populated2.to_sql(name='gbif_fields_returned', con=conn, if_exists='replace')
        print("Summarized fields returned: " + str(datetime.now() - newt))
        """
        # Slower, but more informative method for gbif_fields_returned
        '''
        The method below provides more information on values returned than the
        one above, but is slow.  Can it be improved to be faster?
        '''
        keys = [list(x.keys()) for x in alloccs]
        keys2 = set([])
        for x in keys:
            keys2 = keys2 | set(x)
        dfK = pd.DataFrame(index=keys2, columns=['included(n)', 'populated(n)'])
        dfK['included(n)'] = 0
        dfK['populated(n)'] = 0
        requestsummarytime1 = datetime.now()
        #####################################  START SLOW
        for t in alloccs:
            for y in t.keys():
                dfK.loc[y, 'included(n)'] += 1
                try:
                    int(t[y])
                    dfK.loc[y, 'populated(n)'] += 1
                except:
                    if t[y] == None:
                        pass
                    elif len(t[y]) > 0:
                        dfK.loc[y, 'populated(n)'] += 1
        print("Summarized fields returned: " + str(datetime.now() - requestsummarytime1))
        ######################################  END SLOW
        dfK.sort_index(inplace=True)
        dfK.index.name = 'attribute'
        dfK.to_sql(name='gbif_fields_returned', con=conn, if_exists='replace')

    ############################################################################
    #                         > 100,000 RECORDS (big)
    ############################################################################
    else:
        ########################################################## DOWNLOAD (big)
        ########################################################################
        # Make the data request using the download function.  Results are
        # emailed.
        # First, build a query list.  NoneType values cause problems, so only
        # add arguments if their value isn't NoneType.
        download_filters = ['taxonKey = {0}'.format(gbif_id)]
        if coordinate != None:
            download_filters.append('hasCoordinate = {0}'.format(coordinate))
        if country != None:
            download_filters.append('country = {0}'.format(country))
        if years != None:
            download_filters.append('year >= {0}'.format(years.split(",")[0]))
            download_filters.append('year <= {0}'.format(years.split(",")[1]))
        if months != None:
            download_filters.append('month >= {0}'.format(months.split(",")[0]))
            download_filters.append('month <= {0}'.format(months.split(",")[1]))
        if poly != None:
            download_filters.append('geometry within {0}'.format(poly))
        if geoIssue != None:
            download_filters.append('hasGeospatialIssue = {0}'.format(geoIssue))
        if latRange != None:
            download_filters.append('decimalLatitude >= {0}'.format(latRange.split(",")[0]))
            download_filters.append('decimalLatitude <= {0}'.format(latRange.split(",")[1]))
        if lonRange !=None:
            download_filters.append('decimalLongitude >= {0}'.format(lonRange.split(",")[0]))
            download_filters.append('decimalLongitude <= {0}'.format(lonRange.split(",")[1]))
        bigdown1 = datetime.now()
        d = occurrences.download(download_filters,
                                 pred_type='and',
                                 user = username,
                                 pwd = password,
                                 email = email)

        # Get the value of the download key
        dkey = d[0]

        # Now download the actual zip file containing the Darwin Core files
        # NOTE: The download can take a while to generate and is not immediately
        # available once the download_get command has been issued. Use a
        # while and try loop to make sure the download has succeeded.
        # The zipdownload variable will be a dictionary of the path,
        # the file size, and the download key unique code. It can be used
        # to change the file name, unzip the file, etc.
        print("Downloading Darwin Core Archive zip file for this species .....")
        gotit = None
        while gotit is None:
            try:
                zipdownload = occurrences.download_get(key=dkey, path=outDir)
                gotit = 1
                print("Download complete: " + str(datetime.now() - bigdown1))
            except:
                wait = datetime.now() - bigdown1
                if wait.seconds > 60*45:
                    gotit = 0
                    print("TIMED OUT -- attempting to proceed anyways")
                else:
                    gotit = None

        # Read the "occurrence.txt" file into a Pandas dataframe
        read1 = datetime.now()
        with DwCAReader(outDir + dkey + '.zip') as dwca:
            dfRaw = dwca.pd_read('occurrence.txt', low_memory=False)#, usecols=keeper_keys)

        df0 = dfRaw.filter(items=keeper_keys, axis=1)


        ###########################################  RENAME & DELETE FIELDS (big)
        ########################################################################
        df0.rename(mapper={"id": "occ_id",
                           "decimalLatitude": "latitude",
                           "decimalLongitude": "longitude",
                           "issue": "issues",
                           "eventDate": "occurrenceDate"}, inplace=True, axis='columns')
        df0['coordinateUncertaintyInMeters'].replace(to_replace="UNKNOWN",
                                                     value=np.NaN, inplace=True)
        df0['latitude'] = df0['latitude'].astype(str)
        df0['longitude'] = df0['longitude'].astype(str)
        df0['individualCount'].replace(to_replace="UNKNOWN", value=1,
                                       inplace=True)
        print("Downloaded and loaded records: " + str(datetime.now() - read1))

        ############################  SUMMARY TABLE OF KEYS/FIELDS RETURNED (big)
        ########################################################################
        # Count entries per atrribute(column), reformat as new df with appropriate
        # columns.  Finally, insert into database.
        feather = datetime.now()
        df_populated1 = pd.DataFrame(dfRaw.count(axis=0).T.iloc[1:])
        df_populated1['included(n)'] = len(dfRaw)
        df_populated1['populated(n)'] = df_populated1[0]
        df_populated2 = df_populated1.filter(items=['included(n)', 'populated(n)'], axis='columns')
        df_populated2.index.name = 'attribute'
        df_populated2.to_sql(name='gbif_fields_returned', con=conn, if_exists='replace')
        print("Summarized fields returned: " + str(datetime.now() - feather))

    ################################################# BRING IN NON-GBIF DATA
    ########################################################################
    """
    This functionality is forthcoming.
    """

    ############################################# SUMMARY OF VALUES RETURNED
    ########################################################################
    # Create a table for storing unique attribute values that came back.
    summarytime = datetime.now()
    summary = {'datums': ['WGS84'],
               'issues': set([]),
               'bases': [],
               'institutions': [],
               'collections': [],
               'datasets':[],
               'generalizations': set([]),
               'remarks': set([]),
               'establishment': set([]),
               'IDqualifier': set([]),
               'samplingProtocols': set([])}

    value_counts = {'bases': {},
                      'datums': {'WGS84': 0},
                      'issues': {},
                      'institutions': {},
                      'collections': {},
                      'datasets': {},
                      'samplingProtocols': {}}

    def get_vals(df, column_name):
        '''
        Return a set of unique values from a column
        '''
        stoat = df[column_name].unique()
        stoat = [str(x).split(";") for x in stoat]
        stoat1 = []
        for x in stoat:
            for y in x:
                if y == "" or y == None:
                    stoat1.append('UNKNOWN') # ? Keep?
                else:
                    stoat1.append(y)
        return set(stoat1)

    def set_value_counts(dataframe, groupby, key):
        '''
        Determine how many records there are with each value of an attribute.

        dataframe -- dataframe object to work on.
        groupby -- string column name to group by.
        key -- string key name in value_counts dict to populate a value for.
        '''
        group = dataframe['occ_id'].groupby(dataframe[groupby])
        skua = group.count()
        for x in skua.index:
            value_counts[key][x] = skua[x]

    # datums - ? - couldn't find this info in the table

    summary['issues'] = get_vals(df0, 'issues')
    set_value_counts(df0, 'issues', 'issues')

    summary['bases'] = get_vals(df0, 'basisOfRecord')
    set_value_counts(df0, 'basisOfRecord', 'bases')

    summary['institutions'] = get_vals(df0, 'institutionCode')
    set_value_counts(df0, 'institutionCode', 'institutions')

    summary['collections'] = get_vals(df0, 'collectionCode')
    set_value_counts(df0, 'collectionCode', 'collections')

    summary['datasets'] = get_vals(df0, 'datasetName')
    set_value_counts(df0, 'datasetName', 'datasets')

    try:
        summary['establishment'] = get_vals(df0, 'establishmentMeans')
    except:
        summary['establishment'] = ""

    summary['IDqualifier'] = get_vals(df0, 'identificationQualifier')

    summary['samplingProtocols'] = get_vals(df0, 'samplingProtocol')
    set_value_counts(df0, 'samplingProtocol', 'samplingProtocols')

    # Remove duplicates, make strings for entry into summary table of attributes
    cursor.executescript("""CREATE TABLE unique_values (step TEXT, field TEXT, vals TEXT);""")
    for x in summary.keys():
        vals = str(list(set(summary[x]))).replace('"', '')
        stmt = """INSERT INTO unique_values (step, field, vals)
                  VALUES ("request", "{0}", "{1}");""".format(x, vals)
        cursor.execute(stmt)

    # Store the value summary for the selected fields in a table.
    cursor.executescript("""CREATE TABLE pre_filter_value_counts
                            (attribute TEXT, value TEXT, count INTEGER);""")
    for x in value_counts.keys():
        attribute = value_counts[x]
        for y in value_counts[x].keys():
            z = value_counts[x][y]
            y = y.replace('"', "") # Double quotes occur and throw parsing error.
            frog = """INSERT INTO pre_filter_value_counts (attribute, value, count)
                      VALUES ("{0}", "{1}", "{2}")""".format(x,y,z)
            cursor.execute(frog)
    print("Created summary table of request results: " + str(datetime.now() - summarytime))


    ##########################################  SUMMARIZE SOURCES PRE FILTER
    ########################################################################
    moss = df0.groupby(['institutionCode', 'collectionCode', 'datasetName'])[['occ_id']].size()
    moss.to_sql(name='pre_filter_source_counts', con = conn, if_exists='replace')


    ###############################################  ADD SOME DEFAULT VALUES
    ########################################################################
    if default_coordUncertainty != False:
        df0.fillna(value={'coordinateUncertaintyInMeters': default_coordUncertainty},
                   inplace=True)
    df0.fillna(value={'individualCount': int(1)}, inplace=True)


    ################################################################  FILTER
    ########################################################################
    fiddlertime = datetime.now()
    if filt_coordUncertainty == 1:
        df1 = df0[pd.isnull(df0['coordinateUncertaintyInMeters']) == False]
    if filt_coordUncertainty == 0:
        df1 = df0
    df2 = df1[df1['coordinateUncertaintyInMeters'] <= filt_maxcoord]
    del df1
    df3 = df2[df2['collectionCode'].isin(filt_collection) == False]
    del df2
    df4 = df3[df3['institutionCode'].isin(filt_instit) == False]
    del df3
    df5 = df4[df4['basisOfRecord'].isin(filt_bases) == False]
    del df4
    df7 = df5[df5['samplingProtocol'].isin(filt_sampling) == False]
    del df5
    ''' ISSUES are more complex because multiple issues can be listed per record
    Method used is complex, but hopefully faster than simple iteration over all records
    '''
    df7.fillna(value={'issues': ""}, inplace=True)
    # Format of issues entries differ by method, change json format to email format
    if occ_count < 100000:
        df7['issues'] = [x.replace(', ', ';').replace('[', '').replace(']', '').replace("'", "")
                        for x in df7['issues']]
    unique_issue = list(df7['issues'].unique())
    violations = [x for x in unique_issue if len(set(str(x).split(";")) & set(filt_issues)) != 0] # entries that contain violations
    df8 = df7[df7['issues'].isin(violations) == False] # Records without entries that are violations.
    del df7
    print("Performed post-request filtering: " + str(datetime.now() - fiddlertime))
    newstime = datetime.now()

    # Create any new columns needed
    df8["remarks"] = df8['locality'] + ";" + df8['eventRemarks'] + ";" + df8['locationRemarks'] + ";" + df8['occurrenceRemarks']
    df8["taxon_id"] = taxon_id
    df8["request_id"] = gbif_req_id
    df8["filter_id"] = gbif_filter_id
    df8["retrievalDate"] = datetime.now()
    df8["detection_distance"] = det_dist
    df8["radius_meters"] = df8["detection_distance"] + df8["coordinateUncertaintyInMeters"]
    df8["source"] = "gbif"
    df8["doi_search"] = ""
    df8["weight"] = 10
    df8["weight_notes"] = ""
    df8.drop(labels=["scientificName", "eventRemarks", "locality",
                     "locationRemarks", "institutionID", "occurrenceRemarks"],
                     inplace=True, axis=1)
    print("Calculated new columns, deleted some too: " + str(datetime.now() - newstime))


    #########################################################  HANDLE DUPLICATES
    ############################################################################
    # Find out whether or not to drop duplicates.
    OKsql = """SELECT duplicates_OK FROM gbif_filters
                   WHERE filter_id = '{0}';""".format(gbif_filter_id)
    duplicates_OK = cursor2.execute(OKsql).fetchone()[0]
    conn2.commit()
    conn2.close()
    del cursor2

    if duplicates_OK == "False":
        df9 = drop_duplicates_latlongdate(df8)

    if duplicates_OK == "True":
        df9 = df8.copy()
        print("DUPLICATES ON LATITUDE, LONGITUDE, DATE-TIME INCLUDED")


    ###################################################  INSERT INTO DB (big)
    ########################################################################
    biggin = datetime.now()
    '''  # This is an alternate way to insert records
    sql1 = """INSERT INTO occurrences ('occ_id', 'taxon_id', 'source',
                                       'latitude', 'longitude',
                                       'coordinateUncertaintyInMeters',
                                       'occurrenceDate', 'request_id',
                                       'filter_id', 'generalizations',
                                       'remarks')
              VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);"""
    for x in df9.index:
        insert2 = [df9.loc[x,"id"], taxon_id, df9.loc[x,"source"],
                   df9.loc[x,"decimalLatitude"], df9.loc[x,"decimalLongitude"],
                   df9.loc[x,"coordinateUncertaintyInMeters"],
                   df9.loc[x,"eventDate"], request_id, filter_id,
                   df9.loc[x,"dataGeneralizations"], df9.loc[x,"remarks"]]
        cursor.execute(sql1, [(insert2)])
    conn.commit()
    '''
    df9.to_sql(name='occurrences', con = conn, if_exists='replace',
               chunksize=2000)
    sql_toad = '''SELECT AddGeometryColumn('occurrences', 'geom_xy4326', 4326,
                                           'POINT', 'XY');'''
    cursor.execute(sql_toad)
    print("Inserted records into table: " + str(datetime.now() - biggin))


    ################################## SUMMARY OF VALUES KEPT (FILTER; JSON)
    ########################################################################
    kepttime = datetime.now()
    summary = {'datums': ['WGS84'],
               'issues': set([]),
               'bases': [],
               'institutions': [],
               'collections': [],
               'generalizations': set([]),
               'remarks': set([]),
               'establishment': set([]),
               'IDqualifier': set([])}
    summary['issues'] = get_vals(df9, 'issues')
    summary['bases'] = get_vals(df9, 'basisOfRecord')
    summary['institutions'] = get_vals(df9, 'institutionCode')
    summary['collections'] = get_vals(df9, 'collectionCode')
    try:
        summary['establishment'] = get_vals(df9, 'establishmentMeans')
    except:
        summary['establishment'] = ""
    summary['IDqualifier'] = get_vals(df9, 'identificationQualifier')
    summary['samplingProtocols'] = get_vals(df9, 'samplingProtocol')

    # Remove duplicates, make strings for entry into summary table of attributes
    for x in summary.keys():
        vals = str(list(set(summary[x]))).replace('"', '')
        stmt = """INSERT INTO unique_values (step, field, vals)
                  VALUES ("filter", "{0}", "{1}");""".format(x, vals)
        cursor.execute(stmt)
    print("Summarized unique values retained: " + str(datetime.now() - kepttime))


    ################################################  MAKE POINT GEOMETRY COLUMN
    ############################################################################
    inserttime2 = datetime.now()
    try:
        sql2 = """UPDATE occurrences
                  SET geom_xy4326 = GeomFromText('POINT('||"longitude"||' '||"latitude"||')', 4326);"""
        cursor.execute(sql2)
    except Exception as e:
        print(e)

    ###### EVENTUALLY ADD CODE TO OVERIDE POLYGON GEOMETRY WITH FOOTPRINT
    ################ USE FOOTPRINTWKT HERE
    print("Updated occurrences table geometry column: " + str(datetime.now() - inserttime2))

    #############################################################  BUFFER POINTS
    ############################################################################
    # Buffer the xy points with the coordinate uncertainty
    # in order to create circles.  Create versions in albers and wgs84.  The
    # wgs84 version will be used in plotting with Basemap.  Buffer radius is
    # the sum of detectiondistance from requests.taxa_concepts and
    # coordinate uncertainty in meters here.
    buffertime1 = datetime.now()
    sql_det = """
            UPDATE occurrences
            SET detection_distance = {0};

            UPDATE occurrences
            SET radius_meters = detection_distance + coordinateUncertaintyInMeters;
            """.format(det_dist)#(requestsDB, det_dist)
    cursor.executescript(sql_det)

    sql_buf = """
            /* Transform to albers (5070) and apply buffer */
            ALTER TABLE occurrences ADD COLUMN polygon_5070 BLOB;

            UPDATE occurrences SET polygon_5070 = Buffer(Transform(geom_xy4326,
                                                                    5070),
                                                          radius_meters);

            SELECT RecoverGeometryColumn('occurrences', 'polygon_5070', 5070,
                                         'POLYGON', 'XY');

            /* Transform back to WGS84 so it can be displayed in iPython */
            ALTER TABLE occurrences ADD COLUMN polygon_4326 BLOB;

            UPDATE occurrences SET polygon_4326 = Transform(polygon_5070, 4326);

            SELECT RecoverGeometryColumn('occurrences', 'polygon_4326', 4326,
                                         'POLYGON', 'XY');
    """
    cursor.executescript(sql_buf)
    print("Buffered points: " + str(datetime.now() - buffertime1))

    ###############################################################  EXPORT MAPS
    ############################################################################
    # Export occurrence circles as a shapefile (all seasons)
    try:
        exportSHP(database=spdb, table='occurrences', column='polygon_4326',
                  outFile = outDir + summary_name + '_polygons')
    except Exception as e:
        print('\n Failed to create a point-polygon shapefile -- \n' + str(e))
    conn.commit()
    conn.close()
    print("\nRecords saved in {0}".format(spdb))
