configDir = "T:/Data/"  # Path to folder where you saved your wildlifeconfig file.
working_directory = "T:/temp/"
query_name = 'withRTest'
taxon_info = {'EBIRD_ID': 'Yellow-billed Cuckoo',
    'GBIF_ID': 2496287,
    'ID': 'TestCuckoo',
    'TAXON_EOO': 'POLYGON ((-84.09680233298448 36.69265225442667, '
          '-84.07962135716329 34.5561660300382, -84.07962135716329 '
          '34.5561660300382, -80.25685423694925 34.65515526072436, '
          '-81.15026497965096 36.71331438415306, -84.09680233298448 '
          '36.69265225442667))'}
filter_set = {'bases_omit': '',
 'collection_codes_omit': '',
 'country': 'US',
 'datasets_omit': '',
 'default_coordUncertainty': 1000,
 'duplicates_OK': False,
 'geoissue': False,
 'get_dwca': False,
 'has_coordinate_uncertainty': True,
 'institutions_omit': '',
 'issues': '',
 'lat_range': '27,41',
 'lon_range': '-89,-75',
 'max_coordinate_uncertainty': 10000,
 'months_range': '1,12',
 'name': 'test_filters_1',
 'query_polygon': 'POLYGON ((-82.74809573102132 36.96082629937069, '
                  '-85.0932989306133 35.63154639485496, -81.0987220521874 '
                  '33.56697226279766, -79.4235769096217 36.34054727735634, '
                  '-79.4235769096217 36.34054727735634, -82.74809573102132 '
                  '36.96082629937069))',
 'sampling_protocols_omit': '',
 'use_taxon_geometry': True,
 'years_range': '2000,2020'}

##########################################           def get_GBIF_records(taxon_info, filter_set, working_directory, username, password, email, dwca_download):
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
use_taxon_geometry -- True or False to use geometry saved with taxon concept when
    filtering records.  Defaults to 'True'.
dwca_download -- True or False.  False uses the API, which only works when there are
    fewer than a few 100,000 records.  True uses the download method involving
    your GBIF account and email.  Default is True.  Note: False does not
    provide a download DOI.
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
timestamp = datetime.now()

# TAXON INFO >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
gbif_id = taxon_info["GBIF_ID"]
#det_dist = concept[3]
taxon_polygon =taxon_info["TAXON_EOO"]


#  PREP FILTERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
years = filter_set["years_range"]
print(years)
months = filter_set["months_range"]
print(months)
latRange = filter_set["lat_range"]
print(decimalLatitude)
lonRange = filter_set["lon_range"]
print(decimalLongitude)
coordinate = filter_set["has_coordinate_uncertainty"]
print(coordinate)
geoIssue = filter_set["geoissue"]
print(geoIssue)
country = filter_set["country"]
print(country)
query_polygon = filter_set["query_polygon"]
print(query_polygon)
use_taxon_geometry = filter_set["use_taxon_geometry"]
print(use_taxon_geometry)
dwca_download = filter_set["get_dwca"]
print(dwca_download)

# Sort out geometries
# A geometry could also be stated for the species, assess what to do
# It could also be that user opted not to use species geometry.
if use_taxon_geometry == False:
    taxon_polygon = None
if query_polygon == None and taxon_polygon == None:
    poly = None
elif query_polygon != None and taxon_polygon == None:
    poly = query_polygon
elif query_polygon == None and taxon_polygon != None:
    poly = taxon_polygon
elif query_polygon != None and taxon_polygon != None:
    # Get/use the intersection of the two polygons
    filter_polygon = shapely.wkt.loads(query_polygon)
    sp_polygon = shapely.wkt.loads(taxon_polygon)
    poly_intersection = filter_polygon.intersection(sp_polygon)
    poly = shapely.wkt.dumps(poly_intersection)

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

print("Got request params and sorted out geometry constraints: " + str(datetime.now() - timestamp))
timestamp = datetime.now()


# GET RECORD COUNT >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
occ_count=occ_search["count"]
print(str(occ_count) + " records available")


# API QUERY >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if dwca_download == False:
    # Get records in batches, saving into master list.
    all_jsons = []
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
        all_jsons = all_jsons + occs

    # Load records into a data frame
    dfRaw = pd.DataFrame(columns=keeper_keys)
    insertDict = {}
    for x in keeper_keys:
        insertDict[x] = []
    for x in all_jsons:
        present_keys = list(set(x.keys()) & set(keeper_keys))
        for y in present_keys:
            insertDict[y] = insertDict[y] + [str(x[y])]
        missing_keys = list(set(keeper_keys) - set(x.keys()))
        for z in missing_keys:
            insertDict[z] = insertDict[z] + ["UNKNOWN"]
    insertDF = pd.DataFrame(insertDict)
    df0 = dfRaw.append(insertDF, ignore_index=True, sort=False)
    df0copy = df0.copy() # a copy for gbif_fields_returned below

    # Manage the columns in the data frame
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

    print("Downloaded records: " + str(datetime.now() - timestamp))
    timestamp = datetime.now()

#### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< RESUME HERE !!!!!

    ############################  SUMMARY TABLE OF KEYS/FIELDS RETURNED (SMALL)
    ########################################################################
    # Count entries per attribute(column), reformat as new df with appropriate
    # columns.  Finally, insert into database.
    # NOTE: When pulling from df0copy, only a specified subset of keys are
    # assessed (keeper_keys).  For a more complete picture, all_jsons must be
    # assessed.  That has historically been very slow.
    """ # Fastest, but least informative method for gbif_fields_returned
    newt = datetime.now()
    df0copy.where(df0copy != 'UNKNOWN', inplace=True)
    df_populated1 = pd.DataFrame(df0copy.count(axis=0).T.iloc[1:])
    #df_populated1['included(n)'] = df_populated1[0] # Can this be determined from all_jsons?  Quickly?
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
    keys = [list(x.keys()) for x in all_jsons]
    keys2 = set([])
    for x in keys:
        keys2 = keys2 | set(x)
    dfK = pd.DataFrame(index=keys2, columns=['included(n)', 'populated(n)'])
    dfK['included(n)'] = 0
    dfK['populated(n)'] = 0
    requestsummarytime1 = datetime.now()
    #####################################  START SLOW
    for t in all_jsons:
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
#                   DARWIN CORE ARCHIVE METHOD (email of .zip)
############################################################################
if dwca_download == True:
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

    # Read the relevant files from within the darwin core archive
    read1 = datetime.now()
    with DwCAReader(outDir + dkey + '.zip') as dwca:
        dfRaw = dwca.pd_read('occurrence.txt', low_memory=False)
        citations = dwca.open_included_file('citations.txt').read()
        rights = dwca.open_included_file('rights.txt').read()
        doi = dwca.metadata.attrib["packageId"]

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
if dwca_download == True:
    df8["GBIF_download_doi"] = doi
else:
    df8["GBIF_download_doi"] = "bypassed"
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

############################################################## DWCA METADATA
############################################################################
if dwca_download == True:
    # Store the value summary for the selected fields in a table.
    lizardtime = datetime.now()
    cursor.executescript("""CREATE TABLE GBIF_download_info
                            (download_key TEXT, doi TEXT, citations TEXT,
                             rights TEXT);""")
    cursor.execute('''INSERT INTO GBIF_download_info (doi, citations,
                                                      rights, download_key)
                      VALUES ("{0}", "{1}", "{2}", "{3}")'''.format(doi,
                                                              citations,
                                                              rights,
                                                              dkey))
    print("Stored GBIF Download DOI etc.: " + str(datetime.now() - lizardtime))
