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
##################    is the above in the right place?  should it be after code below? ???????????????????????????????




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
