import pandas as pd
import sqlite3
from datetime import datetime
import wrangler_functions as functions

ebird_data = pd.read_csv("T:/Temp/ebird.csv", index_col=0)
gbif_data = pd.read_csv("T:/Temp/gbif.csv", index_col=0)
working_directory = "T:/Occurrence_Records/"
query_name = 'withRTest1'
filter_set = {'bases_omit': None,
 'collection_codes_omit': None,
 'country': 'US',
 'datasets_omit': None,
 'default_coordUncertainty': 1000,
 'duplicates_OK': False,
 'geoissue': False,
 'get_dwca': True,
 'has_coordinate_uncertainty': True,
 'institutions_omit': None,
 'issues_omit': None,
 'lat_range': '27,41',
 'lon_range': '-89,-75',
 'max_coordinate_uncertainty': 10000,
 'months_range': '9,11',
 'name': 'test_filters_1',
 'query_polygon': 'POLYGON ((-82.74809573102132 36.96082629937069, '
                  '-85.0932989306133 35.63154639485496, -81.0987220521874 '
                  '33.56697226279766, -79.4235769096217 36.34054727735634, '
                  '-79.4235769096217 36.34054727735634, -82.74809573102132 '
                  '36.96082629937069))',
 'sampling_protocols_omit': None,
 'use_taxon_geometry': True,
 'years_range': '2015,2020'}

taxon_info = {'EBIRD_ID': 'Yellow-billed Cuckoo', 'GBIF_ID': 2496287, 'ID': 'TestCuckoo',
'TAXON_EOO': 'POLYGON ((-84.09680233298448 36.69265225442667, '
                          '-84.07962135716329 34.5561660300382, -84.07962135716329 '
                              '34.5561660300382, -80.25685423694925 34.65515526072436, '
                              '-81.15026497965096 36.71331438415306, -84.09680233298448 '
                              '36.69265225442667))'}

def apply_filters(ebird_data, gbif_data, filter_set, taxon_info, working_directory, query_name):
    '''
    Summarizes the values in the data frames, apply filters, summarize what
        values persisted after filtering.  Insert results into the output db.

    PARAMETERS
    ebird_data : a data frame of records from eBird
    gbif_data : a data frame of records from GBIF
    output_database : path to the output database
    filter_set : the filter set dictionary
    taxon_info : the taxon information dictionary

    RETURNS
    filtered_records : a data frame of filtered records.
    '''

    # Create or connect to the database
    output_database = working_directory + query_name + ".sqlite"
    conn = sqlite3.connect(output_database, isolation_level='DEFERRED')
    cursor = conn.cursor()


    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  COMBINE DATA FRAMES
    # Concatenate the gbif and ebird tables
    records3 = pd.concat([ebird_data, gbif_data])


    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SUMMARIZE VALUES
    timestamp = datetime.now()
    # Make a list of columns to summarize values from
    do_not_summarize = ['decimalLatitude', 'decimalLongitude', 'GBIF_download_doi',
                        'coordinateUncertaintyInMeters', 'detection_distance_m',
                        'eventDate', 'eventRemarks', 'filter_set_name',
                        'footprintSRS', 'footprintWKT', 'gbif_id', 'ebird_id',
                        'general_remarks', 'georeferencedBy', 'habitat',
                        'georeferenceRemarks', 'identificationQualifier',
                        'identifiedby', 'identifiedRemarks', 'individualCount',
                        'informationWitheld', 'locality', 'locationAccordingTo',
                        'locationRemarks', 'occurrenceRemarks', 'radius_m',
                        'record_id', 'recordedBy', 'retrieval_date',
                        'taxonConceptID', 'verbatimLocality', 'weight', 'weight_notes']

    # Make a function to do the summarizing
    def summarize_values(dataframe, step):
        """
        Loops through columns and gets a count of unique values.  Packages in a df.
        """
        attributes = []
        summarize = [x for x in dataframe.columns if x not in do_not_summarize]
        for column in summarize:
            value_count = dataframe['record_id'].groupby(dataframe[column]).count()
            value_df = (pd.DataFrame(value_count)
                        .reset_index()
                        .rename({'record_id': step, column: 'value'}, axis=1))
            value_df['attribute'] = column
            value_df = value_df[["attribute", "value", step]]
            if value_df.empty == False:
                attributes.append(value_df)
        result = pd.concat(attributes)
        return result

    # Store summary in a dataframe
    acquired = summarize_values(dataframe=records3, step='acquired')
    print("Summarized values acquired: " + str(datetime.now() - timestamp))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  FILTER
    timestamp = datetime.now()

    # Some filters to be prepped for use
    for x in ['bases_omit', 'collection_codes_omit', 'datasets_omit',
              'institutions_omit', 'issues_omit', 'sampling_protocols_omit']:
        if filter_set[x] == None:
            filter_set[x] = []

    if filter_set['has_coordinate_uncertainty'] == 1:
        records4 = records3[pd.isnull(records3['coordinateUncertaintyInMeters']) == False]
    if filter_set['has_coordinate_uncertainty'] == 0:
        records4 = records3

    records5 = (records4[records4['coordinateUncertaintyInMeters'] <= filter_set['max_coordinate_uncertainty']]
                [lambda x: x['collectionCode'].isin(filter_set['collection_codes_omit']) == False]
                [lambda x: x['institutionID'].isin(filter_set['institutions_omit']) == False]
                [lambda x: x['basisOfRecord'].isin(filter_set['bases_omit']) == False]
                [lambda x: x['samplingProtocol'].isin(filter_set['sampling_protocols_omit']) == False]
                [lambda x: x['datasetName'].isin(filter_set['datasets_omit']) == False]
                )

    ''' ISSUES are more complex because multiple issues can be listed per record
    Method used is complex, but hopefully faster than simple iteration over all records
    '''
    records5.fillna(value={'issues': ""}, inplace=True)
    # Format of issues entries differ by method, change json format to email format
    if filter_set['get_dwca'] == True:
        records5['issues'] = [x.replace(', ', ';').replace('[', '').replace(']', '').replace("'", "")
                              for x in records5['issues']]
    unique_issue = list(records5['issues'].unique())
    violations = [x for x in unique_issue if len(set(str(x).split(";")) & set(filter_set['issues_omit'])) != 0] # entries that contain violations
    records6 = records5[records5['issues'].isin(violations) == False] # Records without entries that are violations.
    print("Performed filtering: " + str(datetime.now() - timestamp))


    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  POPULATE SOME COLUMNS
    if filter_set["default_coordUncertainty"] != None:
        records6.fillna(value={'coordinateUncertaintyInMeters': filter_set["default_coordUncertainty"]},
                        inplace=True)
    records6.fillna(value={'individualCount': int(1)}, inplace=True)
    records6["radius_meters"] = records6["detection_distance_m"] + records6["coordinateUncertaintyInMeters"]
    records6["weight"] = 10
    records6["weight_notes"] = ""


    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  REMOVE SPACE-TIME DUPLICATES
    # Prep some columns by changing data type
    records6 = (records6.astype({'decimalLatitude': 'str', 'decimalLongitude': 'str'})
                .reset_index(drop=True))

    if filter_set["duplicates_OK"] == False:
        records7 = functions.drop_duplicates_latlongdate(records6)

    if filter_set["duplicates_OK"] == True:
        records7 = records6.copy()
        print("DUPLICATES ON LATITUDE, LONGITUDE, DATE-TIME INCLUDED")


    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SPATIAL FILTERING
    # Spatial filtering happens in the get functions (ebird and gbif), not here.


    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SUMMARIZE VALUES AGAIN
    # Store summary in a dataframe
    retained = summarize_values(dataframe=records7, step='retained')

    # Concat acquired and retained data frames
    summary_df = pd.merge(retained, acquired, on=['attribute', 'value'], how='inner')

    # Calculate a difference column
    summary_df['removed'] = summary_df['acquired'] - summary_df['retained']
    summary_df = summary_df[['attribute', 'value', 'acquired', 'removed',
                             'retained']]

    # Save the summary in the output database
    summary_df.to_sql(name='attribute_value_counts', con = conn, if_exists='replace')


    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SAVE
    records6.to_sql(name='occurrence_records', con = conn, if_exists='replace')
    return None

apply_filters(ebird_data, gbif_data, filter_set, taxon_info, working_directory, query_name)
