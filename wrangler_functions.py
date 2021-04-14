# occurrece records table datatypes
attribute_data_types = {'GBIF_download_doi': 'str', 'accessRights': 'str',
             'basisOfRecord': 'str', 'bibliographicCitation': 'str',
             'collectionCode': 'str', 'coordinateUncertaintyInMeters': 'float',
             'dataGeneralizations': 'str', 'datasetName': 'str',
             'decimalLatitude': 'str', 'decimalLongitude': 'str',
             'detection_distance_m': 'int', 'ebird_id': 'str',
             'establishmentMeans': 'str', 'eventDate': 'str', 'eventRemarks': 'str',
             'filter_set_name': 'str', 'footprintSRS': 'str', 'footprintWKT': 'str',
             'gbif_id': 'str', 'general_remarks': 'str', 'georeferenceProtocol': 'str',
             'georeferenceRemarks': 'str', 'georeferenceVerificationStatus': 'str',
             'georeferencedBy': 'str', 'habitat': 'str',
             'identificationQualifier': 'str', 'identifiedBy': 'str',
             'identifiedRemarks': 'str', 'individualCount': 'int',
             'informationWitheld': 'str', 'institutionID': 'str', 'issues': 'str',
             'license': 'str', 'locality': 'str', 'locationAccordingTo': 'str',
             'locationRemarks': 'str', 'modified': 'str', 'occurrenceRemarks': 'str',
             'occurrenceStatus': 'str', 'radius_meters': 'float', 'record_id': 'str',
             'recordedBy': 'str', 'retrieval_date': 'str', 'samplingProtocol': 'str',
             'scientificName': 'str', 'source': 'str', 'taxonConceptID': 'str',
             'taxon_info_name': 'str', 'verbatimLocality': 'str', 'weight': 'int',
             'weight_notes': 'str'}

# Core functions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def build_output_database(output_database):
    """
    Create a database for storing occurrence and taxon concept data.
    The column names that are "camel case" are Darwin Core attributes, whereas
    lower case names containing "_" between words are not.  Not all Darwin
    Core attributes are included.  Only those that could be useful for filtering
    and assigning weights are included.

    Parameters
    ----------
    output_database : Path for sqlite database to create; string.

    Returns
    -------
    Nothing.
    """
    import os
    import sqlite3

    # Delete the database if it already exists
    if os.path.exists(output_database):
        os.remove(output_database)

    # Create or connect to the database
    conn = sqlite3.connect(output_database, isolation_level='DEFERRED')
    cursor = conn.cursor()

    # Create a table for occurrence records.
    sql_cdb = """
            CREATE TABLE IF NOT EXISTS occurrence_records (
                    record_id INTEGER NOT NULL PRIMARY KEY,
                    filter_set_name TEXT NOT NULL,
                    taxon_info_name TEXT NOT NULL,
                    gbif_id TEXT,
                    ebird_id TEXT,
                    source TEXT NOT NULL,
                    retrieval_date TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
                    detection_distance_m INTEGER,
                    radius_meters INTEGER,
                    GBIF_download_doi TEXT,
                    general_remarks TEXT,
                    weight INTEGER DEFAULT 10,
                    weight_notes TEXT,

                    accessRights TEXT,
                    bibliographicCitation TEXT,
                    basisOfRecord TEXT,
                    collectionCode TEXT,
                    coordinateUncertaintyInMeters INTEGER,
                    dataGeneralizations TEXT,
                    datasetName TEXT,
                    decimalLatitude TEXT,
                    decimalLongitude TEXT,
                    establishmentMeans TEXT,
                    eventDate TEXT,
                    eventRemarks TEXT,
                    footprintWKT TEXT,
                    footprintSRS TEXT,
                    georeferencedBy TEXT,
                    georeferenceProtocol TEXT,
                    georeferenceVerificationStatus TEXT,
                    georeferenceRemarks TEXT,
                    habitat TEXT,
                    identifiedBy TEXT,
                    identifiedRemarks TEXT,
                    identificationQualifier TEXT,
                    individualCount INTEGER DEFAULT 1,
                    informationWitheld TEXT,
                    institutionID TEXT,
                    issues TEXT,
                    license TEXT,
                    locality TEXT,
                    locationAccordingTo TEXT,
                    locationRemarks TEXT,
                    modified TEXT,
                    occurrenceStatus TEXT,
                    occurrenceRemarks TEXT,
                    recordedBy TEXT,
                    samplingProtocol TEXT,
                    scientificName TEXT,
                    taxonConceptID INTEGER NOT NULL,
                    verbatimLocality TEXT,
                        FOREIGN KEY (taxonConceptID) REFERENCES taxa(taxonConceptID)
                        ON UPDATE RESTRICT
                        ON DELETE NO ACTION);
    """
    cursor.executescript(sql_cdb)
    return

def get_EBD_records(taxon_info, filter_set, working_directory, EBD_file, query_name):
    '''
    Gets eBird records from a copy of the Ebird Basic Dataset that you acquired.
    Primarily runs R code that uses the Auk package to query the data set in
    an efficient manner.  Some filters can be applied during the query, but
    others have to be applied to the query results.  Date and bounding box
    filters require quite a bit of preparation and conditions.

    Parameters
    ----------
    taxon_info : your taxon concept; dictionary
    filter_set : name of the filter set to apply; dictionary
    working_directory : path to use for table of filtered query results; string
    EBD_file : path to your downloaded copy of the Ebird Basic Dataset; string
    query_name : the name you chose for your query; string

    Returns
    -------
    Data frame of eBird records
    '''
    import rpy2.robjects as robjects
    import rpy2.robjects.packages as rpackages
    from rpy2.robjects.vectors import StrVector
    import pandas as pd
    from datetime import datetime
    import sqlite3
    import geopandas as gpd
    import shapely
    import numpy as np

    # import R's utility package, select a mirror for R packages
    utils = rpackages.importr('utils')
    utils.chooseCRANmirror(ind=1) # select the first mirror in the list

    # R packages to load
    packnames = ('sf', 'auk', 'lubridate', 'tidyverse')
    names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
    if len(names_to_install) > 0:
        utils.install_packages(StrVector(names_to_install))

    # Some file names
    queried_ebd = working_directory + "tmp_ebd.txt"
    processed_ebd = working_directory + query_name + ".csv"
    output_database = working_directory + query_name + '.sqlite'

    # Replace None values in fitler_set with "" to fit R code.
    for x in filter_set.keys():
        if filter_set[x] == None:
            filter_set[x] = ""

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< R CODE
    code = '''
    EBD_file <- "{0}"
    queried_ebd <- "{1}"
    processed_ebd <- "{2}"
    species <- "{3}"
    country <- "{4}"
    months_range <- "{5}"
    years_range <- "{6}"
    lon_range <- "{7}"
    lat_range <- "{8}"
    max_coordinate_uncertainty <- {9}
    taxon_polygon <- "{10}"
    query_polygon <- "{11}"

    # Auk uses filters that are compiled and incorporated into a query. This poses
    # a challenge for dynamic filtering where filters may or may not be used.  We
    # have to set defaults.
    library(auk)
    library(tidyverse)
    library(sf)
    library(lubridate)
    starttime = Sys.time() # Runtime has been 30 min

    # prep dates ---------------------------------------------------------------
    # auk doesn't allow filtering on months AND year with read_ebd; they have to be
    #   done separately, one with auk filters and the other after with dplyr.  I
    #   chose to do the year filtering with auk to minimize size of returned tibble.
    #   This all requires formatting dates as text correctly
    # format start month
    if (months_range != "") {{
      if (as.numeric(strsplit(months_range, ",")[[1]][1]) < 10) {{
        start_month <- paste(c("0", strsplit(months_range, ",")[[1]][1]), collapse="")
        }} else {{
        start_month <- strsplit(months_range, ",")[[1]][1]}}
      start_month <- str_trim(start_month)

      # format end month
      if (as.numeric(strsplit(months_range, ",")[[1]][2]) < 10) {{
        end_month <- paste(c("0", strsplit(months_range, ",")[[1]][2]), collapse="")
        }} else {{
        end_month <- strsplit(months_range, ",")[[1]][2]}}
      end_month <- str_trim(end_month)

      # create vector of ok months for filtering with dplyr
      ok_months <- seq(as.numeric(start_month), as.numeric(end_month))
    }}

    # pull out start and end years
    if (years_range != "") {{
    start_yr <- str_trim(strsplit(years_range, ",")[[1]][1])
    end_yr <- str_trim(strsplit(years_range, ",")[[1]][2])
    }}

    # define data filter according to different scenarios
    if (months_range == "" && years_range == "") {{
      # get all dates
      date_filter <- c("*-01-31", "*-12-31")
    }} else if (months_range != "" && years_range != "") {{
      # get certain months and years.  we have to find the last possible day.
      end_day <- lubridate::days_in_month(paste(c(end_yr, "-", end_month, "-01"),
                                                collapse=""))
      date_filter <- c(paste(c(start_yr, "-", start_month, "-01"), collapse=""),
                       paste(c(end_yr, "-", end_month, "-", end_day), collapse=""))
    }} else if (months_range == "" && years_range != "") {{
      # get all months from certain years
      date_filter <- c(paste(c(start_yr, "-01-01"), collapse=""),
                       paste(c(end_yr, "-12-31"), collapse=""))
    }} else if (months_range != "" && years_range == "") {{
      # get certain months from all years.  we have to find the last possible day.
      yr <- year(today())
      end_day <- lubridate::days_in_month(paste(c(yr, "-", end_month, "-01"),
                                                collapse=""))
      date_filter <- c(paste(c("*-", start_month, "-01"), collapse=""),
                       paste(c("*-", end_month, "-", end_day), collapse=""))
    }}

    # prep bounding box --------------------------------------------------------
    bbox <- NULL
    if (query_polygon == "" && taxon_polygon == "") {{
        bbox <- NULL
    }} else if (query_polygon != "" && taxon_polygon == "") {{
        bbox <- st_bbox(st_as_sfc(query_polygon))
    }} else if (query_polygon == "" && taxon_polygon != "") {{
        bbox <- st_bbox(st_as_sfc(taxon_polygon))
    }} else if (query_polygon != "" && taxon_polygon != "") {{
        # Get/use the intersection of the two polygons
        filter_polygon <- st_as_sfc(query_polygon)
        sp_polygon <- st_as_sfc(taxon_polygon)
        bbox <- st_bbox(st_intersection(filter_polygon, sp_polygon))
    }}

    # prep bounding box vector for filter if lat and lon ranges were provided,
    #   and if other polygons were not
    if (lat_range == "" || lon_range == "") {{
      null_box <- TRUE
    }} else {{
      null_box <- FALSE
    }}

    if (is.null(bbox) == TRUE && null_box == FALSE) {{
      lat_min <- as.numeric(strsplit(lat_range, ",")[[1]][1])
      lat_max <- as.numeric(strsplit(lat_range, ",")[[1]][2])
      lng_min <- as.numeric(strsplit(lon_range, ",")[[1]][1])
      lng_max <- as.numeric(strsplit(lon_range, ",")[[1]][2])
      bbox <- c(lng_min, lat_min, lng_max, lat_max)
      names(bbox) <- c("xmin", "ymin", "xmax", "ymax")
      attr(bbox, "class") = "bbox"
    }}

    # prep country -------------------------------------------------------------
    if (country == "") {{country <- "US"}}

    # prep distance ------------------------------------------------------------
    #   a gps precision for eBird checklists must be assumed, since not given, for
    #   estimation of coordinateUncertaintyInMeters
    EBD_gps_precision <- 10

    # account for gps precision in distance filter.  error could exist on either
    #   end of a straight line path, so double precision when subtracting.
    max_distance <- as.integer(ceiling((max_coordinate_uncertainty-2*EBD_gps_precision)/1000))  ### NOT WORKING

    # query --------------------------------------------------------------------
    records0 <- EBD_file %>%
      # 1. reference file
      auk_ebd() %>%

      # 2. define filters
      auk_species(species=c(species)) %>%
      auk_date(date=date_filter) %>%
      auk_country(country=country) %>%
      auk_bbox(bbox=bbox) %>%
      auk_distance(distance=c(0, max_distance)) %>%

      # 3. run filtering
      auk_filter(file = queried_ebd, overwrite = TRUE) %>%

      # 4. read text file into r data frame
      read_ebd()

    # prep data frame for python -----------------------------------------------
    # add column for eBird species code
    ebird_code <- select(filter(ebird_taxonomy, common_name==species),
                         species_code)[[1]]
    ebd_data <- records0 %>%
                 mutate(eBird_sp_code = ebird_code,
                        retrieval_date = auk_ebd_version(EBD_file)[1][1]) %>%
                 select(eBird_sp_code, global_unique_identifier, checklist_id,
                        project_code, last_edited_date, common_name,
                        observation_count, locality, latitude, longitude,
                        observation_date, observer_id, effort_distance_km,
                        protocol_type, effort_area_ha, trip_comments,
                        species_comments) %>%
                 mutate(effort_distance_m = as.numeric(effort_distance_km)*1000) %>%
                 mutate(coordinateUncertaintyInMeters = effort_distance_m + (2*EBD_gps_precision)) %>%
                 filter(month(observation_date) %in% ok_months) %>%
                 write_csv(processed_ebd)

    endtime = Sys.time()
    print(endtime - starttime)
    '''.format(EBD_file, queried_ebd, processed_ebd, taxon_info['EBIRD_ID'],
    filter_set["country"], filter_set["months_range"], filter_set["years_range"],
    filter_set["lon_range"], filter_set["lat_range"],
    filter_set["max_coordinate_uncertainty"], taxon_info["TAXON_EOO"],
    filter_set["query_polygon"])

    # Run code
    timestamp = datetime.now()
    robjects.r(code)
    print("Ran EBD query with Auk: " + str(datetime.now() - timestamp))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  READ OUTPUT
    records0 = pd.read_csv(processed_ebd)
    '''
    This should eventually be usable (or something similar) to avoid having to write
    the R data frame and then read it in with pandas.  It is supposed to be
    possible to convert robject data frames to pandas dataframes but the rpy2
    available from conda 2.x doesn't actually work.
    # *****************************************************************************
    # Retrieve the filtered ebird data frame  ---- this should work, but doesn't
    rdf = robjects.globalenv['ebd_data']
    # Using a conversion context in which the pandas conversion is
    # added to the default conversion rules, the rpy2 object
    # (an R data frame) is converted to a pandas data frame.
    from rpy2.robjects import pandas2ri
    robjects.pandas2ri.activate() # should automatically convert r dataframe to pandas
    from rpy2.robjects import default_converter
    from rpy2.robjects.conversion import localconverter
    with localconverter(robjects.default_converter + pandas2ri.converter):
        records0 = robjects.conversion.ri2py(rdf)
    '''

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  APPLY SPATIAL FILTER
    timestamp = datetime.now()
    # Make data frame spatial
    gdf = gpd.GeoDataFrame(records0,
                           geometry=gpd.points_from_xy(records0['longitude'],
                                                       records0['latitude']))

    #   It could be that user opted not to use species geometry.
    if filter_set['use_taxon_geometry'] == False:
        EOO = None

    # A geometry could be stated for the species, assess what to do
    EOO = taxon_info["TAXON_EOO"]
    AOI = filter_set["query_polygon"]
    if AOI is None and EOO is None:
        filter_polygon = None
    elif AOI is not None and EOO is None:
        filter_polygon = shapely.wkt.loads(AOI)
    elif AOI is None and EOO is not None:
        filter_polygon = shapely.wkt.loads(EOO)
    elif AOI is not None and EOO is not None:
        # Get/use the intersection of the two polygons in this case
        AOI_polygon = shapely.wkt.loads(AOI)
        EOO_polygon = shapely.wkt.loads(EOO)
        filter_polygon = AOI_polygon.intersection(EOO_polygon)
    print("Calculated the spatial filter polygon: " + str(datetime.now() - timestamp))

    # Find which records have coordinates that fall within the polygon
    timestamp = datetime.now()
    if filter_polygon is not None:
        gdf = gdf[gdf["geometry"].within(filter_polygon)]
    print("Applied spatial filter: " + str(datetime.now() - timestamp))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SUMMARIZE
    timestamp = datetime.now()
    df_populated1 = pd.DataFrame(records0.count(axis=0).T.iloc[1:])
    df_populated1['included(n)'] = len(records0)
    df_populated1['populated(n)'] = df_populated1[0]
    df_populated2 = df_populated1.filter(items=['included(n)', 'populated(n)'], axis='columns')
    df_populated2.index.name = 'attribute'
    conn = sqlite3.connect(output_database, isolation_level='DEFERRED')
    df_populated2.to_sql(name='eBird_fields_returned', con=conn, if_exists='replace')
    print("Summarized fields returned: " + str(datetime.now() - timestamp))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  PREP FOR OUTPUT
    timestamp = datetime.now()
    ''' DELETE THIS  ??????????????????????????????????????????????????????????
    conn = sqlite3.connect(output_database, isolation_level='DEFERRED')
    cursor = conn.cursor()
    schema = conn.execute("PRAGMA table_info (occurrence_records);").fetchall()
    column_names = [x[1] for x in schema]
    '''
    # Rename columns
    gdf = gdf.rename({'eBird_sp_code': 'ebird_id',
             'global_unique_identifier': 'record_id',
             'latitude': 'decimalLatitude',
             'longitude': 'decimalLongitude',
             'observation_count': 'individualCount',
             'observation_date': 'eventDate',
             'project_code': 'collectionCode',
             'protocol_type': 'samplingProtocol',
             'species_comments': 'identifiedRemarks',
             'trip_comments': 'eventRemarks'}, axis=1)
    # Drop columns
    records1 = gdf.filter(list(attribute_data_types.keys()), axis=1)
    # Populate columns
    records1["institutionID"] = "clo"
    records1["collectionCode"] = "EBIRD"
    records1["datasetName"] = "EBD"
    records1["source"] = "eBird"
    records1["basisOfRecord"] = "HUMAN_OBSERVATION"
    records1["GBIF_download_doi"] = "bypassed"

    # Add EBD records to a template dataframe
    schema_df = pd.DataFrame(columns=list(attribute_data_types.keys()))
    records2 = schema_df.combine_first(records1)

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  Results
    records2.to_csv(working_directory + "ebd_data.csv")
    print("Prepared the eBird records for the database: " + str(datetime.now() - timestamp))
    return records2

def get_GBIF_records(taxon_info, filter_set, query_name, working_directory, username, password, email):
    '''
    Retrieves species occurrence records from GBIF.  Filters occurrence
    records, buffers the xy points, and saves them in a database.  Finally,
    exports some Shapefiles.
    Gets species occurrence records from GBIF.  Can accomodate use of the GBIF
    API or Darwin Core Archive download via email.  Some filters can be applied
    during the query, but others have to be applied to the query results.

    Parameters
    ----------
    taxon_info : your taxon concept; dictionary
    filter_set : name of the filter set to apply; dictionary
    query_name : the name you chose for your query; string
    working_directory : path to use for table of filtered query results; string
    username : your GBIF username
    password : your GBIF password
    email : the email account associated with your GBIF account

    Returns
    -------
    Data frame of GBIF occurrence records
    '''
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

    # Some prep
    output_database = working_directory + query_name + '.sqlite'
    conn = sqlite3.connect(output_database, isolation_level='DEFERRED')
    cursor = conn.cursor()

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  TAXON INFO
    gbif_id = taxon_info["GBIF_ID"]
    taxon_polygon =taxon_info["TAXON_EOO"]

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  PREP FILTERS
    years = filter_set["years_range"]
    months = filter_set["months_range"]
    latRange = filter_set["lat_range"]
    lonRange = filter_set["lon_range"]
    coordinate = filter_set["has_coordinate_uncertainty"]
    geoIssue = filter_set["geoissue"]
    country = filter_set["country"]
    query_polygon = filter_set["query_polygon"]
    use_taxon_geometry = filter_set["use_taxon_geometry"]
    dwca_download = filter_set["get_dwca"]

    # Sort out geometries
    # A geometry could also be stated for the species, assess what to do
    # It could also be that user opted not to use species geometry.
    if use_taxon_geometry == False:
        taxon_polygon = None
    if query_polygon is None and taxon_polygon is None:
        poly = None
    elif query_polygon is not None and taxon_polygon is None:
        poly = query_polygon
    elif query_polygon is None and taxon_polygon is not None:
        poly = taxon_polygon
    elif query_polygon is not None and taxon_polygon is not None:
        # Get/use the intersection of the two polygons
        filter_polygon = shapely.wkt.loads(query_polygon)
        sp_polygon = shapely.wkt.loads(taxon_polygon)
        poly_intersection = filter_polygon.intersection(sp_polygon)
        poly = shapely.wkt.dumps(poly_intersection)

    print("Got request params and sorted out geometry constraints: " + str(datetime.now() - timestamp))
    timestamp = datetime.now()

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< GET RECORD COUNT
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
    record_count=occ_search["count"]
    print(str(record_count) + " records available")

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< API QUERY
    if dwca_download == False:
        # Get records in batches, saving into master list.
        all_jsons = []
        batches = range(0, record_count, 300)
        for i in batches:
            batch = occurrences.search(gbif_id,
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
            occs = batch['results']
            all_jsons = all_jsons + occs

        # Load records into a data frame
        dfRaw = pd.DataFrame(columns=attribute_data_types.keys())
        insertDict = {}
        for x in attribute_data_types.keys():
            insertDict[x] = []
        for x in all_jsons:
            present_keys = list(set(x.keys()) & set(attribute_data_types.keys()))
            for y in present_keys:
                insertDict[y] = insertDict[y] + [str(x[y])]
            missing_keys = list(set(attribute_data_types.keys()) - set(x.keys()))
            for z in missing_keys:
                insertDict[z] = insertDict[z] + ["UNKNOWN"]
        insertDF = pd.DataFrame(insertDict)
        df0 = dfRaw.append(insertDF, ignore_index=True, sort=False)
        df0copy = df0.copy() # a copy for gbif_fields_returned below

        df0.drop(["issue", "id"], inplace=True, axis=1)
        df0['coordinateUncertaintyInMeters'].replace(to_replace="UNKNOWN",
                                                     value=np.NaN, inplace=True)
        df0 = df0.astype({'coordinateUncertaintyInMeters': 'float',
                          'decimalLatitude': 'string', 'decimalLongitude': 'string'})
        df0['individualCount'].replace(to_replace="UNKNOWN", value=1,
                                       inplace=True)

        print("Downloaded records: " + str(datetime.now() - timestamp))
        timestamp = datetime.now()

        # Summarize the attributes that were returned
        #    Count entries per attribute(column), reformat as new df with appropriate
        #   columns.  Finally, insert into database.
        #    NOTE: When pulling from df0copy, only a specified subset of keys are
        #    assessed (attribute_data_types.keys()).  For a more complete picture, all_jsons must be
        #   assessed.  That has historically been very slow.
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
        timestamp = datetime.now()
        ####       START SLOW
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
        print("Summarized fields returned: " + str(datetime.now() - timestamp))
        #######       #####  END SLOW
        dfK.sort_index(inplace=True)
        dfK.index.name = 'attribute'

        # Save attribute summary into the output database
        dfK.to_sql(name='gbif_fields_returned', con=conn, if_exists='replace')

        # Prep the records for processing (filtering, storing in output db)
        schema = conn.execute("PRAGMA table_info (occurrence_records);").fetchall()
        column_names = [x[1] for x in schema]
        records0 = pd.DataFrame(columns=column_names)

        # Prep GBIF records
        records1 = (df0.rename({'id': 'record_id'}, axis=1)
                  .filter(records0.columns,axis=1))

        # Add GBIF records to template
        records2 = (records0
                    .combine_first(records1)
                    .fillna({"detection_distance_m": np.nan,
                             "radius_meters": np.nan,
                             "coordinateUncertaintyInMeters": np.nan,
                             "individualCount": np.nan})
                    .astype(attribute_data_types))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< EMAIL QUERY
    if dwca_download == True:
        # Make the data request using the download function.  Results are
        #   emailed.
        # First, build a query list.  NoneType values cause problems, so only
        #   add arguments if their value isn't NoneType.
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
        poly = None # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Remove and debug
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
        timestamp = datetime.now()
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
                zipdownload = occurrences.download_get(key=dkey,
                                                       path=working_directory)
                gotit = 1
                print("Download complete: " + str(datetime.now() - timestamp))
            except:
                wait = datetime.now() - timestamp
                if wait.seconds > 60*45:
                    gotit = 0
                    print("TIMED OUT -- attempting to proceed anyways")
                else:
                    gotit = None

        # Read the relevant files from within the darwin core archive
        timestamp = datetime.now()
        with DwCAReader(working_directory + dkey + '.zip') as dwca:
            dfRaw = dwca.pd_read('occurrence.txt', low_memory=False)
            citations = dwca.open_included_file('citations.txt').read()
            rights = dwca.open_included_file('rights.txt').read()
            doi = dwca.metadata.attrib["packageId"]

        # Record DWCA metadata
        #   Store the value summary for the selected fields in a table.
        timestamp = datetime.now()
        cursor.executescript("""CREATE TABLE GBIF_download_info
                                (download_key TEXT, doi TEXT, citations TEXT,
                                 rights TEXT);""")
        cursor.execute('''INSERT INTO GBIF_download_info (doi, citations,
                                                          rights, download_key)
                          VALUES ("{0}", "{1}", "{2}", "{3}")'''.format(doi,
                                                                  citations,
                                                                  rights,
                                                                  dkey))
        print("Stored GBIF Download DOI etc.: " + str(datetime.now() - timestamp))

        # Summarize the fields returned
        timestamp = datetime.now()
        df_populated1 = pd.DataFrame(dfRaw.count(axis=0).T.iloc[1:])
        df_populated1['included(n)'] = len(dfRaw)
        df_populated1['populated(n)'] = df_populated1[0]
        df_populated2 = df_populated1.filter(items=['included(n)', 'populated(n)'], axis='columns')
        df_populated2.index.name = 'attribute'
        df_populated2.to_sql(name='gbif_fields_returned', con=conn, if_exists='replace')
        print("Summarized fields returned: " + str(datetime.now() - timestamp))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  PREPARE
    timestep = datetime.now()
    # Rename columns
    records1 = dfRaw.rename({"issue": "issues", 'id': 'record_id'}, axis=1)
    # Drop columns
    records1 = records1.filter(items=attribute_data_types.keys(), axis=1)
    # Populate columns
    records1["retrieval_date"] = str(datetime.now())
    if filter_set['get_dwca'] == True:
        records1["GBIF_download_doi"] = doi
    else:
        records1["GBIF_download_doi"] = "bypassed"
    records1["source"] = "GBIF"
    print("Revised field names and values: " + str(datetime.now() - timestamp))


    #records1['coordinateUncertaintyInMeters'].replace(to_replace="UNKNOWN",
    #                                             value=np.NaN, inplace=True)
    #records1['decimalLatitude'] = records1['decimalLatitude'].astype(str)
    #records1['decimalLongitude'] = records1['decimalLongitude'].astype(str)
    #print(records1.individualCount.unique())
    #records1['individualCount'].replace(to_replace="UNKNOWN", value=1, inplace=True)

    # Add GBIF records to template
    records2 = (pd.DataFrame(columns=attribute_data_types.keys())
                .combine_first(records1))

    # Results
    records2.to_csv(working_directory + "gbif_data.csv")
    return records2

def filter_records(ebird_data, gbif_data, filter_set, taxon_info, working_directory, query_name):
    '''
    Summarizes the values in the data frames, apply filters, summarize what
        values persisted after filtering.  Insert results into the output db.

    PARAMETERS
    ----------
    ebird_data : a data frame of records from eBird
    gbif_data : a data frame of records from GBIF
    output_database : path to the output database
    filter_set : the filter set dictionary
    taxon_info : the taxon information dictionary

    RETURNS
    -------
    filtered_records : a data frame of filtered records.
    '''
    import sqlite3
    import pandas as pd
    from datetime import datetime
    import fnmatch

    # Create or connect to the database
    output_database = working_directory + query_name + ".sqlite"
    conn = sqlite3.connect(output_database, isolation_level='DEFERRED')
    cursor = conn.cursor()

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  MANAGE DATA TYPES
    schema = attribute_data_types
    string_atts = {key:value for (key, value) in schema.items() if schema[key] == 'str'}
    ebird_data = ebird_data.astype(string_atts)
    gbif_data = gbif_data.astype(string_atts)

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE EBIRD FROM GBIF
    if ebird_data is not None:
        gbif_data = gbif_data[gbif_data["collectionCode"].str.contains("EBIRD*") == False]

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  COMBINE DATA FRAMES
    if ebird_data is None:
        records3 = gbif_data
    if gbif_data is None:
        records3 = ebird_data
    if gbif_data is not None and ebird_data is not None:
        # Concatenate the gbif and ebird tables
        records3 = pd.concat([ebird_data, gbif_data])

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SUMMARIZE VALUES
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
                        'locationRemarks', 'occurrenceRemarks', 'radius_meters',
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

    # Store value summary in a dataframe
    acquired = summarize_values(dataframe=records3, step='acquired')

    # Summarize sources
    source_df1 = records3[['institutionID', 'collectionCode', 'datasetName', 'record_id']]
    source_summary1 = (source_df1
                       .groupby(by=['institutionID', 'collectionCode', 'datasetName'])
                       .size()
                       .reset_index(name='acquired'))
    print("Summarized values acquired: " + str(datetime.now() - timestamp))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  FILTER
    timestamp = datetime.now()

    # Some filters to be prepped for use
    for x in ['bases_omit', 'collection_codes_omit', 'datasets_omit',
              'institutions_omit', 'issues_omit', 'sampling_protocols_omit']:
        if filter_set[x] == None:
            filter_set[x] = []

    if filter_set['has_coordinate_uncertainty'] == 1:
        records5 = records3[pd.isnull(records3['coordinateUncertaintyInMeters']) == False]
    if filter_set['has_coordinate_uncertainty'] == 0:
        records5 = records3

    records6 = (records5[records5['coordinateUncertaintyInMeters'] <= filter_set['max_coordinate_uncertainty']]
                [lambda x: x['collectionCode'].isin(filter_set['collection_codes_omit']) == False]
                [lambda x: x['institutionID'].isin(filter_set['institutions_omit']) == False]
                [lambda x: x['basisOfRecord'].isin(filter_set['bases_omit']) == False]
                [lambda x: x['samplingProtocol'].isin(filter_set['sampling_protocols_omit']) == False]
                [lambda x: x['datasetName'].isin(filter_set['datasets_omit']) == False]
                )

    ''' ISSUES are more complex because multiple issues can be listed per record
    Method used is complex, but hopefully faster than simple iteration over all records
    '''
    records6.fillna(value={'issues': ""}, inplace=True)
    # Format of issues entries differ by method, change json format to email format
    if filter_set['get_dwca'] == True:
        records6['issues'] = [x.replace(', ', ';').replace('[', '').replace(']', '').replace("'", "")
                              for x in records6['issues']]
    unique_issue = list(records5['issues'].unique())
    violations = [x for x in unique_issue if len(set(str(x).split(";")) & set(filter_set['issues_omit'])) != 0] # entries that contain violations
    records7 = records6[records6['issues'].isin(violations) == False] # Records without entries that are violations.
    print("Performed filtering: " + str(datetime.now() - timestamp))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  POPULATE SOME COLUMNS
    if filter_set["default_coordUncertainty"] != None:
        records7.fillna(value={'coordinateUncertaintyInMeters': filter_set["default_coordUncertainty"]},
                        inplace=True)
    records7.fillna(value={'individualCount': int(1)}, inplace=True)
    records7["detection_distance_m"] = taxon_info["detection_distance_m"]
    records7["radius_meters"] = records7["detection_distance_m"] + records7["coordinateUncertaintyInMeters"]
    #records7.fillna({'radius_meters': 0})
    records7["weight"] = 10
    records7["weight_notes"] = ""
    records7["taxon_id"] = taxon_info["ID"]

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  REMOVE SPACE-TIME DUPLICATES
    # Prep some columns by changing data type
    records7 = (records7
                .astype({'decimalLatitude': 'str', 'decimalLongitude': 'str'})
                .reset_index(drop=True))

    if filter_set["duplicates_OK"] == False:
        records8 = drop_duplicates_latlongdate(records7)

    if filter_set["duplicates_OK"] == True:
        records8 = records7.copy()
        print("DUPLICATES ON LATITUDE, LONGITUDE, DATE-TIME INCLUDED")

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SPATIAL FILTERING
    # Spatial filtering happens in the get functions (ebird and gbif), not here.

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SUMMARIZE VALUES AGAIN
    # Store value summary in a dataframe
    retained = summarize_values(dataframe=records8, step='retained')

    # Concat acquired and retained data frames
    summary_df = pd.merge(retained, acquired, on=['attribute', 'value'], how='inner')

    # Calculate a difference column
    summary_df['removed'] = summary_df['acquired'] - summary_df['retained']
    summary_df = summary_df[['attribute', 'value', 'acquired', 'removed',
                             'retained']]

    # Summarize sources
    source_df2 = records8[['institutionID', 'collectionCode', 'datasetName', 'record_id']]
    source_summary2 = (source_df2
                       .groupby(by=['institutionID', 'collectionCode', 'datasetName'])
                       .size()
                       .reset_index(name='retained'))

    # Concat acquired and retained source summary data frames
    source_summaries = pd.merge(source_summary1, source_summary2,
                                on=['institutionID', 'collectionCode', 'datasetName'],
                                how='inner')
    print(source_summaries)

    # Calculate a difference column
    source_summaries['removed'] = source_summaries['acquired'] - source_summaries['retained']
    source_summaries = source_summaries[['institutionID', 'collectionCode',
                                         'datasetName', 'acquired', 'removed',
                                         'retained']]

    # Save the summaries in the output database
    summary_df.to_sql(name='attribute_value_counts', con = conn, if_exists='replace')
    source_summaries.to_sql(name='sources', con = conn, if_exists='replace')

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SAVE
    # Reformat data to strings and insert into db.
    records7.applymap(str).to_sql(name='occurrence_records', con = conn,
                                  if_exists='replace')
    conn.close()
    return None

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
    df['dup_latPlaces'] = [len(x.split(".")[1]) for x in df['decimalLatitude']]
    df['dup_lonPlaces'] = [len(x.split(".")[1]) for x in df['decimalLongitude']]
    df['dup_OGprec'] = df['dup_latPlaces']
    prec_unequal = df[df['dup_latPlaces'] != df['dup_lonPlaces']]
    for i in prec_unequal.index:
        x = prec_unequal.loc[i]
        if x['dup_latPlaces'] < x['dup_lonPlaces']:
            trim_len = int(x['dup_latPlaces'])
        else:
            trim_len = int(x['dup_lonPlaces'])
        df.loc[i, 'decimalLatitude'] = x['decimalLatitude'][:trim_len + 3]
        df.loc[i, 'decimalLongitude'] = x['decimalLongitude'][:trim_len + 4]
        # Record the resulting precision for reference later
        df.loc[i, 'dup_OGprec'] = trim_len
    df.drop(['dup_latPlaces', 'dup_lonPlaces'], axis=1, inplace=True)

    """
    ########  INITIAL DROP OF DUPLICATES
    Initial drop of duplicates on 'latitude', 'longitude', 'eventDate',
    keeping the first (highest individual count)
    Sort so that the highest individual count is first ## ADD OCCURRENCEDATE BACK IN
    """
    df.sort_values(by=['decimalLatitude', 'decimalLongitude', 'eventDate',
                        'individualCount'],
                    ascending=False, inplace=True, kind='mergesort',
                    na_position='last')

    df.drop_duplicates(subset=['decimalLatitude', 'decimalLongitude', 'eventDate'],
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
        precision : The level of precision (places right of decimal) in decimalLatitude
        and longitude values for the assessment of duplicates.
        df : dataframe to assess and drop duplicates from.  This function works
              'inplace'.
        """
        # Create a df with records from the input df having decimal precision > the
        # precision level being assessed.
        dfLonger = df[df['dup_OGprec'] > precision].copy()
        # Truncate lat and long values
        dfLonger['decimalLatitude'] = [x[:precision + 3] for x in dfLonger['decimalLatitude']]
        dfLonger['decimalLongitude'] = [x[:precision + 4] for x in dfLonger['decimalLongitude']]

        # Create a df with records having the precision being
        # investigated
        dfShorter1 = df[df['dup_OGprec'] == precision]

        # Find records in dfShorter1 with latitude, longitude, date combo
        # existing in dfLonger and append to list of duplis
        dfduplis = pd.merge(dfShorter1, dfLonger, how='inner',
                            on=['decimalLatitude', 'decimalLongitude', 'eventDate'])
        dups_ids = dfduplis['record_id_x']
        for d in dups_ids:
            duplis.append(d)

    # Drop latitude longitude duplicates at lower decimal precisions
    for p in precisions:
        drop_duplicates(p, df)

    # Drop rows from the current main df that have been identified as duplicates.
    df2 = df[df['record_id'].isin(duplis) == False].copy()

    # Drop excess columns
    df2.drop(['dup_OGprec'], inplace=True, axis=1)

    duptime = datetime.now() - startduptime
    print(str(initial_length - len(df2)) + " duplicate records dropped: {0}".format(duptime))
    return df2

# Helper functions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def get_GBIF_code(name, rank='species'):
    """
    Returns the GBIF species code for a scientific name.

    Example: gbifcode = getGBIFcode(name = "Dendroica cerulea")
    """
    from pygbif import species
    key = species.name_backbone(name = name, rank='species')['usageKey']
    return key

def get_record_details(key):
    """
    Returns a dictionary holding all GBIF details about the record.

    Example: details = getRecordDetails(key = 1265907957)
    """
    from pygbif import occurrences
    details = occurrences.get(key = key)
    return details

def generate_shapefile(database, output_file, footprints=True):
    '''
    Exports a shapefile of species occurrence records from a wildlife wrangler
        output SQLite database.

    PARAMETERS
    database : the sqlite database to use.
    output_file : Path (and name) of the file to be created.
    footprints : True creates a map of record footprints (buffered coordinates),
        whereas False returns records as points base on coordinates provided.

    OUTPUT
    shapefile saved to disk
    '''
    from datetime import datetime
    import sqlite3
    import pandas as pd
    import geopandas as gpd
    timestamp = datetime.now()

    # Get the record coordinates as a data frame
    records = (pd.read_sql("""SELECT * FROM occurrence_records;""",
                           con=sqlite3.connect(database))
                          .astype({'decimalLongitude': 'float',
                                   'decimalLatitude': 'float',
                                   'radius_meters': 'float'}))

    # Make a point geometry from coordinates
    gdf = gpd.GeoDataFrame(records,
                           geometry=gpd.points_from_xy(records['decimalLongitude'],
                                                       records['decimalLatitude']))

    # Set the coordinate reference system
    gdf.crs={'init' :'epsg:4326'}

    # Save points OR
    if footprints == False:
        gdf.to_file(outFile)

    # Save polygons (footprints)
    if footprints == True:
        # Reproject states and record coordinates to facilitate buffering
        footprints = gdf.to_crs(epsg=5070)

        # Buffer points with appropriate radii for record footprints
        footprints['footprint']=footprints.apply(lambda x: x.geometry.buffer(x.radius_meters), axis=1)
        footprints.set_geometry(col='footprint', inplace=True, drop=True)
        footprints.to_file(outFile)

    print("Exported shapefile: " + str(datetime.now() - timestamp))
    return

def ccw_wkt_from_shapefile(shapefile, out_txt):
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
