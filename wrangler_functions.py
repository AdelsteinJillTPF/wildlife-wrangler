# occurrece records table datatypes
output_schema = {"GBIF_download_doi": "str",
                 "accessRights": "str",
                 "basisOfRecord": "str",
                 "bibliographicCitation": "str",
                 "collectionCode": "str",
                 "coordinatePrecision": "float",
                 "coordinateUncertaintyInMeters": "float",
                 "dataGeneralizations": "str",
                 "datasetName": "str",
                 "decimalLatitude": "str",
                 "decimalLongitude": "str",
                 "detection_distance_m": "int",
                 "ebird_id": "str",
                 "effort_distance_m": "int",
                 "establishmentMeans": "str",
                 "eventDate": "str",
                 "eventRemarks": "str",
                 "filter_set_name": "str",
                 "footprintSRS": "str",
                 "footprintWKT": "str",
                 "gbif_id": "str",
                 "general_remarks": "str",
                 "geodeticDatum": "str",
                 "georeferenceProtocol": "str",
                 "georeferenceRemarks": "str",
                 "georeferenceVerificationStatus": "str",
                 "georeferencedBy": "str",
                 "gps_accuracy_m": "int",
                 "habitat": "str",
                 "identificationQualifier": "str",
                 "identifiedBy": "str",
                 "identifiedRemarks": "str",
                 "individualCount": "int",
                 "informationWitheld": "str",
                 "institutionID": "str",
                 "issues": "str",
                 "license": "str",
                 "locality": "str",
                 "locationAccordingTo": "str",
                 "locationRemarks": "str",
                 "modified": "str",
                 "nominal_xy_precision": "float",
                 "occurrenceRemarks": "str",
                 "occurrenceStatus": "str",
                 "organismQuantity": "str",
                 "organismQuantityType": "str",
                 "radius_m": "float",
                 "record_id": "str",
                 "recordedBy": "str",
                 "retrieval_date": "str",
                 "samplingProtocol": "str",
                 "samplingEffort": "str",
                 "scientificName": "str",
                 "source": "str",
                 "taxonConceptID": "str",
                 "taxon_info_name": "str",
                 "verbatimLocality": "str",
                 "weight": "int",
                 "weight_notes": "str"}

# Core functions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def build_output_database(output_database):
    """
    Create a database for storing occurrence and taxon concept data.
    The column names that are "camel case" are Darwin Core attributes, whereas
    lower case names containing "_" between words are not.  Only those Darwin
    core attributes that could be useful for filtering and assigning weights
    are included.

    Parameters
    ----------
    output_database : Path for sqlite database to create; string.

    Returns
    -------
    Nothing.
    """
    import os
    import sqlite3
    import pandas as pd

    # Delete the database if it already exists
    if os.path.exists(output_database):
        os.remove(output_database)

    # Create or connect to the database
    conn = sqlite3.connect(output_database)

    # Create a table for occurrence records.
    df = (pd.DataFrame(columns=output_schema.keys())
            .astype(output_schema)
            .to_sql(name='occurrence_records', con=conn, if_exists='replace'))

    conn.close()
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
    import pandas as pd
    from datetime import datetime
    import sqlite3
    import geopandas as gpd
    import shapely
    from shapely.wkt import dumps, loads
    import numpy as np
    import rpy2.robjects as robjects
    import rpy2.robjects.packages as rpackages
    from rpy2.robjects.vectors import StrVector

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
    for x in taxon_info.keys():
        if taxon_info[x] == None:
            taxon_info[x] = ""

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
    #   This all requires formatting dates as text correctly.
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
    print(ok_months)

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
    max_distance <- as.integer(ceiling((max_coordinate_uncertainty-(2*EBD_gps_precision))/1000))  ### NOT WORKING????
    print(max_distance)

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

    # It could be that user opted not to use species geometry.
    if filter_set['use_taxon_geometry'] == False:
        EOO = None

    # A geometry could be stated for the species, assess what to do
    # Replace "" values in fitler_set with None to fit Python code.
    for x in filter_set.keys():
        if filter_set[x] == "":
            filter_set[x] = None
    for x in taxon_info.keys():
        if taxon_info[x] == "":
            taxon_info[x] = None

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
    records1 = gdf.filter(list(output_schema.keys()), axis=1)
    # Populate columns
    records1["institutionID"] = "clo"
    records1["collectionCode"] = "EBIRD"
    records1["datasetName"] = "EBD"
    records1["source"] = "eBird"
    records1["basisOfRecord"] = "HUMAN_OBSERVATION"
    records1["GBIF_download_doi"] = "bypassed"
    records1["occurrenceStatus"] = "PRESENT"
    records1 = (records1
                .fillna({"effort_distance_m": 0, "gps_accuracy_m": 10})
                .replace({"individualCount": {"X": 1}}))

    # Add EBD records to a template dataframe
    schema_df = pd.DataFrame(columns=list(output_schema.keys()))
    records2 = schema_df.combine_first(records1)

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  Results
    print("Prepared the eBird records for processing: " + str(datetime.now() - timestamp))
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
    EOO = taxon_info["TAXON_EOO"]
    AOI = filter_set["query_polygon"]

    # It could be that user opted not to use species geometry.
    if filter_set['use_taxon_geometry'] == False:
        EOO = None

    # A geometry could be stated for the species, assess what to do
    if AOI is None and EOO is None:
        filter_polygon = None
    elif AOI is not None and EOO is None:
        filter_polygon = AOI
    elif AOI is None and EOO is not None:
        filter_polygon = EOO
    elif AOI is not None and EOO is not None:
        # Get/use the intersection of the two polygons in this case
        AOI_polygon = shapely.wkt.loads(AOI)
        EOO_polygon = shapely.wkt.loads(EOO)
        intersection = AOI_polygon.intersection(EOO_polygon)
        # Make the polygon's outer ring counter clockwise
        if intersection.exterior.is_ccw == False:
            print("Reordered filter polygon coordinates")
            intersection = shapely.geometry.polygon.orient(intersection, sign=1.0)
            # Get the well-known text version of the polygon
            filter_polygon = shapely.wkt.dumps(intersection)
        else:
            filter_polygon = shapely.wkt.dumps(intersection)

    #print(shapely.wkt.loads(filter_polygon).exterior.is_ccw)
    print("Prepared filter set and sorted out geometry constraints: " + str(datetime.now() - timestamp))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< GET RECORD COUNT
    timestamp = datetime.now()
    # First, find out how many records there are that meet criteria
    occ_search = occurrences.search(gbif_id,
                                    year=years,
                                    month=months,
                                    decimalLatitude=latRange,
                                    decimalLongitude=lonRange,
                                    hasGeospatialIssue=geoIssue,
                                    hasCoordinate=True,
                                    country=country,
                                    geometry=filter_polygon)
    record_count=occ_search["count"]
    print(str(record_count) + " records available")

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< API QUERY
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
                                       hasCoordinate=True,
                                       country=country,
                                       geometry=filter_polygon)
            occs = batch['results']
            all_jsons = all_jsons + occs

        # Get a list of keys that were returned
        api_keys = set([])
        for j in all_jsons:
            api_keys = api_keys | set(j.keys())

        # Load json records into a dataframe, via a dictionary
        insertDict = {}
        for k in list(api_keys):
            insertDict[k] = []
        for j in all_jsons:
            present_keys = list(set(j.keys()) & api_keys)
            for prk in present_keys:
                insertDict[prk] = insertDict[prk] + [str(j[prk])]
            missing_keys = list(api_keys - set(j.keys()))
            for mik in missing_keys:
                insertDict[mik] = insertDict[mik] + ["UNKNOWN"]
        dfRaw = pd.DataFrame(insertDict).rename({"occurrenceID": "record_id"}, axis=1)
        print("Downloaded records: " + str(datetime.now() - timestamp))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< EMAIL QUERY
    if dwca_download == True:
        timestamp = datetime.now()
        '''
        Request data using the download function.  Results are emailed.
        '''
        # First, build a query list.  NoneType values cause problems, so only
        #   add arguments if their value isn't NoneType.
        download_filters = ['taxonKey = {0}'.format(gbif_id)]
        download_filters.append('hasCoordinate = True')
        if country is not None:
            download_filters.append('country = {0}'.format(country))
        if years is not None:
            download_filters.append('year >= {0}'.format(years.split(",")[0]))
            download_filters.append('year <= {0}'.format(years.split(",")[1]))
        if months is not None:
            download_filters.append('month >= {0}'.format(months.split(",")[0]))
            download_filters.append('month <= {0}'.format(months.split(",")[1]))
        if filter_polygon is not None:
            download_filters.append("geometry = {0}".format(filter_polygon))
        if geoIssue is not None:
            download_filters.append('hasGeospatialIssue = {0}'.format(geoIssue))
        if latRange is not None:
            download_filters.append('decimalLatitude >= {0}'.format(latRange.split(",")[0]))
            download_filters.append('decimalLatitude <= {0}'.format(latRange.split(",")[1]))
        if lonRange  is not None:
            download_filters.append('decimalLongitude >= {0}'.format(lonRange.split(",")[0]))
            download_filters.append('decimalLongitude <= {0}'.format(lonRange.split(",")[1]))

        # Get the value of the download key
        try:
            d = occurrences.download(download_filters,
                                     pred_type='and',
                                     user = username,
                                     pwd = password,
                                     email = email)
            dkey = d[0]
        except Exception as e:
            print(e)
            print(download_filters)

        # Now download the actual zip file containing the Darwin Core files
        # NOTE: The download can take a while to generate and is not immediately
        # available once the download_get command has been issued. Use a
        # while and try loop to make sure the download has succeeded.
        # The zipdownload variable will be a dictionary of the path,
        # the file size, and the download key unique code. It can be used
        # to change the file name, unzip the file, etc.
        print("Downloading Darwin Core Archive zip file for this species .....")
        timestamp = datetime.now()
        gotit = False
        while gotit == False:
            try:
                zipdownload = occurrences.download_get(key=dkey,
                                                       path=working_directory)
            except:
                pass

            try:
                # Read the relevant files from within the darwin core archive
                with DwCAReader(working_directory + dkey + '.zip') as dwca:
                    dfRaw = dwca.pd_read('occurrence.txt', low_memory=False)
                    doi = dwca.metadata.attrib["packageId"]
                    try:
                        citations = dwca.open_included_file('citations.txt').read()
                    except Exception as e:
                        citations = "Failed"
                        print(e)
                    try:
                        rights = dwca.open_included_file('rights.txt').read()
                    except Exception as e:
                        rights = "Failed"
                        print(e)
                    gotit = True
                print("Download complete: " + str(datetime.now() - timestamp))
            except:
                wait = datetime.now() - timestamp
                if wait.seconds > 60*45:
                    gotit = True
                    print("TIMED OUT")
                else:
                    pass

        # Record DWCA metadata
        #   Store the value summary for the selected fields in a table.
        timestamp = datetime.now()
        cursor.executescript("""CREATE TABLE GBIF_download_info
                                (download_key TEXT, doi TEXT, citations TEXT,
                                 rights TEXT);""")
        cursor.execute('''INSERT INTO GBIF_download_info (doi, download_key)
                          VALUES ("{0}", "{1}")'''.format(doi, dkey))

        try:
            cursor.execute('''INSERT INTO GBIF_download_info (citations)
                              VALUES ("{0}")'''.format(citations))
        except Exception as e:
            print(e)
            cursor.execute('''INSERT INTO GBIF_download_info (citations)
                              VALUES ("Failed")''')
        try:
            cursor.execute('''INSERT INTO GBIF_download_info (rights)
                              VALUES ("{0}")'''.format(rights))
        except Exception as e:
            print(e)
            cursor.execute('''INSERT INTO GBIF_download_info (rights)
                              VALUES ("Failed")''')

        print("Stored GBIF Download DOI etc.: " + str(datetime.now() - timestamp))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SUMMARIZE
    timestamp = datetime.now()
    if dwca_download == False: # We don't want to count the "UNKNOWNS" we added.
        df_raw2 = dfRaw.replace({"UNKNOWN": np.nan})
        df_populated1 = pd.DataFrame(df_raw2.count(axis=0).T.iloc[1:])
        df_populated1['included(n)'] = df_populated1[0]
        df_populated1['populated(n)'] = df_populated1[0]

    if dwca_download == True:
        df_raw2 = dfRaw.copy()
        df_populated1 = pd.DataFrame(df_raw2.count(axis=0).T.iloc[1:])
        df_populated1['included(n)'] = len(dfRaw)
        df_populated1['populated(n)'] = df_populated1[0]
    df_populated2 = df_populated1.filter(items=['included(n)', 'populated(n)'],
                                         axis='columns')
    df_populated2.index.name = 'attribute'
    df_populated2.to_sql(name='gbif_fields_returned', con=conn, if_exists='replace')
    print("Summarized fields returned: " + str(datetime.now() - timestamp))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  PREPARE
    timestep = datetime.now()
    # Rename columns
    records1 = dfRaw.rename({"issue": "issues", "id": "record_id"}, axis=1)
    # Drop columns
    records1 = records1.filter(items=output_schema.keys(), axis=1)
    # Populate columns
    records1["retrieval_date"] = str(datetime.now())
    if filter_set["get_dwca"] == True:
        records1["GBIF_download_doi"] = doi
    else:
        records1["GBIF_download_doi"] = "bypassed"
    records1["source"] = "GBIF"

    # Add GBIF records to template; replace and fillna to support astype()
    records2 = (pd.DataFrame(columns=output_schema.keys())
                .combine_first(records1)
                # this replace is needed for API method
                .replace({"coordinateUncertaintyInMeters": {"UNKNOWN": np.nan},
                          "radius_m": {"UNKNOWN": np.nan},
                          "coordinatePrecision": {"UNKNOWN": np.nan},
                          "nominal_xy_precision": {"UNKNOWN": np.nan},
                          "individualCount": {"UNKNOWN": 1},
                          "weight": {"UNKNOWN": 10},
                          "detection_distance_m": {"UNKNOWN": 0}})
                .fillna({"coordinateUncertaintyInMeters": 0,
                         "radius_m": 0,
                         "individualCount": 1,
                         "weight": 10,
                         "detection_distance_m": 0,
                         "effort_distance_m": 0,
                         "coordinate_precision": 1,
                         "gps_accuracy_m": 30})
                .astype(output_schema))
    print("Prepared GBIF records for processing: " + str(datetime.now() - timestamp))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  Results
    return records2

def process_records(ebird_data, gbif_data, filter_set, taxon_info, working_directory, query_name):
    '''
    Summarizes the values in the data frames, populates some fields,
        apply filters, summarize what values persisted after filtering.  Insert
        results into the output db.

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
    import numpy as np
    from datetime import datetime
    timestamp = datetime.now()

    # Create or connect to the database
    output_database = working_directory + query_name + ".sqlite"
    conn = sqlite3.connect(output_database, isolation_level='DEFERRED')
    cursor = conn.cursor()

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  MANAGE DATA TYPES
    schema = output_schema
    string_atts = {key:value for (key, value) in schema.items() if schema[key] == 'str'}

    if ebird_data is not None:
        ebird_data = ebird_data.astype(string_atts)
    if gbif_data is not None:
        gbif_data = gbif_data.astype(string_atts)

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE EBIRD FROM GBIF
    if gbif_data is not None:
        if ebird_data is not None:
            gbif_data = gbif_data[gbif_data["collectionCode"].str.contains("EBIRD*") == False]

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  COMBINE DATA FRAMES
    if ebird_data is None:
        df_unfiltered = gbif_data
    if gbif_data is None:
        df_unfiltered = ebird_data
    if gbif_data is not None and ebird_data is not None:
        # Concatenate the gbif and ebird tables
        df_unfiltered = pd.concat([ebird_data, gbif_data])
    print("Prepared data frames for processing: " + str(datetime.now() - timestamp))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SUMMARIZE VALUES
    timestamp = datetime.now()
    # Make a list of columns to summarize values from
    do_not_summarize = ['decimalLatitude', 'decimalLongitude', 'GBIF_download_doi',
                        'coordinateUncertaintyInMeters', 'detection_distance_m',
                        'eventDate', 'eventRemarks', 'filter_set_name',
                        'footprintSRS', 'footprintWKT', 'gbif_id', 'ebird_id',
                        "effort_distance_m", 'general_remarks', 'georeferencedBy',
                        'habitat', 'georeferenceRemarks', 'identificationQualifier',
                        'identifiedBy', 'identifiedRemarks', 'individualCount',
                        'informationWitheld', 'locality', 'locationAccordingTo',
                        'locationRemarks', "modified", 'occurrenceRemarks',
                        'radius_m', 'record_id', 'recordedBy', 'retrieval_date',
                        'taxonConceptID', 'verbatimLocality', 'weight',
                        'weight_notes']

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
            value_df["attribute"] = column
            value_df = value_df[["attribute", "value", step]]
            if value_df.empty == False:
                attributes.append(value_df)
        result = pd.concat(attributes)
        return result

    # Store value summary in a dataframe
    acquired = summarize_values(dataframe=df_unfiltered, step='acquired')

    # Summarize sources
    source_df1 = df_unfiltered[['institutionID', 'collectionCode', 'datasetName', 'record_id']]
    source_summary1 = (source_df1
                       .groupby(by=['institutionID', 'collectionCode', 'datasetName'])
                       .size()
                       .reset_index(name='acquired'))
    print("Summarized values acquired: " + str(datetime.now() - timestamp))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  POPULATE SOME COLUMNS
    timestamp = datetime.now()
    df_unfiltered.fillna(value={'individualCount': int(1)}, inplace=True)
    df_unfiltered["weight"] = 10
    df_unfiltered["weight_notes"] = ""
    df_unfiltered["taxon_id"] = taxon_info["ID"]
    df_unfiltered["gbif_id"] = taxon_info["GBIF_ID"]
    df_unfiltered["ebird_id"] = taxon_info["EBIRD_ID"]
    df_unfiltered["detection_distance_m"] = taxon_info["detection_distance_m"]
    df_unfiltered["filter_set_name"] = filter_set["name"]

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  COORDINATE PRECISION
    '''In WGS84, coordinate precision is limited by longitude and varies across
    latitudes and number of digits provided.  Thus, coordinates have a nominal
    precision that may limit values.   Populate a column for this...'''

    # Trim decimal length to 5 digits (lat and long).  Anything more is false precision.
    df_unfiltered["decimalLatitude"] = df_unfiltered["decimalLatitude"].apply(lambda x: x.split(".")[0] + "." + x.split(".")[1][:5])
    df_unfiltered["decimalLongitude"] = df_unfiltered["decimalLongitude"].apply(lambda x: x.split(".")[0] + "." + x.split(".")[1][:5])

    # Calculate the number of digits for latitude and longitude
    df_unfiltered['digits_latitude'] = [len(x.split(".")[1]) for x in df_unfiltered['decimalLatitude']]
    df_unfiltered['digits_longitude'] = [len(x.split(".")[1]) for x in df_unfiltered['decimalLongitude']]

    # Longitude precision
    digitsX = {1: 10, 2: 100, 3: 1000, 4: 10000, 5: 100000}
    df_unfiltered["temp"] = df_unfiltered["decimalLatitude"].apply(lambda x: 111321 * np.cos(float(x) * np.pi/180))
    df_unfiltered["temp2"] = df_unfiltered["digits_longitude"].apply(lambda x: digitsX[x])
    df_unfiltered["nominal_x_precision"] = 100 * df_unfiltered["temp"]/df_unfiltered["temp2"]# decimal moved based on digits.
    # Latitude precision
    digitsY = {1: 11112.0, 2: 1111.2, 3: 111.1, 4: 11.1, 5: 1.1} # Lookup for latitude precision
    df_unfiltered["nominal_y_precision"] = df_unfiltered["digits_latitude"].apply(lambda x: digitsY[x])

    # Put the larger of the two nominal precisions in a column
    df_unfiltered["nominal_xy_precision"] = np.where(df_unfiltered["nominal_y_precision"] > df_unfiltered["nominal_x_precision"], df_unfiltered["nominal_y_precision"], df_unfiltered["nominal_x_precision"])

    # Clean up
    df_unfiltered.drop(["temp", "temp2", "digits_latitude", "digits_longitude", "nominal_x_precision", "nominal_y_precision"], axis=1, inplace=True)

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< BUFFER RADIUS
    '''
    Calculate a buffer distance from various parameters for the
    point-radius method.  Compilation has to differ with data source and
    whether the user chose to use a default coordinate uncertainty.  Components
    of radius may include coordinateUncertaintyInMeters, coordinatePrecision,
    GPS_accuracy_m, effort_distance_m, detection_distance_m.

    Records are broken apart by source (GBIF, GBIF/EOD, EBD), processed,
    and then concatenated in order to account for all conditions.

    If footprintWKT is provided, it will be used by spatial_output instead
    of point buffering.
    '''
    # Records from GBIF with coordinate uncertainty (georeferenced)
    georef = df_unfiltered[df_unfiltered["coordinateUncertaintyInMeters"] > 0.0].copy()
    if georef.empty == False:
        georef.fillna({"coordinatePrecision":1.1}, inplace=True)
        georef["gps_accuracy_m"] = np.where(georef["eventDate"].apply(lambda x: datetime.strptime(x, "%Y-%m-%dT%H:%M:%S").year) < 2000, 100, 30)
        georef["radius_m"] = georef["coordinateUncertaintyInMeters"]
        print("Number of georeferenced GBIF records: " + str(len(georef)))

    # Records from GBIF without coordinate uncertainty
    gbif_nogeo = df_unfiltered[(df_unfiltered["coordinateUncertaintyInMeters"] == 0.0) & (df_unfiltered["collectionCode"].str.contains("EBIRD*") == False)].copy()
    if gbif_nogeo.empty == False:
        gbif_nogeo["gps_accuracy_m"] = np.where(gbif_nogeo["eventDate"].apply(lambda x: datetime.strptime(x, "%Y-%m-%dT%H:%M:%S").year) < 2000, 100, 30)
        if filter_set["default_coordUncertainty"] is not None:
            print("Applying default coordinate uncertainties for GBIF records")
            gbif_nogeo.fillna({"coordinatePrecision":1.1}, inplace=True)
            gbif_nogeo["radius_m"] = filter_set["default_coordUncertainty"]

        if filter_set["default_coordUncertainty"] is None:
            print("Approximating coordinate uncertanties for GBIF records")
            gbif_nogeo.fillna({"coordinatePrecision":1.1}, inplace=True)
            gbif_nogeo["radius_m"] = gbif_nogeo["gps_accuracy_m"] + gbif_nogeo["detection_distance_m"] + gbif_nogeo["effort_distance_m"]

    # Records from EBD
    ebd_geo = df_unfiltered[df_unfiltered["source"] == "eBird"].copy()
    if ebd_geo.empty == False:
        ebd_geo.fillna({"coordinatePrecision":1.1}, inplace=True)
        ebd_geo["gps_accuracy_m"] = np.where(ebd_geo["eventDate"].apply(lambda x: datetime.strptime(x, "%Y-%m-%d").year) < 2000, 100, 30)
        ebd_geo["radius_m"] = ebd_geo["effort_distance_m"] + ebd_geo["gps_accuracy_m"] + ebd_geo["detection_distance_m"]

    # Records from EOD (via GBIF)
    eod_nogeo = df_unfiltered[(df_unfiltered["source"] == "GBIF") & (df_unfiltered["collectionCode"].str.contains("EBIRD*") == True)].copy()
    if eod_nogeo.empty == False:
        eod_nogeo.fillna({"coordinatePrecision":1.1}, inplace=True)
        eod_nogeo["gps_accuracy_m"] = np.where(eod_nogeo["eventDate"].apply(lambda x: datetime.strptime(x, "%Y-%m-%dT%H:%M:%S").year) < 2000, 100, 30)
        eod_nogeo["effort_distance_m"] = 8047 # eBird best practices allows distance up to 5 mi length.
        eod_nogeo["radius_m"] = eod_nogeo["effort_distance_m"] + eod_nogeo["gps_accuracy_m"] + eod_nogeo["detection_distance_m"]

    # Concat df's if necessary
    if filter_set['has_coordinate_uncertainty'] == True:
        df_unfiltered2 = georef

    to_concat = []
    for x in [gbif_nogeo, georef, eod_nogeo, ebd_geo]:
        if x.empty == False:
            to_concat.append(x)

    if len(to_concat) > 1:
        df_unfiltered2 = pd.concat(to_concat)
    if len(to_concat) == 1:
        df_unfiltered2 = to_concat[0]

    # Where coordinate precision is poor, overwrite the radius to be the precision.
    df_unfiltered2["radius_m"] = np.where(df_unfiltered2["nominal_xy_precision"] > df_unfiltered2["radius_m"], df_unfiltered2["nominal_xy_precision"], df_unfiltered2["radius_m"])
    df_unfiltered2["radius_m"] = np.where(df_unfiltered2["coordinatePrecision"] > df_unfiltered2["radius_m"], df_unfiltered2["coordinatePrecision"], df_unfiltered2["radius_m"])

    # Test to make sure that no records were lost in the previous steps
    if len(df_unfiltered2) != len(df_unfiltered):
        print("AN ERROR OCCURRED !!!!!!!!!!!!!")
    else:
        print("Prepared records and calculated radii:" + str(datetime.now() - timestamp))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  FILTER
    timestamp = datetime.now()
    # Some filters to be prepped for use
    for x in ['bases_omit', 'collection_codes_omit', 'datasets_omit',
              'institutions_omit', 'issues_omit', 'sampling_protocols_omit']:
        if filter_set[x] == None:
            filter_set[x] = []

    df_filter2 = (df_unfiltered2[df_unfiltered2['radius_m'] <= filter_set['max_coordinate_uncertainty']]
                [lambda x: x['collectionCode'].isin(filter_set['collection_codes_omit']) == False]
                [lambda x: x['institutionID'].isin(filter_set['institutions_omit']) == False]
                [lambda x: x['basisOfRecord'].isin(filter_set['bases_omit']) == False]
                [lambda x: x['samplingProtocol'].isin(filter_set['sampling_protocols_omit']) == False]
                [lambda x: x['datasetName'].isin(filter_set['datasets_omit']) == False]
                [lambda x: x['occurrenceStatus'] != "ABSENT"]
                )

    # Case where user demands records had coordinate uncertainty
    if filter_set['has_coordinate_uncertainty'] == True:
        df_filter2 = df_filter2[df_filter2["coordinateUncertaintyInMeters"] > 0]

    ''' ISSUES are more complex because multiple issues can be listed per record
    Method used is complex, but hopefully faster than simple iteration over all records
    '''
    df_filter2.fillna(value={'issues': ""}, inplace=True)
    # Format of issues entries differ by method, change json format to email format
    if filter_set['get_dwca'] == True:
        df_filter2['issues'] = [x.replace(', ', ';').replace('[', '').replace(']', '').replace("'", "")
                              for x in df_filter2['issues']]
    unique_issue = list(df_filter2['issues'].unique())
    violations = [x for x in unique_issue if len(set(str(x).split(";")) & set(filter_set['issues_omit'])) != 0] # entries that contain violations
    df_filter3 = df_filter2[df_filter2['issues'].isin(violations) == False] # Records without entries that are violations.
    print("Performed filtering: " + str(datetime.now() - timestamp))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  REMOVE SPACE-TIME DUPLICATES
    # Prep some columns by changing data type
    df_filter3 = (df_filter3
                  .astype({'decimalLatitude': 'str', 'decimalLongitude': 'str'})
                  .reset_index(drop=True))
    if filter_set["duplicates_OK"] == False:
        df_filterZ = drop_duplicates_latlongdate(df_filter3)

    if filter_set["duplicates_OK"] == True:
        df_filterZ = df_filter3.copy()
        print("DUPLICATES ON LATITUDE, LONGITUDE, DATE-TIME INCLUDED")

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SPATIAL FILTERING
    # Spatial filtering happens in the get functions (ebird and gbif), not here.

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SUMMARIZE VALUES AGAIN
    timestamp = datetime.now()
    # Store value summary in a dataframe
    if df_filterZ.empty == False:
        retained = summarize_values(dataframe=df_filterZ, step='retained')

    if df_filterZ.empty == True:
        retained = acquired.copy().drop(["acquired"], axis=1)
        retained["retained"] = 0

    # Concat acquired and retained data frames
    summary_df = pd.merge(retained, acquired, on=['attribute', 'value'], how='inner')

    # Calculate a difference column
    summary_df['removed'] = summary_df['acquired'] - summary_df['retained']
    summary_df = summary_df[['attribute', 'value', 'acquired', 'removed',
                             'retained']]

    # Summarize sources
    if df_filterZ.empty == False:
        source_df2 = df_filterZ[['institutionID', 'collectionCode', 'datasetName', 'record_id']]
        source_summary2 = (source_df2
                           .groupby(by=['institutionID', 'collectionCode', 'datasetName'])
                           .size()
                           .reset_index(name='retained'))

    if df_filterZ.empty == True:
        print(source_summary1)
        source_summary2 = source_summary1.copy().drop(["acquired"], axis=1)
        source_summary2["retained"] = 0

    # Concat acquired and retained source summary data frames
    source_summaries = pd.merge(source_summary1, source_summary2,
                                on=['institutionID', 'collectionCode', 'datasetName'],
                                how='inner')

    # Calculate a difference column
    source_summaries['removed'] = source_summaries['acquired'] - source_summaries['retained']
    source_summaries = source_summaries[['institutionID', 'collectionCode',
                                         'datasetName', 'acquired', 'removed',
                                         'retained']]

    # Save the summaries in the output database
    summary_df.to_sql(name='attribute_value_counts', con = conn, if_exists='replace')
    source_summaries.to_sql(name='sources', con = conn, if_exists='replace')
    print("Saved summary of filtering results: " + str(datetime.now() - timestamp))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SAVE
    # Reformat data to strings and insert into db.
    df_filterZ.replace("nan", pd.NA).applymap(str).to_sql(name='occurrence_records', con = conn,
                                  if_exists='replace')
    conn.close()
    return None

def nominal_precisions(longitude, latitude, produce):
    '''
    Calculates the nominal precisions based on an WGS84 coordinates.
    Method is based on information from wikipedia page on latitude and posts
    at https://gis.stackexchange.com/questions/8650/measuring-accuracy-of-latitude-and-longitude
    https://wiki.openstreetmap.org/wiki/Precision_of_coordinates

    PARAMETERS
    ----------
    latitude : decimal degrees (EPSG:4326) latitude as string.
    longitude : decimal degrees (EPSG:4326) longitude as string.
    produce : 'longitude', 'latitude', or 'both'

    RETURNS
    -------
    x : uncertainty in longitude (meters) as float.
    y : uncertianty in latitude (meters) as float.

    EXAMPLE
    -------
    x, y = nominal_precisions("-93.455", "26.3455", produce="both")
    '''
    import numpy as np

    lat = latitude.split(".")
    long = longitude.split(".")

    # Longitude
    digitsX = {1: 10, 2: 100, 3: 1000, 4: 10000, 5: 100000}
    x =(111321 * np.cos(float(latitude) * np.pi/180))/digitsX[len(long[1])] # decimal gets moved based on digits.

    # Latitude
    digitsY = {1: 11112.0, 2: 1111.2, 3: 111.1, 4: 11.1, 5: 1.1} # Lookup for latitude precision
    y = digitsY[len(lat[1])]

    if produce == "both":
        return x, y
    if produce == "longitude":
        return x
    if produce == "latitude":
        return y

def drop_duplicates_latlongdate(df):
    '''
    Function to find and remove duplicate occurrence records within the
    wildlife wrangler workflow.  When duplicates exist, the record with the
    higher decimal precision is kept, and if precision values are equal, then the
    record with the smallest radius_m is retained. Accounts for existence
    of records with a mix of decimal precision in latitude and longitude
    values. The process is a little complex.  The first data frame is cleaned
    up by dropping duplicates based on which record has smaller buffer radius.
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
    #df = df.astype(output_schema)

    """
    ############ RECTIFY UNEQUAL LAT-LONG PRECISION
    First, trim decimal length in cases where decimal length differs between
    latitude and longitude values, result is equal latitude and longitude
    length.  Record the trimmed decimal precision in a temp column for use later
    as a record to "verbatim" precision.
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
    keeping the first (lowest radius_m)
    Sort so that the lowest radius_m is first
    """
    df = (df
          .sort_values(by=['decimalLatitude', 'decimalLongitude', 'eventDate',
                           'radius_m'],
                       ascending=True, kind='mergesort', na_position='last')
          .drop_duplicates(subset=['decimalLatitude', 'decimalLongitude',
                                   'eventDate'],
                           keep='first'))

    """
    #########  FIND IMPRECISE DUPLICATES
    Get a list of "verbatim" precisions that are present in the data to loop through.
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
    df2.drop(columns=['dup_OGprec'], axis=1, inplace=True)

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

def spatial_output(database, make_file, mode, output_file=None, epsg=4326):
    '''
    Creates a shapefile of species occurrence records from a wildlife wrangler
        output SQLite database.

    PARAMETERS
    ----------
    database : the sqlite database to use.
    make_file : whether to save a shapefile.  False just returns a geodatframe.
    mode : three options: 1) "points" merely creates points from record coordinates.
        2) "footprints" uses the point-radius and shape methods to map polygons.
        3) "random" creates the footprints and then picks a random point within
        each footprint polygon.
    output_file : Path (and name) of the file to be created.
    epsg : integer epsg number for the desired geographic projection for results.

    OUTPUT
    ------
    geopandas data frame and a shapefile saved to disk
    '''
    from datetime import datetime
    import sqlite3
    import pandas as pd
    import numpy as np
    import geopandas as gpd
    import random
    from shapely.geometry import Point
    from shapely import wkt

    timestamp = datetime.now()

    # Get the record coordinates as a data frame
    records = (pd.read_sql("""SELECT * FROM occurrence_records;""",
                           con=sqlite3.connect(database))
                           .astype({'decimalLongitude': 'float',
                                    'decimalLatitude': 'float',
                                    'radius_m': 'float'}))

    # Make a point geometry from coordinates
    gdf = gpd.GeoDataFrame(records,
                           geometry=gpd.points_from_xy(records['decimalLongitude'],
                                                       records['decimalLatitude']))

    # Set the coordinate reference system and reproject to Albers
    gdf.crs={"init" : "epsg:4326"}
    #gdf = gdf.to_crs(epsg=5070)

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  POINTS
    # If user just wants shapefile of given coordinates (point method)
    if mode == "points":
        if epsg != 4326:
            gdf = gdf.to_crs(epsg=epsg)

        # Save points
        if make_file == True:
            gdf.to_file(output_file)
            print("Generated points")

        out = gdf

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  FOOTPRINTS
    else:
        # Point-radius method -- needs to be EPSG 5070
        gdf = gdf.to_crs(epsg=5070)
        gdf["footprint"] = gdf.apply(lambda x: x.geometry.buffer(x.radius_m), axis=1)
        gdf.set_geometry(col="footprint", inplace=True, drop=True)
        gdf = gdf.to_crs(epsg=4326)
        gdf["footprint"] = gdf["geometry"]

        # Use the footprintWKT when provided, be sure to reproject (shape method)
        if list(gdf['footprintWKT'].unique()) != ["<NA>"]:
            print("Some footprintWKT were provided.")
            # Separate data frame for records with footprintWKT, load WKT
            ftp = gdf[gdf["footprintWKT"] != "<NA>"].copy()

            # Points aren't appropriate footprintWKT for this framework - remove.
            ftp = ftp[~ftp["footprintWKT"].str.contains("POINT")].copy()

            if ftp.empty == False:
                ftp["footprint"] = ftp["footprintWKT"].apply(lambda x: wkt.loads(str(x)))

                # Handle geometry
                ftp.set_geometry(col="footprint", inplace=True, drop=False)
                ftp.crs={"init" : "epsg:4326"}
                ftp = ftp.copy()[["record_id", "footprint"]]

                # Combine the two data frame's geoseries, with preference to footprint
                gdf.update(ftp)

        # Reset the geometry column
        gdf.set_geometry(col='footprint', inplace=True, drop=True)

        if mode != "random":
            print("Generated footprints")

        # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  RANDOM POINTS
        # If the user wants a random point selected from within each footprint
        if mode == "random":
            # Function to generate random points
            def generate_random(number, polygon):
                points = []
                minx, miny, maxx, maxy = polygon.bounds
                while len(points) < number:
                    pnt = Point(random.uniform(minx, maxx), random.uniform(miny, maxy))
                    if polygon.contains(pnt):
                        points.append(pnt)
                return points

            # Generate random points
            gdf["random_points"] = gdf["geometry"].apply(lambda x: generate_random(1, x)[0])
            gdf.drop(["geometry"], axis=1, inplace=True)
            gdf.set_geometry(col='random_points', inplace=True, drop=True)
            print("Generated random points")

        # Reproject
        if epsg != 4326:
            gdf = gdf.to_crs(epsg=epsg)

        # Save points
        #gdf.drop(["geometry"], axis=1, inplace=True)
        if make_file == True:
            gdf.to_file(output_file)

        out = gdf

    print("Finished: " + str(datetime.now() - timestamp))
    return out

def CONUS_bbox():
    """
    Returns the bounding box of the conterminous U.S. as a tuple.
    """
    bbox = (-171.791110603, 18.91619, -66.96466, 71.3577635769)

def verify_results(database):
    '''
    Compares the occurrence record attributes to the filters that were supposed to be applied.


    PARAMETERS
    ---------
    database: path to a wrangler output database.  Like ""Z:/Occurrence_Records/test1.sqlite""

    RESULTS
    -------
    prints messages if tests are failed. No output indicates all tests were passed.
    '''
    import sqlite3
    import pandas as pd
    import geopandas as gpd
    import shapely

    # Connect to a database
    conn = sqlite3.connect(database)
    cursor = conn.cursor()

    # Get the taxon concept ----------------------------------------------------
    taxon_concept = (pd.read_sql(sql="SELECT * FROM taxon_concept;", con=conn)
                     .rename({"index": "key", "0": "value"}, axis=1)
                     .set_index("key"))

    # Get the filter set that was applied --------------------------------------
    filter_set = (pd.read_sql(sql="SELECT * FROM filter_set", con=conn)
                  .rename({"index": "parameter", "0": "value"}, axis=1)
                  .set_index("parameter"))

    # Presence only ------------------------------------------------------------
    presabs = [x[0] for x in conn.execute("SELECT DISTINCT occurrenceStatus FROM occurrence_records;").fetchall()]
    if presabs != ["PRESENT"] and presabs != ['nan', 'PRESENT']:
        print("!!! Failed occurrenceStatus test: " + str(presabs))

    # DWCA archive -------------------------------------------------------------
    gbif_yn = [x[0] for x in conn.execute("SELECT DISTINCT source FROM occurrence_records;").fetchall()]
    if 'GBIF' in gbif_yn:
        if filter_set.loc["get_dwca", "value"] == "True":
            try:
                conn.execute("SELECT download_key FROM GBIF_download_info;").fetchall()
            except:
                print("!!! Failed DWCA download key test")

    # Latitude range -----------------------------------------------------------
    lat_range = filter_set.loc["lat_range", "value"]
    if lat_range != "None":
        lat_range = [float(x) for x in lat_range.split(",")]
        max_lat = float(conn.execute("SELECT MAX(decimalLatitude) FROM occurrence_records;").fetchall()[0][0])
        min_lat = float(conn.execute("SELECT MIN(decimalLatitude) FROM occurrence_records;").fetchall()[0][0])
        if (min_lat < lat_range[0] or max_lat > lat_range[1]):
            print("!!! Failed test of decimalLatitude values")

    # Longitude range ----------------------------------------------------------
    lon_range = filter_set.loc["lon_range", "value"]
    if lon_range != "None":
        lon_range = [float(x) for x in lon_range.split(",")]
        max_lon = float(conn.execute("SELECT MAX(decimalLongitude) FROM occurrence_records;").fetchall()[0][0])
        min_lon = float(conn.execute("SELECT MIN(decimalLongitude) FROM occurrence_records;").fetchall()[0][0])
        if (min_lon < lon_range[0] or max_lon > lon_range[1]):
            print("!!! Failed test of decimalLongitude values")

    # Years range --------------------------------------------------------------
    yrs_range = filter_set.loc["years_range", "value"]
    if yrs_range != "None":
        yrs_range = [float(x) for x in yrs_range.split(",")]
        max_yr = float(conn.execute("SELECT MAX(strftime('%Y', eventDate)) FROM occurrence_records;").fetchall()[0][0])
        min_yr = float(conn.execute("SELECT MIN(strftime('%Y', eventDate)) FROM occurrence_records;").fetchall()[0][0])
        if (min_yr < yrs_range[0] or max_yr > yrs_range[1]):
            print("!!! Failed test of year (eventDate) values")

    # Months range -------------------------------------------------------------
    months_range = filter_set.loc["months_range", "value"]
    if months_range != "None":
        months_range = [float(x) for x in months_range.split(",")]
        max_month = float(conn.execute("SELECT MAX(strftime('%m', eventDate)) FROM occurrence_records;").fetchall()[0][0])
        min_month = float(conn.execute("SELECT MIN(strftime('%m', eventDate)) FROM occurrence_records;").fetchall()[0][0])
        # Months range could be like 1,12
        if months_range[0] < months_range[1]:
            if (min_month < months_range[0] or max_month > months_range[1]):
                print("!!! Failed test of month (eventDate) values")
        # Months range could be like 11,3
        if months_range[0] > months_range[1]:
            no_months = list(range(months_range[1] + 1, months_range[0]) -1)
            months = conn.execute("SELECT DISTINCT strftime('%m', eventDate) FROM occurrence_records;").fetchall()
            months = [int(x[0]) for x in months]
            if len(set(months) & set(no_months)) != 0:
                print("!!! Failed test of month (eventDate) values")

    # eBird ID -----------------------------------------------------------------
    ebd_yn = [x[0] for x in conn.execute("SELECT DISTINCT source FROM occurrence_records;").fetchall()]
    ebirdid = conn.execute("SELECT DISTINCT ebird_id FROM occurrence_records;").fetchall()[0][0]
    out_ebirdid = taxon_concept.loc["EBIRD_ID", "value"]
    if 'eBird' in ebd_yn:
        if ebirdid != out_ebirdid:
            print("!!! Failed test of ebird_id")

    # GBIF ID ------------------------------------------------------------------
    gbifid = conn.execute("SELECT DISTINCT gbif_id FROM occurrence_records;").fetchall()[0][0]
    out_gbifid = taxon_concept.loc["GBIF_ID", "value"]
    if gbifid != out_gbifid:
        print("!!! Failed test of gbif_id")

    # Maximum coordinate uncertainty -------------------------------------------
    mcu_out = float(conn.execute("SELECT MAX(coordinateUncertaintyInMeters) FROM occurrence_records;").fetchall()[0][0])
    dd = float(taxon_concept.loc["detection_distance_m", "value"])
    mcu = float(filter_set.loc["max_coordinate_uncertainty", "value"])
    if mcu_out > mcu:
        print("!!! Failed test of maximum coordinate uncertainty. Max = " + str(mcu_out))

    # Various ------------------------------------------------------------------
    multilists = {"issues_omit": "issues", "sampling_protocols_omit": "samplingProtocol", "bases_omit": "basisOfRecord",
                  "datasets_omit": "datasetName", "collection_codes_omit": "collectionCode", "institutions_omit": "institutionID"}

    def test_multilist(attribute):
        '''Test attribute values that involove multilists'''
        values = set([])
        data = [x[0].split("; ") for x in conn.execute("SELECT DISTINCT {0} FROM occurrence_records;".format(multilists[attribute])).fetchall()]
        for x in data:
            values = set(x) | values
        invalid = set(filter_set.loc[attribute, "value"].replace("'", "").replace("[", "").replace("]", "").split(", "))
        if len(invalid & values) >=1:
            print("!!! Failed test of {0}.".format(multilists[attribute]))

    for x in multilists.keys():
        test_multilist(x)

    # Duplicates ---------------------------------------------------------------
    if filter_set.loc["duplicates_OK", "value"] == "False":
        records = pd.read_sql("SELECT decimalLatitude, decimalLongitude, eventDate FROM occurrence_records;", con=conn)
        if len(records[records.duplicated() == True]) >= 1:
            print("!!! Failed test for duplicates")

    # Test spatial parameters --------------------------------------------------
    records2 = (pd.read_sql("""SELECT * FROM occurrence_records;""", con=conn)
                .astype({'decimalLongitude': 'float',
                         'decimalLatitude': 'float',
                         'radius_m': 'float'}))

    # Make a point geometry from coordinates
    gdf = gpd.GeoDataFrame(records2,
                           geometry=gpd.points_from_xy(records2['decimalLongitude'],
                                                       records2['decimalLatitude']))

    # Set the coordinate reference system
    gdf.crs={'init' :'epsg:4326'}

    # Test species extent of occurrence polygon
    if filter_set.loc["use_taxon_geometry", "value"] == "True":
        if taxon_concept.loc["TAXON_EOO", "value"] != "None":
            poly = shapely.wkt.loads(taxon_concept.loc["TAXON_EOO", "value"])
            gdf_eoo = gdf[~gdf["geometry"].within(poly)]
            if len(gdf_eoo) >= 1:
                print("!!! Failed species extent of occurrence test.")

    # Test query polygon
    if filter_set.loc["query_polygon", "value"] != "None":
        poly = shapely.wkt.loads(filter_set.loc["query_polygon", "value"])
        gdf_query = gdf[~gdf["geometry"].within(poly)]
        if len(gdf_query) >= 1:
            print("!!! Failed query polygon test.")

    conn.close()
    return
