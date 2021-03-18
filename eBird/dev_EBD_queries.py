working_directory = "T:/temp/"
EBD_file = "T:/data/eBird/ebd_US_relDec-2020.txt"
queried_ebd = working_directory + "ebd_queried.txt"
processed_ebd = working_directory + "ebd_filtered.csv"
species = "Yellow-billed Cuckoo"
country = "US"
months_range = "3,6"
years_range = "2018,2020"
lon_range = "-89,-75"
lat_range = "27,41"
max_coordinate_uncertainty = 10000
collection_codes_omit = "EBIRD_ATL_VA"
sampling_protocols_omit = "Historical", "Stationary"
taxon_polygon = "POLYGON ((-84.09680233298448 36.69265225442667, -84.07962135716329 34.5561660300382, -84.07962135716329 34.5561660300382, -80.25685423694925 34.65515526072436, -81.15026497965096 36.71331438415306, -84.09680233298448 36.69265225442667))"
query_polygon = "POLYGON ((-82.74809573102132 36.96082629937069, -85.0932989306133 35.63154639485496, -81.0987220521874 33.56697226279766, -79.4235769096217 36.34054727735634, -79.4235769096217 36.34054727735634, -82.74809573102132 36.96082629937069))"

def get_EBD_records(EBD_file, queried_ebd, processed_ebd, species, country,
                    months_range, years_range, lon_range, lat_range,
                    max_coordinate_uncertainty, collection_codes_omit,
                    sampling_protocols_omit, taxon_polygon, query_polygon):
    """
    Gets eBird records from a copy of the Ebird Basic Dataset that you acquired.
    Primarily runs R code that uses the auk package to query the data set in
    an efficient manner.  Some filters can be applied during the query, but
    others have to be applied to the query results.  Date and bounding box
    filters require quite a bit of preparation and conditions.

    Parameters
    ----------
    EBD_file : path to your downloaded copy of the Ebird Basic Dataset; string
    queried_ebd : path to use for an intermediate table of query results; string
    processed_ebd : path to use for table of filtered query results; string
    species : eBird common name for the species; string
    country : two letter country code; like "US"; string
    months_range : start and end months of interest; like "3,6"; string
    years_range : start and end years of interest; like "2000,2005"; string
    lon_range : range of desired longitudes for bounding box filtering; like "-89,-75"; string
    lat_range : range of desired latitudes for bonding box filtering; like "27,41"; string
    max_coordinate_uncertainty : maximum traveling count distance in m; integer
    collection_codes_omit :
    sampling_protocols_omit :
    taxon_polygon :
    query_polygon :

    Returns
    -------
    Data frame of filtered ebird records
    """
    import rpy2.robjects as robjects
    import rpy2.robjects.packages as rpackages
    from rpy2.robjects.vectors import StrVector
    import pandas as pd

    # import R's utility package, select a mirror for R packages
    utils = rpackages.importr('utils')
    utils.chooseCRANmirror(ind=1) # select the first mirror in the list

    # R packages to load
    packnames = ('sf', 'auk', 'lubridate', 'tidyverse')
    names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
    if len(names_to_install) > 0:
        utils.install_packages(StrVector(names_to_install))

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
    collection_codes_omit <- "{10}"
    sampling_protocols_omit <- "{11}"
    taxon_polygon <- "{12}"
    query_polygon <- "{13}"

    # Auk uses filters that are compiled and incorporated into a query. This poses
    # a challenge for dynamic filtering where filters may or may not be used.  We
    # have to set defaults.
    library(auk)
    library(tidyverse)
    library(sf)
    library(lubridate)
    starttime = Sys.time() # Runtime has been 30 min

    # prep dates -------------------------------------------------------------------
    # auk doesn't allow filtering on months AND year with read_ebd; they have to be
    #   done separately, one with auk filters and the other after with dplyr.  I
    #   chose to do the year filtering with auk to minimize size of returned tibble.
    #   This all requires formatting dates as text correctly
    # format start month
    if (months_range != "") {{
      if (length(strsplit(months_range, ",")[[1]][1]) == 1) {{
        start_month <- paste(c("0", strsplit(months_range, ",")[[1]][1]), collapse="")
        }} else {{
        start_month <- strsplit(months_range, ",")[[1]][1]}}
      start_month <- str_trim(start_month)

      # format end month
      if (length(strsplit(months_range, ",")[[1]][2]) == 1) {{
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

    # prep bounding box ------------------------------------------------------------
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

    # prep country -----------------------------------------------------------------
    if (country == "") {{country <- "US"}}

    # prep distance ----------------------------------------------------------------
    #   a gps precision for eBird checklists must be assumed, since not given, for
    #   estimation of coordinateUncertaintyInMeters
    EBD_gps_precision <- 10

    # account for gps precision in distance filter.  error could exist on either
    #   end of a straight line path, so double precision when subtracting.
    max_distance <- as.integer(ceiling((max_coordinate_uncertainty-2*EBD_gps_precision)/1000))  ### NOT WORKING

    # query -----------------------------------------------------------------
    ebd_data_0 <- EBD_file %>%
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

    # prep data frame for python ----------------------------------------
    # add column for eBird species code
    ebird_code <- select(filter(ebird_taxonomy, common_name==species),
                         species_code)[[1]]
    ebd_data <- ebd_data_0 %>%
                 mutate(eBird_sp_code = ebird_code,
                        retrieval_date = Sys.time()) %>%
                 select(eBird_sp_code, global_unique_identifier, checklist_id,
                        project_code, last_edited_date, common_name,
                        observation_count, locality, latitude, longitude,
                        observation_date, observer_id, effort_distance_km,
                        protocol_type, effort_area_ha, trip_comments,
                        species_comments) %>%
                 mutate(effort_distance_m = as.numeric(effort_distance_km)*1000) %>%
                 mutate(coordinateUncertaintyInMeters = effort_distance_m + (2*EBD_gps_precision)) %>%
                 #filter(coordinateUncertaintyInMeters <= max_coordinate_uncertainty) %>%  This should be true, test so.
                 filter(month(observation_date) %in% ok_months) %>%
                 filter(!project_code %in% c(collection_codes_omit)) %>%
                 filter(!protocol_type %in% c(sampling_protocols_omit)) %>%
                 write_csv(processed_ebd)

    endtime = Sys.time()
    print(endtime - starttime)
    '''.format(EBD_file, queried_ebd, processed_ebd, species, country, months_range,
    years_range, lon_range, lat_range, max_coordinate_uncertainty, collection_codes_omit,
    sampling_protocols_omit, taxon_polygon, query_polygon)

    # Run code
    robjects.r(code)

    # Read output data frame
    ebd_data = pd.read_csv(processed_ebd)
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
        ebd_data = robjects.conversion.ri2py(rdf)
    print(type(ebd_data))
    '''
    return ebd_data

ebird_data = get_EBD_records(EBD_file, queried_ebd, processed_ebd, species, country,
                             months_range, years_range, lon_range, lat_range,
                             max_coordinate_uncertainty, collection_codes_omit,
                             sampling_protocols_omit, taxon_polygon, query_polygon)
