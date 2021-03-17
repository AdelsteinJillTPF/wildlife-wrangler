# Auk uses filters that are compiled and incorporated into a query. This poses
# a challenge for dynamic filtering where filters may or may not be used.  We
# have to set defaults.
library(auk)
library(tidyverse)
library(sf)
library(lubridate)
starttime = Sys.time() # Runtime has been 30 min

# path to the ebird data file --------------------------------------------------
#input_file <- system.file("extdata/ebd-sample.txt", package = "auk")
input_file <- "T:/data/eBird/ebd_US_relDec-2020.txt"
queried_ebd <- "T:/temp/ebd_queried.txt"
processed_ebd <- "T:/temp/ebd_filtered.csv"

# parameters -------------------------------------------------------------------
species <- "Yellow Warbler"
country <- "US"
months_range <- "3,6"
years_range <- "2018,2020"
lon_range <- "-89,-75"
lat_range <- "27,41"
max_coordinate_uncertainty <- 10000
collection_codes_omit <- c("EBIRD")
sampling_protocols_omit <- c("Historical", "Stationary")
taxon_polygon <- "POLYGON ((-84.09680233298448 36.69265225442667, -84.07962135716329 34.5561660300382, -84.07962135716329 34.5561660300382, -80.25685423694925 34.65515526072436, -81.15026497965096 36.71331438415306, -84.09680233298448 36.69265225442667))"
query_polygon = "POLYGON ((-82.74809573102132 36.96082629937069, -85.0932989306133 35.63154639485496, -81.0987220521874 33.56697226279766, -79.4235769096217 36.34054727735634, -79.4235769096217 36.34054727735634, -82.74809573102132 36.96082629937069))"

# prep dates -------------------------------------------------------------------
# auk doesn't allow filtering on months AND year with read_ebd; they have to be
#   done separately, one with auk filters and the other after with dplyr.  I
#   chose to do the year filtering with auk to minimize size of returned tibble.
#   This all requires formatting dates as text correctly
# format start month
if (months_range != "") {
  if (length(strsplit(months_range, ",")[[1]][1]) == 1) {
    start_month <- paste(c("0", strsplit(months_range, ",")[[1]][1]), collapse="")
    } else {
    start_month <- strsplit(months_range, ",")[[1]][1]}
  start_month <- str_trim(start_month)

  # format end month
  if (length(strsplit(months_range, ",")[[1]][2]) == 1) {
    end_month <- paste(c("0", strsplit(months_range, ",")[[1]][2]), collapse="")
    } else {
    end_month <- strsplit(months_range, ",")[[1]][2]}
  end_month <- str_trim(end_month)

  # create vector of ok months for filtering with dplyr
  ok_months <- seq(as.numeric(start_month), as.numeric(end_month))
}

# pull out start and end years
if (years_range != "") {
start_yr <- str_trim(strsplit(years_range, ",")[[1]][1])
end_yr <- str_trim(strsplit(years_range, ",")[[1]][2])
}

# define data filter according to different scenarios
if (months_range == "" && years_range == "") {
  # get all dates
  date_filter <- c("*-01-31", "*-12-31")
} else if (months_range != "" && years_range != "") {
  # get certain months and years.  we have to find the last possible day.
  end_day <- lubridate::days_in_month(paste(c(end_yr, "-", end_month, "-01"),
                                            collapse=""))
  date_filter <- c(paste(c(start_yr, "-", start_month, "-01"), collapse=""),
                   paste(c(end_yr, "-", end_month, "-", end_day), collapse=""))
} else if (months_range == "" && years_range != "") {
  # get all months from certain years
  date_filter <- c(paste(c(start_yr, "-01-01"), collapse=""),
                   paste(c(end_yr, "-12-31"), collapse=""))
} else if (months_range != "" && years_range == "") {
  # get certain months from all years.  we have to find the last possible day.
  yr <- year(today())
  end_day <- lubridate::days_in_month(paste(c(yr, "-", end_month, "-01"),
                                            collapse=""))
  date_filter <- c(paste(c("*-", start_month, "-01"), collapse=""),
                   paste(c("*-", end_month, "-", end_day), collapse=""))
}

# prep bounding box ------------------------------------------------------------
bbox <- NULL
if (query_polygon == "" && taxon_polygon == "") {
    bbox <- NULL
} else if (query_polygon != "" && taxon_polygon == "") {
    bbox <- st_bbox(st_as_sfc(query_polygon))
} else if (query_polygon == "" && taxon_polygon != "") {
    bbox <- st_bbox(st_as_sfc(taxon_polygon))
} else if (query_polygon != "" && taxon_polygon != "") {
    # Get/use the intersection of the two polygons
    filter_polygon <- st_as_sfc(query_polygon)
    sp_polygon <- st_as_sfc(taxon_polygon)
    bbox <- st_bbox(st_intersection(filter_polygon, sp_polygon))
}

# prep bounding box vector for filter if lat and lon ranges were provided,
#   and if other polygons were not
if (lat_range == "" || lon_range == "") {
  null_box <- TRUE
} else {
  null_box <- FALSE
}

if (is.null(bbox) == TRUE && null_box == FALSE) {
  lat_min <- as.numeric(strsplit(lat_range, ",")[[1]][1])
  lat_max <- as.numeric(strsplit(lat_range, ",")[[1]][2])
  lng_min <- as.numeric(strsplit(lon_range, ",")[[1]][1])
  lng_max <- as.numeric(strsplit(lon_range, ",")[[1]][2])
  bbox <- c(lng_min, lat_min, lng_max, lat_max)
}

# prep country -----------------------------------------------------------------
if (country == "") {country <- "US"}

# prep distance ----------------------------------------------------------------
#   a gps precision for eBird checklists must be assumed, since not given, for
#   estimation of coordinateUncertaintyInMeters
EBD_gps_precision <- 10

# account for gps precision in distance filter.  error could exist on either
#   end of a straight line path, so double precision when subtracting.
max_distance <- (max_coordinate_uncertainty-(2*EBD_gps_precision))/1000

# query -----------------------------------------------------------------
ebd_data <- input_file %>%
  # 1. reference file
  auk_ebd() %>%

  # 2. define filters
  auk_species(species=c(species)) %>%
  auk_date(date=date_filter) %>%
  auk_country(country=country) %>%
  auk_bbox(bbox=bbox) %>%
  auk_distance(distance=c(0, max_distance) %>%

  # 3. run filtering
  auk_filter(file = queried_ebd, overwrite = TRUE) %>%

  # 4. read text file into r data frame
  read_ebd()

# display -------------------------------------------------------------
View(ebd_data)

# prep data frame for wrangler ----------------------------------------
# add column for eBird species code
# get eBird species code
ebird_code <- select(filter(ebird_taxonomy, common_name==species),
                     species_code)[[1]]
ebd_data2 <- ebd_data %>%
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
             filter(!project_code %in% collection_codes_omit) %>%
             filter(!protocol_type %in% sampling_protocols_omit) %>%
             write_csv(processed_ebd)
View(ebd_data2)

endtime = Sys.time()
print(endtime - starttime)


# tests ------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
max(ebd_data2$coordinateUncertaintyInMeters) <= max_coordinate_uncertainty
max(ebd_data2$effort_distance_m)

library(maps)
demo(nc, ask = FALSE, echo = FALSE)
usa = st_as_sf(map('usa', plot = FALSE, fill = TRUE))
sf <- st_as_sf(ebd_data2, coords=c("longitude", "latitude"), crs=4326)
bbox2 <- st_as_sfc(bbox) %>% st_set_crs(4326)
ggplot() +
  geom_sf(data=usa) +
  geom_sf(data=st_as_sfc(query_polygon, crs=4326), aes()) +
  geom_sf(data=st_as_sfc(taxon_polygon, crs=4326), aes()) +
  geom_sf(data=sf, aes()) +
  geom_sf(data=bbox2, aes())
