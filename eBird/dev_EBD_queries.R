library(auk)
library(tidyverse)
starttime = Sys.time() # Runtime has been 30 min

# retrieve filter sets from parameters db -------------------------------
# ebd needs its own tables for requests and filtering?

# path to the ebird data file -------------------------------------------
#input_file <- system.file("extdata/ebd-sample.txt", package = "auk")
input_file <- "T:/data/eBird/ebd_US_relDec-2020.txt"

# parameters ------------------------------------------------------------
queried_ebd <- "T:/temp/ebd_queried.txt"
processed_ebd <- "T:/temp/ebd_filtered.csv"
species <- c("Swainson's Warbler")
country <- "US"
date <- c("2015-01-01", "2021-12-31")
bbox <- ""#c(lng_min, lat_min, lng_max, lat_max)
distance <- c(0, 10)

# query -----------------------------------------------------------------
ebd_data <- input_file %>%
  # 1. reference file
  auk_ebd() %>%
  # 2. define filters
  auk_species(species=species) %>%
  auk_date(date=date) %>%
  auk_country(country=country) %>%
  auk_bbox(bbox=bbox) %>%
  auk_distance(distance=distance) %>%

  # 3. run filtering
  auk_filter(file = queried_ebd, overwrite = TRUE) %>%
  # 4. read text file into r data frame
  read_ebd()

# display -------------------------------------------------------------
View(ebird_data)

# prep data frame for wrangler ----------------------------------------
# add column for eBird species code
# get eBird species species code
ebird_code <- select(filter(ebird_taxonomy, common_name==species),
                     species_code)[[1]]
ebd_data2 <- ebd_data %>%
             mutate(eBird_sp_code = ebird_code,
                    retrieval_date = Sys.time()) %>%
             select(eBird_sp_code, global_unique_identifier, checklist_id,
                    last_edited_date, common_name, observation_count, locality,
                    latitude, longitude, observation_date, observer_id,
                    effort_distance_km, effort_area_ha, trip_comments,
                    species_comments) %>%
             mutate(effort_distance_m = effort_distance_km/1000) %>%
             write_csv(processed_ebd)
View(ebd_data2)

endtime = Sys.time()
print(endtime - starttime)
