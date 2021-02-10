import pandas as pd

# List of columns to use
cols = ['checklist_id', 'global_unique_identifier', 'last_edited_date',
       'common_name', 'scientific_name', 'observation_count', 'locality',
       'latitude', 'longitude', 'observation_date', 'observer_id',
       'sampling_event_identifier', 'project_code', 'effort_distance_km',
       'group_identifier', 'approved', 'reviewed', 'trip_comments',
       'species_comments', 'eBird_code']

# Read in csv, rename columns,
df = (pd.read_csv("T:/temp/ebd_filtered.csv", usecols=cols))

# Add some columns (and populated cells)
df["basisOfRecord"] = "HumanObservation"
df["datasetName"] = "EBird_Basic_Dataset"
df["institutionCode"] = "CLO"
df["samplingProtocol"] = "birdwatching"
df["remarks"] = df["trip_comments"] + df["species_comments"]
df["source"] = "eBird"

# Map for renaming columns to fit wrangler schema
rename_dict = {"ebird_sp_code": "taxonConceptID",
               "last_edited_date": "modified",
               "observer_id": "recordedBy",
               "locality": "Locality",
               "project_code": "collectionCode",
               "effort_distance_m": "coordinateUncertaintyInMeters",
               "latitude": "decimalLatitude",
               "longitude": "decimalLongitude",
               "observation_date": "eventDate",
               "global_unique_identifier": "record_id",
               "obseration_count": "individualCount"}
df.rename(rename_dict, axis=1, inplace=True)
