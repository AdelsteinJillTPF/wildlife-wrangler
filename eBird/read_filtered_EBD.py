import pandas as pd

# Read in csv
df = pd.read_csv("T:/temp/ebd_filtered.csv")
# Add some columns (and populated cells)
df["basisOfRecord"] = "HumanObservation"
df["datasetName"] = "EBird_Basic_Dataset"
df["institutionCode"] = "CLO"
df["samplingProtocol"] = "birdwatching"
df["remarks"] = df["trip_comments"] + df["species_comments"]
df["source"] = "eBird"

# Rename columns to fit wrangler schema
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
               "obseration_count": "individualCount",
               ""
               }
