configDir = "T:/Data/"  # Path to folder where you saved your wildlifeconfig file.
query_name = "geomgbifdebug"

import pandas as pd
import numpy as np
from datetime import datetime
import sys
sys.path.append(configDir)
import wranglerconfig as config
codeDir = config.codeDir
sys.path.append(codeDir)
import wrangler_functions as functions
from pygbif import occurrences
from dwca.read import DwCAReader
working_directory = config.workDir
output_db = working_directory + query_name + '.sqlite'
username = config.gbif_username
password = config.gbif_password
email = config.gbif_email
EBD_file = config.EBD_file


download_filters = ['taxonKey = 2496287', 'hasCoordinate = True', 'country = US',
                    'year >= 2011', 'year <= 2020', 'month >= 6', 'month <= 8',
                    'decimalLatitude >= 30', 'decimalLatitude <= 35',
                    'decimalLongitude >= -89', 'decimalLongitude <= -75']

# Get the value of the download key
d = occurrences.download(download_filters,
                         pred_type='and',
                         user = username,
                         pwd = password,
                         email = email)
dkey = d[0]
print(dkey)

print("Downloading Darwin Core Archive zip file for this species .....")
timestamp = datetime.now()
gotit = False
while gotit == False:
    try:
        zipdownload = occurrences.download_get(key=dkey, path=working_directory)
    except Exception as e:
        print(e)

    try:
        # Read the relevant files from within the darwin core archive
        with DwCAReader(working_directory + dkey + '.zip') as dwca:
            dfRaw = dwca.pd_read('occurrence.txt', low_memory=False)
            citations = dwca.open_included_file('citations.txt').read()
            rights = dwca.open_included_file('rights.txt').read()
            doi = dwca.metadata.attrib["packageId"]
            gotit = True
        print("Download complete: " + str(datetime.now() - timestamp))
    except Exception as e:
        wait = datetime.now() - timestamp
        if wait.seconds > 60*45:
            gotit = True
            print("TIMED OUT")
        else:
            pass

records1 = dfRaw.rename({"issue": "issues", "id": "record_id"}, axis=1)
# Drop columns
records1 = records1.filter(items=functions.output_schema.keys(), axis=1)
# Populate columns
records1["retrieval_date"] = str(datetime.now())
records1["GBIF_download_doi"] = doi
records1["source"] = "GBIF"

records2 = (pd.DataFrame(columns=functions.output_schema.keys())
            .combine_first(records1)
            #.replace({"coordinateUncertaintyInMeters": {"UNKNOWN": np.nan},
            #               "radius_m": {"UNKNOWN": np.nan},
            #               "individualCount": {"UNKNOWN": 1},
            #               "weight": {"UNKNOWN": 10},
            #               "detection_distance_m": {"UNKNOWN": 0}})
            .fillna({"coordinateUncertaintyInMeters": 0,
                     "radius_m": 0,
                     "individualCount": 1,
                     "weight": 10,
                     "detection_distance_m": 0,
                     "effort_distance_m": 0,
                     "coordinate_precision": 1,
                     "gps_accuracy_m": 30})
            .astype(functions.output_schema)
            )


records3 = records2[["coordinateUncertaintyInMeters", "radius_m", "individualCount", "weight", "detection_distance_m"]]


#####################################################################
####################     FROM GITHUB ISSUE     ######################
#####################################################################
download_filters2 = ["taxonKey = 2496287",
'geometry = "POLYGON((-82.77 36.97, -85.07 35.67, -81.07 33.57, -79.74 36.37, -79.47 36.37, -82.77 36.97))"']

download = occurrences.download(download_filters2,
                                pred_type='and',
                                user = username,
                                pwd = password,
                                email = email)
