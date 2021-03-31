configDir = "T:/Data/"  # Path to folder where you saved your wildlifeconfig file.
working_directory = "T:/Occurrence_Records/"
output_db = "T:/Occurrence_Records/withRTest.sqlite"
query_name = 'withRTest2'
taxon_info = {'EBIRD_ID': 'Yellow-billed Cuckoo',
    'GBIF_ID': 2496287,
    'ID': 'TestCuckoo',
    'TAXON_EOO': 'POLYGON ((-84.09680233298448 36.69265225442667, '
          '-84.07962135716329 34.5561660300382, -84.07962135716329 '
          '34.5561660300382, -80.25685423694925 34.65515526072436, '
          '-81.15026497965096 36.71331438415306, -84.09680233298448 '
          '36.69265225442667))'}
filter_set = {'bases_omit': '',
 'collection_codes_omit': '',
 'country': 'US',
 'datasets_omit': '',
 'default_coordUncertainty': 1000,
 'duplicates_OK': False,
 'geoissue': False,
 'get_dwca': True,
 'has_coordinate_uncertainty': True,
 'institutions_omit': '',
 'issues': '',
 'lat_range': '27,41',
 'lon_range': '-89,-75',
 'max_coordinate_uncertainty': 10000,
 'months_range': '1,12',
 'name': 'test_filters_1',
 'query_polygon': 'POLYGON ((-82.74809573102132 36.96082629937069, '
                  '-85.0932989306133 35.63154639485496, -81.0987220521874 '
                  '33.56697226279766, -79.4235769096217 36.34054727735634, '
                  '-79.4235769096217 36.34054727735634, -82.74809573102132 '
                  '36.96082629937069))',
 'sampling_protocols_omit': '',
 'use_taxon_geometry': True,
 'years_range': '2000,2020'}

import sqlite3
import pprint
import json
import pandas as pd
import numpy as np
pd.set_option('display.width', 600)
pd.set_option('display.max_colwidth', 30)
pd.set_option('display.max_rows', 150)
from IPython.display import Image
from pygbif import occurrences
import matplotlib.pyplot as plt
import os
from datetime import datetime
t1 = datetime.now()
import sys
sys.path.append(configDir)
import wranglerconfig as config
codeDir = config.codeDir
sys.path.append(codeDir)
import wrangler_functions as functions
working_directory = config.workDir
output_db = working_directory + query_name + '.sqlite'
username = config.gbif_username
password = config.gbif_password
email = config.gbif_email
EBD_file = config.EBD_file
print("Notebook run " + str(t1))
print(output_db)


gbif_data = get_GBIF_records(taxon_info, filter_set, working_directory, username, password, email)
gbif_data.to_csv("T:/Temp/gbif.csv")
