filter_set = {'bases_omit': None,
             'collection_codes_omit': None,
             'country': 'US',
             'datasets_omit': None,
             'default_coordUncertainty': 1000,
             'duplicates_OK': False,
             'geoissue': None,
             'get_dwca': True,
             'has_coordinate_uncertainty': False,
             'institutions_omit': None,
             'issues_omit': None,
             'lat_range': '27,41',
             'lon_range': '-89,-75',
             'max_coordinate_uncertainty': 10000,
             'months_range': '5,6',
             'name': 'test_filters_1',
             'query_polygon': 'POLYGON ((-82.74809573102132 36.96082629937069, '
                              '-85.0932989306133 35.63154639485496, -81.0987220521874 '
                              '33.56697226279766, -79.4235769096217 36.34054727735634, '
                              '-79.4235769096217 36.34054727735634, -82.74809573102132 '
                              '36.96082629937069))',
             'sampling_protocols_omit': None,
             'use_taxon_geometry': False,
             'years_range': '2016,2017'}

taxon_info = {'EBIRD_ID': 'Yellow-billed Cuckoo',
                 'GBIF_ID': 2496287,
                 'ID': 'TestCuckoo',
                 'TAXON_EOO': 'POLYGON ((-84.09680233298448 36.69265225442667, '
                              '-84.07962135716329 34.5561660300382, -84.07962135716329 '
                              '34.5561660300382, -80.25685423694925 34.65515526072436, '
                              '-81.15026497965096 36.71331438415306, -84.09680233298448 '
                              '36.69265225442667))',
                 'detection_distance_m': 200}


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
config_path = "T:/Data/"
import sys
sys.path.append(config_path)
import wranglerconfig as config
sys.path.append(config.codeDir)
import wrangler_functions as functions

# Define some variables
working_directory = config.workDir
username = config.gbif_username
password = config.gbif_password
email = config.gbif_email
EBD_file = config.EBD_file

po0ly = "POLYGON ((-82.74809573102132 36.96082629937069, -85.0932989306133 35.63154639485496, -81.0987220521874 33.56697226279766, -79.4235769096217 36.34054727735634, -79.4235769096217 36.34054727735634, -82.74809573102132 36.96082629937069))"
po0ly2 = shapely.wkt.loads(po0ly)
po0ly2.exterior.is_ccw



from pygbif import occurrences # I'm using version 0.5.0

download_filters = ["taxonKey = 2496287",
                    "geometry = 'POLYGON ((-82.7 36.9, -85.0 35.6, -81.0 33.5, -79.4 36.3, -79.4 36.3, -82.7 36.9))'"]

download = occurrences.search(download_filters,
                                pred_type='and',
                                user = username,
                                pwd = password,
                                email = email)


from pygbif import utils
x = 'POLYGON ((-82.7 36.9, -85.0 35.6, -81.0 33.5, -79.4 36.3, -79.4 36.3, -82.7 36.9))'
utils.wkt_rewind(x)
