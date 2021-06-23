# Import modules
import os, shutil
from datetime import datetime
from pygbif import occurrences as occ
from pygbif import species
import pandas as pd
from dwca.read import DwCAReader

t0 = datetime.now()

downDir = 'T:/Occurrence_Records/downloads/'

sciName = 'Wilsonia citrina'  # Small number of records
#sciName = 'Buteo lagopus'  # >200k records

# First use the species module to get the taxonKey for a species scientific name
print("\nWorking on the following species:", sciName)
tkey = species.name_backbone(name = sciName, rank='species')['usageKey']

# Now use the download method to get a download key for that species
# It returns a dictionary of results containing request parameters
print("Getting the download key .....")
res = occ.download(['taxonKey = {0}'.format(tkey), 'hasCoordinate = TRUE', 'country = US'],
                    user='gapper',
                    pwd='metspirates',
                    email='pythonprocessing@gmail.comu')

# Get the value of the download key
dkey = res[0]

# Now download the actual zip file containing the Darwin Core files
print("Attempting to download the Darwin Core Archive zip file for this species .....")

success = False
while success == False:
    try:
        zipdownload = occ.download_get(key=dkey,path=downDir)
    except:
        pass

    try:
        with DwCAReader(downDir + dkey + '.zip') as dwca:
            print(' Reading occurrence records into Pandas dataframe ....')
            dfOcc = dwca.pd_read('occurrence.txt', parse_dates=True)
            success = True
    except Exception as e:
        print(e)
        pass
print('   There are', len(dfOcc), 'records for this species\n')
