'''

    Sample scripting to use pygbif's occurrence module download request.
    This will download a zip file of Darwin Core Archive files. The file
    "occurrence.txt" is the file with species records.
    Additionally, installing the dwca-reader package will facilitate
    reading the Darwin Core Archive occurrence file into a Pandas
    dataframe for further manipulation.

    You must have a GBIF account to access the downloading ability. Once
    you have established an account, you will need to use the username
    and password as parameters in the download function. These can either
    be supplied directly in the code, can be set in your user environment
    variable for your OS, or be created in a separate config file.

    PYTHON 3.6


'''

# Import modules
import os, shutil
from datetime import datetime
from pygbif import occurrences as occ
from pygbif import species
import pandas as pd
from dwca.read import DwCAReader

pd.set_option('display.max_columns', 10)
pd.set_option('display.max_rows', 100)

t0 = datetime.now()

workDir = 'C:/Data/USGS Analyses/GAP-Habitat-Map-Assessment/'
downDir = workDir + 'downloads/'

'''
# Make temporary directory for downloads
#  remove it if it already exists
if os.path.exists(downDir):
    shutil.rmtree(downDir)
    os.mkdir(downDir)
else:
    os.mkdir(downDir)
'''


sciName = 'Wilsonia citrina'  # Small number of records
#sciName = 'Buteo lagopus'  # >200k records

# First use the species module to get the taxonKey for a species scientific name
print("\nWorking on the following species:", sciName)
tkey = species.name_backbone(name = sciName, rank='species')['usageKey']

# Now use the download method to get a download key for that species
# It returns a dictionary of results containing request parameters
print("Getting the download key .....")
res = occ.download(['taxonKey = {0}'.format(tkey), 'hasCoordinate = TRUE', 'country = US'],
                    user='pythonprocessing@gmail.com',
                    pwd='@re!!',
                    email='mo@ncsu.edu')

# Get the value of the download key
dkey = res[0]

# Now download the actual zip file containing the Darwin Core files
'''
 NOTE:
     The download can take a while to generate and is not immediately
     available once the download_get command has been issued. Use a
     while and try loop to make sure the download has succeeded.
     The zipdownload variable will be a dictionary of the path,
     the file size, and the download key unique code. It can be used
     to change the file name, unzip the file, etc.
'''

print("Attempting to download the Darwin Core Archive zip file for this species .....")
gotit = None
while gotit is None:
    try:
        zipdownload = occ.download_get(key=dkey,path=downDir)
        gotit = 1
    except:
        pass

t1 = datetime.now()

print("\n\n+++++ Download time was", t1 - t0, '+++++')

'''
    Now try reading the "occurrence.txt" file into a Pandas dataframe
    using the DWCA Reader Python package.
    pip install python-dwca-reader

'''
print('\n' + '='*55)
with DwCAReader(downDir + dkey + '.zip') as dwca:
    print(' Reading occurrence records into Pandas dataframe ....')
    dfOcc = dwca.pd_read('occurrence.txt', parse_dates=True)

print('   There are', len(dfOcc), 'records for this species\n')

########## Start manipulating the occurrences dataframe ##########

## Get only a handful of columns
dfMap = dfOcc[["decimalLatitude", "decimalLongitude",
               "coordinateUncertaintyInMeters",
               "day", "month", "year"]]

## Get only unique x-y locations based on day, month, and year
# NOTE: Here is a way to count the number of missing records across all columns
#dfMap.isna().sum()
dfMapUnique = dfMap.drop_duplicates(["decimalLatitude", "decimalLongitude",
                                     "coordinateUncertaintyInMeters",
                                     "day", "month", "year"])

## Drop rows with null value for coordinate uncertainty
dfMap2 = dfMapUnique.dropna(subset=['coordinateUncertaintyInMeters'])

## Get only the records where coordinate uncertainty < 10k and not 0.0
dfMap3 = dfMap2[(dfMap2['coordinateUncertaintyInMeters']<10000) &
                (dfMap2['coordinateUncertaintyInMeters']!=0.0)]


# Remove the downloaded file
#print('!!! Deleting downloaded zip file !!!')
#os.remove(downDir + dkey + '.zip')
