import pandas as pd
import numpy as np

df_unfiltered = pd.read_csv("T:/Temp/dev_acc_prec.csv", usecols=["decimalLatitude", "decimalLongitude"],
                 dtype={"decimalLatitude": str, "decimalLongitude": str})

# Trim decimal length to 5 digits (lat and long)
df_unfiltered["decimalLatitude"] = df_unfiltered["decimalLatitude"].apply(lambda x: x[:8])
df_unfiltered["decimalLongitude"] = df_unfiltered["decimalLongitude"].apply(lambda x: x[:9])

# Calculate the number of digits for latitude and longitude
df_unfiltered['digits_latitude'] = [len(x.split(".")[1]) for x in df_unfiltered['decimalLatitude']]
df_unfiltered['digits_longitude'] = [len(x.split(".")[1]) for x in df_unfiltered['decimalLongitude']]

# Remove records with fewer than 2 digits
df_unfiltered = df_unfiltered[(df_unfiltered["digits_latitude"] > 1) & (df_unfiltered["digits_longitude"] > 1)]

# Calculate nominal precisions (meters)
digitsX = {1: 10, 2: 100, 3: 1000, 4: 10000, 5: 100000}
df_unfiltered["temp"] = df_unfiltered["decimalLatitude"].apply(lambda x: (111321 * np.cos(float(x) * np.pi/180)))
df_unfiltered["temp2"] = df_unfiltered["digits_longitude"].apply(lambda x: digitsX[x])
df_unfiltered["nominal_x_precision"] = df_unfiltered["temp"]/df_unfiltered["temp2"]# decimal moved based on digits.

digitsY = {1: 11112.0: 1111.2, 3: 111.1, 4: 11.1, 5: 1.1} # Lookup for latitude precision
df_unfiltered["nominal_y_precision"] = df_unfiltered["digits_latitude"].apply(lambda x: digitsY[x])

# Put the larger of the two nominal precisions in a column
df_unfiltered["nominal_xy_precision"] = np.where(df_unfiltered["nominal_y_precision"] > df_unfiltered["nominal_x_precision"], df_unfiltered["nominal_y_precision"], df_unfiltered["nominal_x_precision"])

df_unfiltered.drop(["temp", "temp2", "digits_latitude", "digits_longitude", "nominal_x_precision", "nominal_y_precision"], axis=1, inplace=True)

print(df_unfiltered.head(30))




# Record df length before removing duplicates
initial_length = len(df)
#df = df.astype(output_schema)

################################################################################
################################################################################
################################################################################
"""
############ RECTIFY UNEQUAL LAT-LONG PRECISION
First, trim decimal length in cases where decimal length differs between
latitude and longitude values, result is equal latitude and longitude
length.  Record the trimmed decimal precision in a temp column for use later
as a record to "native" precision.
"""
df['dup_latPlaces'] = [len(x.split(".")[1]) for x in df['decimalLatitude']]
df['dup_lonPlaces'] = [len(x.split(".")[1]) for x in df['decimalLongitude']]
df['dup_OGprec'] = df['dup_latPlaces']
prec_unequal = df[df['dup_latPlaces'] != df['dup_lonPlaces']]
for i in prec_unequal.index:
    x = prec_unequal.loc[i]
    if x['dup_latPlaces'] < x['dup_lonPlaces']:
        trim_len = int(x['dup_latPlaces'])
    else:
        trim_len = int(x['dup_lonPlaces'])
    df.loc[i, 'decimalLatitude'] = x['decimalLatitude'][:trim_len + 3]
    df.loc[i, 'decimalLongitude'] = x['decimalLongitude'][:trim_len + 4]
    # Record the resulting precision for reference later
    df.loc[i, 'dup_OGprec'] = trim_len
df.drop(['dup_latPlaces', 'dup_lonPlaces'], axis=1, inplace=True)

"""
########  INITIAL DROP OF DUPLICATES
Initial drop of duplicates on 'latitude', 'longitude', 'eventDate',
keeping the first (highest individual count)
Sort so that the highest individual count is first ## ADD OCCURRENCEDATE BACK IN
"""
df.sort_values(by=['decimalLatitude', 'decimalLongitude', 'eventDate',
                    'individualCount'],
                ascending=False, inplace=True, kind='mergesort',
                na_position='last')

df.drop_duplicates(subset=['decimalLatitude', 'decimalLongitude', 'eventDate'],
                   keep='first', inplace=True)

"""
#########  FIND IMPRECISE DUPLICATES
Get a list of "native" precisions that are present in the data to loop through.
Next, iterate through this list collecting id's of records that need to be
removed from the main df.
"""
# Get list of unique precisions.  Order is important: descending.
precisions = list(set(df['dup_OGprec']))
precisions.sort(reverse=True)
# The highest precisions listed at this point has already been done: drop it.
precisions = precisions[1:]

# List for collecting records that are duplicates
duplis = []

# The precision-specific duplicate testing happens repeatedly, so make it a
# function.
def drop_duplicates(precision, df):
    """
    Function to find undesirable duplicates at a particular decimal precision.

    Parameters
    ----------
    precision : The level of precision (places right of decimal) in decimalLatitude
    and longitude values for the assessment of duplicates.
    df : dataframe to assess and drop duplicates from.  This function works
          'inplace'.
    """
    # Create a df with records from the input df having decimal precision > the
    # precision level being assessed.
    dfLonger = df[df['dup_OGprec'] > precision].copy()
    # Truncate lat and long values
    dfLonger['decimalLatitude'] = [x[:precision + 3] for x in dfLonger['decimalLatitude']]
    dfLonger['decimalLongitude'] = [x[:precision + 4] for x in dfLonger['decimalLongitude']]

    # Create a df with records having the precision being
    # investigated
    dfShorter1 = df[df['dup_OGprec'] == precision]

    # Find records in dfShorter1 with latitude, longitude, date combo
    # existing in dfLonger and append to list of duplis
    dfduplis = pd.merge(dfShorter1, dfLonger, how='inner',
                        on=['decimalLatitude', 'decimalLongitude', 'eventDate'])
    dups_ids = dfduplis['record_id_x']
    for d in dups_ids:
        duplis.append(d)

# Drop latitude longitude duplicates at lower decimal precisions
for p in precisions:
    drop_duplicates(p, df)

# Drop rows from the current main df that have been identified as duplicates.
df2 = df[df['record_id'].isin(duplis) == False].copy()

# Drop excess columns
df2.drop(columns=['dup_OGprec'], axis=1, inplace=True)
