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

target_df = pd.DataFrame(columns=['record_id', 'request_id', 'filter_id', 'source', 'retrieval_date',
       'detection_distance', 'radius_meters', 'GBIF_download_doi',
       'general_remarks', 'weight', 'weight_notes', 'accessRights',
       'bibliographicCitation', 'basisOfRecord', 'collectionCode',
       'coordinateUncertaintyInMeters', 'dataGeneralizations', 'datasetName',
       'decimalLatitude', 'decimalLongitude', 'establishmentMeans',
       'eventDate', 'eventRemarks', 'footprintWKT', 'footprintSRS',
       'georeferencedBy', 'georeferenceProtocol',
       'georeferenceVerificationStatus', 'georeferenceRemarks', 'habitat',
       'identifiedBy', 'identifiedRemarks', 'identificationQualifier',
       'individualCount', 'informationWitheld', 'institutionCode', 'issues',
       'license', 'locationAccordingTo', 'locationRemarks', 'modified',
       'occurrenceStatus', 'occurrenceRemarks', 'recordedBy',
       'samplingProtocol', 'taxonConceptID', 'verbatimLocality'])


def validate_occurrences_df(dataframe):
    """
    Checks whether a data frame matches the schema of occurrences in an output
    database occurrences table.

    PARAMETERS
    dataframe : input dataframe to assess the format of

    RETURNED
    True or False
    """
    # Are the columns names correct?

    # Are the data types of each column correct?
