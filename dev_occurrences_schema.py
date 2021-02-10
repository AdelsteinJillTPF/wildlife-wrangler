'''
Dev script for work on designing a schema for occurrences table that is created

Preserved Darwin Core names as much as possible
'''
spdb = "T:/temp/build_output_db_test.sqlite"
# def build_output_db(spdb)
"""
Description: Create a database for storing occurrence and taxon concept
data.
"""
spdb = spdb
# Delete the database if it already exists
if os.path.exists(spdb):
    os.remove(spdb)

# Create or connect to the database
conn = sqlite3.connect(spdb, isolation_level='DEFERRED')
#conn.enable_load_extension(True)
#conn.execute('SELECT load_extension("mod_spatialite")')
cursor = conn.cursor()

## Make database spatial
#conn.executescript('''SELECT InitSpatialMetaData(1);''')
conn.commit()

################################################# Create tables
sql_cdb = """
        /* Create a table for occurrence records, WITH GEOMETRY */
        CREATE TABLE IF NOT EXISTS occurrences (
                occ_id INTEGER NOT NULL PRIMARY KEY,
                taxon_id INTEGER NOT NULL,
                basisOfRecord TEXT,
                issues TEXT,
                collectionCode TEXT,
                institutionCode TEXT,
                datasetName TEXT,
                identificationQualifier TEXT,
                source TEXT NOT NULL,
                request_id TEXT NOT NULL,
                filter_id TEXT NOT NULL,
                latitude TEXT,
                longitude TEXT,
                coordinateUncertaintyInMeters INTEGER,
                occurrenceDate TEXT,
                retrievalDate TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
                individualCount INTEGER DEFAULT 1,
                dataGeneralizations TEXT,
                remarks TEXT,
                detection_distance INTEGER,
                radius_meters INTEGER,
                footprintWKT TEXT,
                weight INTEGER DEFAULT 10,
                weight_notes TEXT,
                GBIF_download_doi TEXT,
                    FOREIGN KEY (taxon_id) REFERENCES taxa(taxon_id)
                    ON UPDATE RESTRICT
                    ON DELETE NO ACTION);
"""
cursor.executescript(sql_cdb)
