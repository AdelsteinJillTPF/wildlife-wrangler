"""
N. Tarr
February 19, 2021

This script updates a copy of the wildlife-wrangler.sqlite so that it
accomodates incorporation of the eBird Basic Dataset.

A new table for eBird filtering parameters is added.
"""
import sqlite3

# Identify the parameter sqlite database to update.
db = "T:/code/wildlife-wrangler/wildlife-wrangler_TEMPLATE.sqlite"

# Create or connect to the database
conn = sqlite3.connect(db)

sql="""
/*Create a table for requests/queries of EBD records. ------------------------*/
CREATE TABLE IF NOT EXISTS "EBD_requests" (
"request_id"	TEXT NOT NULL UNIQUE,
"source"	TEXT DEFAULT 'EBD',
"start_date"   TEXT,
"end_date"  TEXT,
"country"  TEXT DEFAULT 'US',
"lat_range" TEXT,
"lon_range" TEXT,
"months_range" TEXT,
"max_distance" TEXT,
"creator"	TEXT,
"notes"	TEXT,
PRIMARY KEY("request_id")
);

/* Insert the example EBD request --------------------------------------------*/
INSERT INTO EBD_requests (request_id, start_date, end_date, lat_range,
                          lon_range, months_range, max_distance, creator, notes)
VALUES ("example1", "2015-01-01", "2020-12-31", "27,41", "-91,-75", "3,5", "10",
        "N. Tarr", "An example")
"""
conn.executescript(sql)

sql="""
/*Insert record for table description ----------------------------------------*/
INSERT INTO table_descriptions (table_name, description)
VALUES ("EBD_requests", "Stores and documents filter sets used during
queries of a copy of the the EBird Basic Dataset via the auk package
in R.")
"""
conn.executescript(sql)

sql="""
/*Insert column descriptions--------------------------------------------------*/
INSERT INTO column_descriptions (table_name, column_name, description)
VALUES
    ("EBD_requests", "request_id", "Unique code for the filter set"),
    ("EBD_requests", "source", "eBird"),
    ("EBD_requests", "start_date", "Earliest date to include.  Should be text
            formatted like 2015-02-23"),
    ("EBD_requests", "end_date", "Most recent date to include.  Should be
            text like 2021-04-11"),
    ("EBD_requests", "country", "Defaults to US, but could be excluded."),
    ("EBD_requests", "lat_range", "Latitude range like 27,41"),
    ("EBD_requests", "lon_range", "Longitude range like -90,-71"),
    ("EBD_requests", "months_range", "Range of months to include like 3,5"),
    ("EBD_requests", "max_distance", "Maximum allowable distance for a
            traveling count in kilometers."),
    ("EBD_requests", "creator", "Name of the creator"),
    ("EBD_requests", "notes", "Any relevant notes about the filter set.");
"""
conn.executescript(sql)
conn.close()


sql = """
DROP TABLE EBD_requests;

DELETE FROM table_descriptions WHERE table_name = "EBD_requests";

DELETE FROM column_descriptions WHERE table_name = "EBD_requests";
"""
# Create or connect to the database
conn = sqlite3.connect(db)
conn.executescript(sql)
conn.close()
