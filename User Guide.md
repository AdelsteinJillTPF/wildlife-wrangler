## Installation
1.  Use git to clone the master repo from github to you computer.
2.  Build a conda environment. Running the following code in a conda shell
    should create a suitable environment named "wrangler":
    a. "conda create -n wrangler python=3.6 pandas jupyter basemap-data-hires notebook numpy shapely"
    b. "conda activate wrangler"
    c. "pip install pygbif python-dwca-reader sciencebasepy"
3.  Install spatialite.  Spatialite is a spatial extension for SQLite.  SQLite
    is included in Python, but not the spatialite extension.  See the
    included file "spatialite PC install notes.txt" for more instructions.
4.  Copy the wildlife-wrangler_TEMPLATE.sqlite file to your computer, and
    delete "_TEMPLATE" from the copied file's name ("wildlife-wrangler.sqlite").
5.  Copy the wranglerconfig_TEMPLATE.txt file to your computer and delete
    "_TEMPLATE" from the copied file's name ("wranglerconfig.txt").
6.  Set up an account on GBIF.org and enter your credentials into
    your copy of wranglerconfig.txt.
7.  Download DB Browser for SQLite (https://sqlitebrowser.org/) or
    another application for interacting with SQLite database.

## Using the Wildlife Wrangler
The Wildlife Wrangler is a set of tools to facilitate queries and filtering of
species occurrence records.  While it automates many tasks, there is still
a bit of work that must be done by the user.

### General Overview of the Workflow
The user enters species concepts and unique sets of filtering parameters in their copy of wildlife-wrangler.sqlite, then fills out a small portion of a
Jupyter Notebook with a species concept of interest, filter parameters sets to use, and relevant directories.  When the notebook is
run, code that is stored in "wrangler_functions.py" and the notebook itself retrieves records from GBIF, filters out unsuitable records, creates an output database where
suitable records are stored along with documentation and summaries of record attributes before and after filtering, and performs various summaries of data attributes within the notebook.  Thus, the primary results of running the wrangler are 1) the notebook with documentation and data summaries and 2) the output (SQLite) database containing suitable records.  

### Detailed Instructions
1.  Open your copy of "__widlife-wrangler.sqlite__".  
2.  In the "__species_concepts__" table, enter in a species concept by
    entering a unique species code of your choosing in "species_id", the
    corresponding gbif species id in "gbif_id", "common_name", and "scientific_name".  GBIF species id codes can be retrieved from their
    website, or with the "getGBIFcode" function available in wildlife_functions.py.  Finally, enter an estimate of the maximum
    distance from an observer at which a species could be detected with common surveying methods in "detection_distance_meters".  For example, a shrew should have a value close to 0, whereas a wren should have a value closer to 100 m.  All other fields in this table may be helpful, but are not required.
3.  In the "__gbif_requests__" table, enter a unique code for your filter
    set in "request_id".  Fill out each field, but note many fields have defaults.  These will be criteria for a first stage of filtering: only records meeting these criteria will be requested from the API.
4.  In the "__gbif_filters__" table, enter a unique code for your
    post-request filter set in "fiter_id". Fill out all other fields but note that defaults are present for some fields.
5.  Copy "__report_TEMPLATE.ipynb__" to a location outside of the wrangler
    repo, say to your project directory.  Rename the notebook to whatever
    you like.  Using a name with the species code request_id, and filter_id is helpful.
6.  In conda, activate your wrangler environment.  Open Jupyter Notebook
    and navigate to your renamed copy of "report_TEMPLATE.ipynb".  Fill out
    the first two cells of the notebook and run it.  Run time can range from
    a few seconds to several hours.
7.  When you have completed a query/notebook, you can export the notebook as
    an html file and archive it for reference later.  

### Where to Find Help
*  The wrangler uses the [pygbif](https://pygbif.readthedocs.io/en/latest/)  package, so its documentation can help
   explain the request step.
*  Within "wildlife-wrangler.sqlite", all tables and columns are explained
   in the "table_definitions" and "column_definitions" tables.  EXample
   entries are included for "species_concepts", "gbif_requests", and "gbif_filters".
*  Example report notebooks are provided in the examples folder.
*  [Conda](https://docs.conda.io/projects/conda/en/latest/index.html)
*  [Jupyter Notebook](https://jupyter.org/)
*  [SQLite](https://www.sqlite.org/index.html)
*  [Spatialite](https://www.gaia-gis.it/fossil/libspatialite/index)
