## Installation on PC
1.  Use [Git](https://git-scm.com/) to clone the master repo from github to your computer.
2.  Build a [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment. Running the following lines of code (individually) in a conda shell should create a suitable environment named "wrangler":
    a. "conda create -n wrangler python=3.6 pandas jupyter basemap-data-hires notebook numpy shapely matplotlib=3.2.2"
    b. "conda activate wrangler"
    c. "pip install pygbif python-dwca-reader"
3.  Install [Spatialite](https://www.gaia-gis.it/fossil/libspatialite/index).  Spatialite is a spatial extension for [SQLite](https://www.sqlite.org/index.html).  
    SQLite is included in Python, but not the spatialite extension.  See the included file "spatialite PC install notes.txt" for more instructions.
4.  Copy the wildlife-wrangler_TEMPLATE.sqlite file to your computer, and
    delete "TEMPLATE" from the copied file's name ("wildlife-wrangler.sqlite").
5.  Copy the wranglerconfig_TEMPLATE.txt file to your computer and delete
    "TEMPLATE" from the copied file's name ("wranglerconfig.txt").  Fill out
    this file.  Create a folders named "Input" and "Output"
    within your designated working directory ("workDir").  Use "/" at the
    end of paths (i.e. "C:/Data/").
6.  Set up an account on GBIF.org and enter your credentials into
    your copy of wranglerconfig.txt.
7.  Download [DB Browser for SQLite](https://sqlitebrowser.org/) or
    another application for interacting with SQLite database.
8.  Run one of the notebook documents in the "Tests" folder.

## Using the Wildlife Wrangler
The Wildlife Wrangler is a set of tools to facilitate queries and filtering of species occurrence records.  While it automates many tasks, there is still a bit of work that must be done by the user.

#### General Overview of the Workflow
The user enters species information and unique sets of filtering parameters in their copy of wildlife-rangler.sqlite, then fills out a small portion of a [Jupyter Notebook](https://jupyter.org/) document with codes for a taxon concept and filter parameter sets to use and relevant directories.  When the notebook document is run, code that is stored in "wrangler_functions.py" and the notebook document itself retrieves records from GBIF, filters out unsuitable records, creates an output database where suitable records are stored along with documentation and summaries of record attributes before and after filtering, and performs various summaries of data attributes within the notebook document.  Thus, the primary results of running the wrangler are 1) the notebook document with documentation and data summaries and 2) the output (SQLite) database containing suitable records.  

#### Key Components of the Framework
*  __wildlife-wrangler.sqlite__ -- a centralized place to store filtering criteria and species definitions.  Saving filter sets (criteria) and taxa concepts as unique items in a database makes it much easier to explore combinations of species definitions and filtering criteria.  For example, if you want to use the same criteria for 20 species, you can call the same criteria each of the 20 times with just the codes.  This setup was chosen with the running of hundreds of queries over time in mind.
*  __report.ipynb__ -- this is where you control/run the wrangler.  It's kind of like a beefed-up form and report all in one.  Once you have species definitions and filter sets entered into wildlife-wrangler.sqlite, you can copy report notebook documents to create and run occurrence record queries/requests/downloads.
*  __wranglerconfig.py__ -- this is a .py file where you store some personal information that you wouldn't want saved in the notebook document: your email address and password for your GBIF account, which is needed in order to requests large downloads from GBIF.
*  __wrangler_functions.py__ -- a python module containing the meat and and potatoes functions of the wrangler.  DO NOT CHANGE!  Much of the necessary code is kept here to avoid having a thousand lines of code in the report.ipynb.  You can call some functions from this by importing the module in ipython (i.e., "import wrangler_functions as wranglers").  That can be handy for using the "getGBIFcode" function.

#### Detailed Instructions
1.  Open your copy of "__widlife-wrangler.sqlite__".  
2.  In the "__taxa_concepts__" table, enter in a taxon by entering a unique species code of your choosing in "taxon_id", the corresponding gbif species id in "gbif_id", "common_name", and "scientific_name".  GBIF species id codes can be retrieved from their website, or with the "getGBIFcode" function available in wildlife_functions.py.  Finally, enter an estimate of the maximum distance from an observer at which a species could be detected with common surveying methods in "detection_distance_meters".  For example, a shrew should have a value close to 0, whereas a wren should have a value closer to 100 m.  All other fields in this table may be helpful, but are not required.  
3.  In the "__gbif_requests__" table, enter a unique code for your filter set in "request_id".  Fill out each field, but note many fields have defaults.    These will be criteria for a first stage of filtering: only records meeting these criteria will be requested from the API.
4.  In the "__gbif_filters__" table, enter a unique code for your post-request filter set in "fiter_id". Fill out all other fields but note that defaults are present for some fields.  
5.  Copy "__report_TEMPLATE.ipynb__" to a location outside of the wrangler repo, say to your project directory.  Rename the notebook document to whatever you like.  Using a name with the species code request_id, and filter_id is helpful.
6.  In conda, activate your wrangler environment.  Open Jupyter Notebook and navigate to your renamed copy of "report_TEMPLATE.ipynb".  Fill out the first two cells of the notebook document and run it.  Run time can range from a few seconds to several hours.
7.  When you have completed a query/notebook document, you can export the notebook document as an html file and archive it for reference later.  

#### Where to Find Help
*  Within "wildlife-wrangler.sqlite", all tables and columns are explained in the "table_definitions" and "column_definitions" tables.  Example entries are included for "taxa_concepts", "gbif_requests", and "gbif_filters".
*  Example report notebook documents are provided in the examples folder.
*  The wrangler uses the [pygbif](https://pygbif.readthedocs.io/en/latest/) package, so its documentation can help explain the request step. These are the GBIF fields currently used to answer key questions about records:
   * What? -- "id", "gbifID", "individualCount", "identificationQualifier"

   * When? -- "eventDate", "retrievalDate"

   * Where? -- "coordinateUncertaintyInMeters", "decimalLatitude", "decimalLongitude", "footprintWKT", "geodeticDatum"

   * Who provided? -- "collectionCode", "institutionCode", "datasetName"

   * How obtained? -- "basisOfRecord", "samplingProtocol", "establishmentMeans",

   * Issues, notes, comments -- "issue" or "issues", "locality", "eventRemarks", "locationRemarks", "occurrenceRemarks"
