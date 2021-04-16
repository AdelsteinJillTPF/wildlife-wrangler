## Installation for Windows 10
1.  Use [Git](https://git-scm.com/) to clone the master repo from github to your computer.
2.  Build a [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment. Running the following lines of code (individually) in a conda shell should create a suitable environment named "wrangler":
    a. "conda create -n wrangler python=3.6 ipython rpy2 pandas jupyter notebook numpy shapely matplotlib r-base geopandas descartes"
    b. "conda activate wrangler"
    c. "pip install pygbif python-dwca-reader"
3.  Copy the __wranglerconfig_TEMPLATE.txt__ file to your     computer and delete "TEMPLATE" from the copied file's name ("__wranglerconfig.txt__").  Fill out this file.  Create a folders named "Input" and "Output" within your designated working directory ("workDir").  Use "/" at the end of paths (i.e. "C:/Data/").
4.  Set up an account on GBIF.org and enter your credentials into your copy of wranglerconfig.txt.
5.  Download [DB Browser for SQLite](https://sqlitebrowser.org/) or another application for interacting with SQLite database.
6.  Run one of the notebook documents in the "Tests" folder.

## Using the Wildlife Wrangler
The Wildlife Wrangler ("the wrangler") is a set of tools to facilitate queries and filtering of species occurrence records.  While it automates many tasks, there is still a bit of work that must be done by the user.

#### General Overview of the Workflow
The user begins by creating a copy of the __report_TEMPLATE.ipynb__, which they then customize by adding a taxon concept, filter parameter sets to use, and paths to relevant directories.  When the customized notebook document is run, code that is stored in __"wrangler_functions.py"__ and the notebook document itself retrieves records from the Global Biodiversity Information Facility ([GBIF](https://gbif.org)) and/or eBird Basic Dataset ([eBird](https://ebird.org/science/use-ebird-data/)), filters out unsuitable records, creates an output database where suitable records are stored along with documentation and summaries of record attributes before and after filtering.  Additionally, various summaries of data attributes within the notebook document are performed.  Thus, the primary results of running the code are 1) the notebook document with documentation and data summaries and 2) the output SQLite database containing suitable records.  

#### Key Components of the Framework
*  __Filter sets__ -- Saving filter sets (criteria) and taxa concepts as unique items (JSON files) makes it much easier to explore combinations of taxa definitions and filtering criteria.  For example, if you want to use the same criteria for 20 species, you can call the same criteria each of the 20 times with just the codes.  This setup was chosen with the running of hundreds of queries over time in mind.
*  __Report.ipynb__ -- this is where you control/run the wrangler.  The file is a combined form and report all in one.  You can copy the report notebook document, rename it, fill it out, and run it to perform queries/downloads.
*  __wranglerconfig.py__ -- this is a .py file where you store some personal information that you wouldn't want saved in the notebook document: your email address and password for your GBIF account, which is needed in order to request downloads from GBIF.  A file path to an eBird Basic Dataset is saved here as well.
*  __wrangler_functions.py__ -- a python module containing the core functions of the wrangler.  DO NOT CHANGE THIS!  Much of the necessary code is kept here to avoid having a thousand lines of code in the report.ipynb.  You can call some helper functions from this by importing the module in ipython (i.e., "import wrangler_functions as wranglers").  That can be handy for using the "get_GBIF_code" and "generate_shapefile" functions.

#### Detailed Instructions
1.  Copy "__report_TEMPLATE.ipynb__" to a location outside of the wrangler repo, say to your project directory.  Rename the notebook document to whatever you like.
2.  In conda, activate your wrangler environment.  Open Jupyter Notebook and navigate to your renamed copy of __"report_TEMPLATE.ipynb"__.  
3.  Fill out the notebook document and run it.  A JSON file of the taxon information and filter parameters ("filter set") will be saved when the notebook is run.  Alternatively, you can specify existing JSON's to use.  NOTE: Run time can range from a few seconds to several hours.
4.  When you have completed a query/notebook document, you can export the notebook document as an html file and archive it for reference later.  The html versions can be exported within Jupyter Notebook or by adapting and running the following command line code: jupyter nbconvert --to html --TemplateExporter.exclude_input=True --output-dir="C:/YourFolder/" NOTEBOOK_NAME_HERE.ipynb
5.  Open the output database in DB Browser or elsewhere and adjust the record weights as desired.

#### Where to Find Help
*  Example report notebook documents are provided in the examples folder.
*  The wrangler uses the [pygbif](https://pygbif.readthedocs.io/en/latest/) and [Auk]("https://cornelllabofornithology.github.io/auk/") packages, so their documentations can help explain some aspects. These are the GBIF fields currently used to answer key questions about records:
   * What? -- "id", "gbifID", "individualCount", "identificationQualifier"

   * When? -- "eventDate", "retrievalDate"

   * Where? -- "coordinateUncertaintyInMeters", "decimalLatitude", "decimalLongitude", "footprintWKT", "geodeticDatum"

   * Who provided? -- "collectionCode", "institutionID", "datasetName"

   * How obtained? -- "basisOfRecord", "samplingProtocol", "establishmentMeans", "source"

   * Issues, notes, comments -- "issue" or "issues", "locality", "eventRemarks", "locationRemarks", "occurrenceRemarks"
