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
The user begins by creating a copy of the __Query_Form_TEMPLATE.ipynb__, which they then customize by adding a taxon concept, filter parameter sets to use, and paths to relevant directories.  When the customized notebook document is run, code that is stored in __"wrangler_functions.py"__ and the notebook document itself retrieves records from the Global Biodiversity Information Facility ([GBIF](https://gbif.org)) and/or eBird Basic Dataset ([eBird](https://ebird.org/science/use-ebird-data/)), filters out unsuitable records, creates an output database where suitable records are stored along with documentation and summaries of record attributes before and after filtering.  Additionally, various summaries of data attributes within the notebook document are performed.  Thus, the primary results of running the code are 1) the notebook document with documentation and data summaries and 2) the output SQLite database containing suitable records.  

#### Key Components of the Framework
*  __Filter sets__ -- Saving filter sets (criteria) and taxa concepts as unique items (JSON files) makes it much easier to explore combinations of taxa definitions and filtering criteria.  For example, if you want to use the same criteria for 20 species, you can call the same criteria each of the 20 times with just the codes.  This setup was chosen with the running of hundreds of queries over time in mind.
*  __Report.ipynb__ -- this is where you control/run the wrangler.  The file is a combined form and report all in one.  You can copy the report notebook document, rename it, fill it out, and run it to perform queries/downloads.
*  __wranglerconfig.py__ -- this is a .py file where you store some personal information that you wouldn't want saved in the notebook document: your email address and password for your GBIF account, which is needed in order to request downloads from GBIF.  A file path to an eBird Basic Dataset is saved here as well.
*  __wrangler_functions.py__ -- a python module containing the core functions of the wrangler.  DO NOT CHANGE THIS!  Much of the necessary code is kept here to avoid having a thousand lines of code in the report.ipynb.  You can call some helper functions from this by importing the module in ipython (i.e., "import wrangler_functions as wranglers").  That can be handy for using the "get_GBIF_code" and "generate_shapefile" functions.

#### Important Processes
__Incorporating Locational Uncertainty__ -- Species occurrence records represent events that occurred at specific locations, and thus, they have spatial extents and positions.  The geographical boundaries of observation events are defined by the observation method used and/or the observer’s behavior.  For spatial analyses and summaries, records must be assigned to geographic extents, and that process is referred to as georeferencing.  The geographic extents of observation events vary in size and shape and are rarely precisely mapped.  Usually, they must be approximated with information that has been associated with a record, such as spatial coordinates and text describing a location (Chapman and Wieczorek 2020).  

  Given that the geographic boundaries of many occurrence records in GBIF and eBird are not precisely mapped, we have to identify strategies for approximating them, as well as the sources of uncertainty.  In Chapman and Wieczorek’s (2020) description of best practices for georeferencing, they identify methods for georeferencing and the components of geographical uncertainties regarding the locations of observations.  The Wildlife Wrangler employs logic and methods that are meant to match the terminology, concepts, and methods presented by Chapman and Wieczorek (2020) to the best extent possible, with the goal of identifying or approximating a spatial representation of a record's entire location, including all uncertainties involved.  

Chapman and Wieczorek (2020) identify several key components for georeferencing that are addressed by the Wildlife Wrangler.

  * Georeferencing method -- Records are most commonly georeferenced by identifying a buffer radius length that, when applied to a record's coordinates, includes the entirety of the spatial extent of the event ("point-radius method").  However, the spatial extents of events are sometimes georeferenced as polygons instead ("shape method"), and such polygons can be stored in the Darwin Core "footprintWKT" field.  
  * Coordinate precision -- The Darwin Core schema includes a field for this named "coordinatePrecision", but it is rarely provided for many species.  However, the limits in precision of a pair of geographic coordinates can be calculated based upon the number of digits right of the decimal and their position on the WGS84 datum.  That precision may not agree with reported "coordinatePrecision" or "coordinateUncertaintyInMeters".  The Wildlife Wrangler calculates the best possible precision for a given set of coordinates and stores it in a field named "nominal_xy_precision".  Missing values in "coordinatePrecision" are replaced with 1.1 m.
  * GPS accuracy -- The accuracy and precision of coordinates influence locational uncertainty.  Low precision limits accuracy, but high precision does not guarantee high accuracy.  Therefore, GPS accuracy is important to consider when precision is high, which is the case for many records.  Many GPS units report a value they call accuracy, but it is not true accuracy.  True accuracy can only be determined by comparing reported positions to truth.  Accuracy is not included in GBIF/Darwin Core, and thus must be incorporated as an assumption.  Chapman and Wieczorek (2020) discuss this issue and recommend an assumed accuracy of 30 m for readings after 2000 and 100 m for readings from before 2000.  The Wrangler uses those as default accuracy values in a field named "gps_accuracy_m".
  * Distance traveled by observer -- Some occurrence records are derived from transects or traveling counts and the distance covered by the observer during the observation event adds to uncertainty about the location of individuals recorded.  This information is not provided by GBIF, but the eBird Basic Dataset does include it, and the Wildlife Wrangler stores it in a field named "effort_distance_m".  EBird Observation Dataset records, which are retrieved from GBIF, are assigned a default value of 8,407 m (5 miles) because most eBird records are from traveling counts that can be up to that length under the eBird Best Practices.  Zero is assigned to GBIF records.   
  * Detection distance -- Many species, especially birds, are observed at considerable distances from the observer's recorded location.  The maximum amount of uncertainty from this component can differ across species.  Detection distance is not provided, so the user must rely on an estimate of the maximum detection distance for a study species.  That value is stored in the "detection_distance_m" field.
  * Coordinate uncertainty -- For the point-radius method, the appropriate buffer radius for a record is stored in the Darwin Core field named "coordinateUncertaintyInMeters".  A well-georeferenced record would have a coordinate uncertainty value that encompasses all of the previously listed components.  In reality, this field is often empty and the user is left to approximate or estimate a value.  The user can opt to provide a standard, default coordinate uncertainty value to use when "coordinateUncertaintyInMeters" is empty.


The following processes are relevant to incorporating and accounting for uncertainties in the locations of recorded individuals.  

  1. Decimal degrees values are truncated to five or fewer digits because few GPS devices have precision that truly exceeds that.
  2. When "footprintWKT" is provided, the shape method is used and the record is mapped with the polygon that it delineates. Otherwise, the point-radius method is used.
  3. When the point-radius method must be used and a "coordinateUncertaintyInMeters" value is provided, that value is used for the buffer radius.  When no "coordinateUncertaintyInMeters" value is provided, a point-buffer radius is approximated by summing values from other fields and is stored in the "radius_m" field.  Which fields are used to calculate the radius length depends upon the data source.  The following table summarizes which fields are summed for each data source.

  | Source | coordinateUncertaintyInMeters | GPS_accuracy_m | effort_distance_m | detection_distance_m |
  | ------ | ----------------------------- | -------------- | ----------------- | -------------------- |
  | GBIF w/ uncertainty | X |  |  |  |
  | GBIF w/o default uncertainty | | X | X (0) | X |
  | GBIF w/ default uncertainty | X | | | |
  | EBD | | X | X | X |
  | EOD (from GBIF) | | X | X (8,047)| X |

  4. When the nominal coordinate precision is larger than the approximated buffer radius, the radius is replaced with the precision.  When "coordinatePrecision" is greater than both of those values, it is used instead.   

In summary, mapping occurrence records involves approximating non-point geometries, and the relevant information is often not provided.  However, the Wildlife Wrangler enables a user to customize approximations in broad strokes that are explicitly documented.  Choices in filter and taxon concept parameters determine how records will be mapped.  Approximations of locational uncertatinty can be fine tuned by adjusting relevant columns for individual records in the output .sqlite database.

__Removing Duplicate Records__ -- Queries commonly include duplicates based on the latitude, longitude, and date fields.  The user can opt to keep or exclude duplicates.  If they choose to exclude them, a multi-step process is triggered to account for two major issues.  One, the values of latitude and longitude for a record may have different numbers of digits to the right of the decimal (i.e., latitude has eight decimals but longitude has six).  Two, not all records have the same number of digits to the right of the decimal for latitude and longitude (i.e., one record may have two for latitude and longitude while another has 12).  The process used is as follows:
  1. Latitude and longitude values of each record are truncated to the shorter of the two in cases where they differ.
  2. If duplicates occur after that step, then the one with the largest individual count is kept, or the sorted first record if individual counts are the same.
  3. Records are identified that are a duplicate of a record with higher precision (e.g. (10.123, -10.123) would be flagged as a duplicate of (10.1234, -10.1234)).
  4. Duplicates are again removed.

#### Detailed Instructions
1.  Copy "__Query_Form_TEMPLATE.ipynb__" to a location outside of the wrangler repo, say to your project directory.  Rename the notebook document to whatever you like.
2.  In conda, activate your wrangler environment.  Open Jupyter Notebook and navigate to your renamed copy of __"Query_Form_TEMPLATE.ipynb"__.  
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

## References
Chapman, A.D. & Wieczorek, J.R. (2020) Georeferencing Best Practices. Copenhagen: GBIF Secretariat. https://doi.org/10.15468/doc-gg7h-s853
