# The Wildlife Wrangler
Maintained by Nathan M. Tarr (github.com/nmtarr)

## Purpose
The abundance of wildlife occurrence datasets that are currently accessible can be valuable for efforts such as species distribution modeling and range delineation.  However, the task of downloading and cleaning occurrence records is often complex due to errors and uncertainties that are present in datasets.  This repository provides a framework for collecting and filtering occurrence data that are freely available through API's and data requests.  

## Framework
Data is requested from occurrence dataset API's and filtered according to species- and request-specific parameters.  Filtered occurrence records are saved in a database.  The details of species-concepts and filter parameter sets are stored in a database for use and reference.  Additionally, Jupyter Notebook documents are created that describe the filtered datasets for the sake of documentation and for decision making about filter parameterization.

## Features
This framework is designed to have certain features that provide summaries that can be interpreted at face value with high confidence or trusted when feeding into analyses and evaluations.
* __Automation__ -- The volume of data and species involved necessitates the processes be automated. Automation also reduces subjectivity in decision making, enables thorough documentation, and ensures repeatability.  However, some aspects of data wrangling and cleaning are unavoidably analog.

* __Detailed parameterization__ -- Data requests and record filters can be parameterized on a per-species and per-event basis. Rules do not have to be applied generally across large numbers of species or evaluations.

* __Open source__ -- Processes are coded in Python and SQL and use sqlite3, a built-in Python package, with an extension for spatial queries.

* __Transparency__ -- Summaries of occurrence data and models using empirical data include subjectivity in the form of choices regarding parameters and rules for handling/filtering/cleaning data.  This framework is meant to provide a way to document those choices.  

* __Geospatial processing__ -- Some geospatial operations are performed including re-projecting and buffering points.  A shapefile of buffered points is also created as output.

* __Data summaries__ -- Summaries of the attributes of occurrence record datasets that are returned by queries are created.

* __Spatial filters__ -- Queries can be limited to within geometries (polygons), and spatial restrictions can be assigned to species concepts (i.e., extent of occurrence polygon).  The user can also specify a continent and/or country within which to return records.

* __Duplicates__ -- Queries commonly include duplicates based on the latitude, longitude, and eventDate fields.  The user can opt to keep or exclude duplicates.  If they choose to exclude them, a multi-step process is triggered to account for two major issues.  One, the values of latitude and longitude for a record may have different numbers of digits to the right of the decimal (i.e., lat has eight decimals and lon has six).  Two, not all records have the same number of digits to the right of the decimal for latitude and longitude (i.e, one record may have two for lat and long while another has 12).  The process used is as follows: 1) latitude and longitude values of each record are truncated to the shorter of the two in cases where they differ, 2) if duplicates occur after that step, then the one with the largest individual count is kept, or the first if individual counts are the same, 3) records are identified that are a duplicate of a record with higher precision (e.g. (10.123, -10.123) would be flagged as a duplicate of (10.1234, -10.1234)) and then 4) they are removed.

* __Sources__ -- GBIF aggregates records from many datasets, collections, and institutions.  The user can specify collections and institutions to omit.  

* __Wildlife-centric__ -- The framework addresses several filter criteria that are especially relevant to species-occurrence records for studies of wildlife distributions and habitat-associations.
  * _Occurrence year_ -- Species' distributions and taxonomy can change over time, which warrants the ability to select records from within user-defined time periods.
  * _Occurrence month_ -- Relevant for investigations of migratory species and individual seasons.
  * _Coordinate uncertainty and issues_ -- Geographic coordinates vary immensely among records and some records have known issues that limit their value.  The framework allows the user to exclude records on the basis of reported issues and coordinate uncertainty.  It is also possible to filter out records that do not have an associated coordinate.
  * _Basis of record and sampling protocols_ -- Datasets such as GBIF include a variety of types of records, such as preserved specimens and fossil records.  Additionally, different sampling protocols may have been employed.  The user can choose which types to filter out.
  * _Detection distance_ -- Detection distance is the distance between an observer and an individual recorded.  It can affect the utility and potential scale of analyses of data because it adds to the locational uncertainty resulting from gps precision and observer movement during transects and traveling counts.  Different taxa may be sampled with different methods, and the uncertainty surrounding the exact locations of individuals recorded can vary among methods.  For example, small mammals that are captured in traps can confidently be assigned to the trap's location, but loud-singing birds detected in an auditory survey could be hundreds of meters away from the observer and thus, the coordinate associated with the record.  Some researchers choose not to address this, while others structure analyses around it.  The framework allows for either strategy, and anything in between, and the decision is documented. A default max detection distance can be specified for the species-concept and then over-ridden during query execution.  
  * _Locational uncertainty of recorded individuals_ -- As mentioned in *detection distance* and *coordinate uncertainty and issues*, although records are recorded as x,y coordinates, there are varying degrees of uncertainty regarding the exact locations of individuals recorded ("locational uncertainty").  Locational uncertainty is the sum of gps precision and the maximum possible detection distance of the species during the survey or sampling event.  This information is rarely attributed to records, so researchers must make assumptions or best guesses.  The framework provides a means of documenting and explaining those choices.  Additionally, points are buffered according to locational uncertainty of each record: the buffer radius is the reported or assumed coordinate uncertainty plus a user-provided value for maximum detection distance.  Additionally, records from transects or traveling counts may need their length added to given or assumed coordinate uncertainty from gps or other source of geolocation.   
  * _Occurrence data sets are dynamic_ -- Some datasets enable data contributors to go back and edit attributes of occurrence records.  In addition, historic records may be added that change the set of records associated with a past time period.  That is to say that a query of years past run today may be different than the same query run tomorrow.  This represents a challenge for provenance.  The method to handle this here is to document data request parameters and post-request filter sets as uniquely identifiable objects that are stored and documented in a database ('wildlife-wrangler.sqlite').  Records are linked to the filter sets used to acquire them.
  * _Instability of taxonomic lists and species concepts_ -- Taxonomic classifications are constantly being scrutinized and are revised annually.  As a result, separate projects may define a species differently.  In many cases, this is not a big problem because the species is unique and easily identifiable and taxonomic changes regard names only.  More problematic are cases where genetic studies have identified species that are nearly identical physiologically and were once considered a single species.  Such cases often reveal geographic patterns in species occurrence that may be used by some observers as a basis for identification of individuals in the field, thus introducing circularity. For example, a bird watcher decides an individual's identity in part based on which species is supposed to occur in the area *or* a species' range is revised and eBird changes records of the old species concept from the area where it is now known not to occur to the correct species.  Dynamic species-concepts are a daunting challenge, but can hopefully be handled with scrutinizing taxonomic matches.  In this framework, species-definitions are documented in a table within 'wildlife-wrangler.sqlite' with columns for the years the concept was valid and a geometry column where a polygon of potential occurrence can be recorded as Well-Known Text (WKT).  At present, the user must do some work to fill out the species-concepts table with species of interest before running queries.
  * _Data sensitivity regarding poaching risk_ -- More a challenge than a feature is the issue of occurrence records for many species, especially herpetofauna, being "sensitive data" due to the threat of poaching.  Individuals of rare and imperiled species are regularly captured in the wild for sale on black markets, sometimes providing hundreds of dollars per individual.  Poachers are known to review and interpret scientific data to determine the locations of populations that they can collect from.  This creates a serious challenge for those managing and summarizing locational data.  Managers of some databases "fuzz" or buffer points to coarsen the information to a level that isn't useful to poachers, which also limits the usefulness of data for conservation assessments.  Not only does this issue create limitations on the user's end, it also creates restrictions on how data providers can/should share and store records.  Data in the hands of federal agencies is subject to Freedom of Information Act requests.  This issue may shape aspects of the framework.  For example, it may be necessary to make this technology deployable to individual users and include a portal for data sets that are available to the user but remain unavailable from public facing or federal data sets.


## Recent changes (April 10, 2020)
* Added ability to limit requests to within geometries.
* Added ability to specify a limiting polygon for a species concept.
* Occurrence record database (output) now includes column for weight that users can use to omit or devalue undesirable records.
* Handling requests of > 200,000 records.
* Changed data filtering process from python dictionary based to pandas dataframe based to improve speed and shorten code.
* The older albers projection (EPSG: 102008) abandoned in favor of EPSG: 5070.
* Generalization of framework to better facilitate multiple users.  Parameters.sqlite replaced with a template (wildlife-wrangler_TEMPLATE.sqlite) that user can build from.  Each user will access a local copy of wildlife-wrangler.sqlite filled out for their needs.
* Added config file template.  Config file is necessary to avoid sharing your email password when performing large queries of GBIF.
* Species concept and filter sets are now documented in output databases.
* Improved handling of duplicates (see section on duplicates above for info).


## Coming soon
* Ability to incorporate bird records directly from a copy of the eBird EBD that user has downloaded.
* Making species level geometry filtering optional if polygon is present in species concepts table.
* Incorporating GBIF fields "dataGeneralizations", "georeferenceRemarks", and "informationWitheld".
* Overriding polygon geometry columns in output database is a "footprintWKT" value was provided.
* Incorporating species definition start and end dates.

## Inputs
Data is gathered from catalogs and databases via API's, so there are few inputs.  However, the 'wildlife-wrangler.sqlite' database is needed, which includes tables for species definitions, data request parameters, and post-request filtering criteria.

GBIF is currently the only dataset currently used but others can/will be added later including eBird.

## Outputs
On a per-species, per-query basis
* A database of filtered species-occurrence records with documentation.  The format supports display in QGIS and other GIS.
* Notebook documents that describe decisions made and the data acquired.

## Constraints
* Queries returning > 5,000,000 records may fail.
* Setup of spatialite can be difficult.
* Thorough and accurate specification of species definition is difficult and very time-consuming.
* Processing speed is limited in some cases by lack of spatial indexing, because setup of spatialite with spatial indexing enabled is very difficult.

## Dependencies
Python 3 and numerous packages including sqlite3 with the spatialite extension are needed.  Running the following code in a conda shell should create a suitable environment named "wrangler":
1. "conda create -n wrangler python=3.6 pandas jupyter basemap-data-hires notebook numpy shapely"
2. "conda activate wrangler"
3. "pip install pygbif python-dwca-reader sciencebasepy"

See the included spatialite install notes for how to install spatialite on windows 10.  

## Code
All code is included in this repository.  Runtimes of discrete tasks made grouping code into separate functions preferable.  

## USGS Approved Software Disclaimer
This software has been approved for release by the U.S. Geological Survey (USGS).
Although the software has been subjected to rigorous review, the USGS reserves
the right to update the software as needed pursuant to further analysis and
review. No warranty, expressed or implied, is made by the USGS or the
U.S. Government as to the functionality of the software and related material
nor shall the fact of release constitute any such warranty. Furthermore, the
software is released on condition that neither the USGS nor the U.S. Government
shall be held liable for any damages resulting from its authorized
or unauthorized use.
