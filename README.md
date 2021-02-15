# The Wildlife Wrangler
__Maintainer__
Nathan M. Tarr (nmtarr@ncsu.edu, [orcid:0000-0003-2925-8948](https://orcid.org/0000-0003-2925-8948)), github.com/nmtarr)


__Contributors__
Alexa McKerrow (amckerrow@usgs.gov, [orcid:0000-0002-8312-2905](https://orcid.org/0000-0002-8312-2905))

Matthew Rubino (mjrubino@ncsu.edu, [orcid:0000-0003-0651-3053](https://orcid.org/0000-0003-0651-3053))

## USGS Software Release Information
IP-120038
https://doi.org/10.5066/P98K7E93


## Purpose
The abundance of wildlife occurrence datasets that are currently accessible can be valuable for efforts such as species distribution modeling and range delineation.  However, the task of downloading and cleaning occurrence records is often complex due to errors and uncertainties that are present in datasets.  This repository provides a framework for collecting and filtering occurrence data that are freely available through API's and data requests.  

## Framework
Data is requested from occurrence dataset API's and filtered according to species- and request-specific parameters.  Filtered occurrence records are saved in a database.  The details of taxa concepts and filter parameter sets are stored in a database for use and reference.  Additionally, Jupyter Notebook documents are created that describe the filtered datasets for the sake of documentation and for decision making about filter parameterization.

## Features
This framework is designed to have certain features that provide summaries that can be interpreted at face value with high confidence or trusted when feeding into analyses and evaluations.
* __Automation__ -- The volume of data and species involved necessitates the processes be automated. Automation also reduces subjectivity in decision making, enables thorough documentation, and ensures repeatability.  However, some aspects of data wrangling and cleaning are unavoidably analog.

* __Open source components__ -- Processes are coded in Python and SQL and use sqlite3, a built-in Python package, with an extension for spatial queries.

* __Detailed parameterization__ -- Data requests and record filters can be parameterized on a per-species and per-event basis. Rules do not have to be applied generally across large numbers of species or evaluations.

* __Transparency__ -- Summaries of occurrence data and models using empirical data include subjectivity in the form of choices regarding parameters and rules for handling/filtering/cleaning data.  This framework is meant to provide a way to document those choices.  

* __Geospatial processing__ -- Some geospatial operations are performed including re-projecting and buffering points.  A shapefile of buffered points is also created as output.

* __Data summarization__ -- Summaries of the attributes of occurrence record datasets that are returned by queries are created.

* __Spatial filtering__ -- Queries can be limited to within geometries (polygons), and spatial restrictions can be assigned to species definitions (i.e., extent of occurrence polygon).  The user can also specify a continent and/or country within which to return records.

* __Source filtering__ -- Global Biodiversity Information Facility ([GBIF](https://gbif.org)) aggregates records from many datasets, collections, and institutions.  The user can specify collections and institutions to omit.  

* __Duplicate handling__ -- Queries commonly include duplicates based on the latitude, longitude, and eventDate fields.  The user can opt to keep or exclude duplicates.  If they choose to exclude them, a multi-step process is triggered to account for two major issues.  One, the values of latitude and longitude for a record may have different numbers of digits to the right of the decimal (i.e., lat has eight decimals and lon has six).  Two, not all records have the same number of digits to the right of the decimal for latitude and longitude (i.e., one record may have two for lat and long while another has 12).  The process used is as follows: 1) latitude and longitude values of each record are truncated to the shorter of the two in cases where they differ, 2) if duplicates occur after that step, then the one with the largest individual count is kept, or the first if individual counts are the same, 3) records are identified that are a duplicate of a record with higher precision (e.g. (10.123, -10.123) would be flagged as a duplicate of (10.1234, -10.1234)) and then 4) they are removed.

* __Wildlife-centric design__ -- The framework addresses several filter criteria that are especially relevant to species occurrence records for studies of wildlife distributions and habitat associations.
  * _Occurrence year_ -- Species' true distributions and taxonomic concepts can change over time, which warrants a
  way to select records from within user-defined time periods.
  * _Occurrence month_ -- Relevant for investigations of migratory species and individual seasons.
  * _Coordinate uncertainty_ -- Geographic coordinates vary immensely among records.  The framework allows the user to exclude records on the basis of coordinate uncertainty.  It is also possible to filter out records that do not have an associated coordinate.
  * _Known issues_ -- Some records have known issues that limit their value.  The framework enables users to exclude
  records based on issues identified for records by providers.
  * _Basis of record and sampling protocols_ -- Datasets accessed through GBIF include a variety of types of records, such as preserved specimens and fossil records.  Additionally, different sampling protocols may have been employed.  The user can choose which types to filter out.
  * _Detection distance_ -- Detection distance is the distance between an observer and an individual animal that is recorded.  It can affect the utility and possible scales of analyses of data because it adds to the locational uncertainty associated with the precision of the Global Positioning System (GPS) reading and observer movement during transects and traveling counts.  Different taxa may be sampled with different methods, and the uncertainty surrounding the exact locations of individuals recorded can vary among methods.  For example, small mammals that are captured in traps can confidently be assigned to the trap's location, but loud-singing birds detected in an auditory survey could be hundreds of meters away from the observer and thus, the coordinate associated with the record.  Some researchers choose not to address this, while others structure analyses around it.  The framework allows for either strategy, and anything in between, and the decision can be documented. A default max detection distance can be specified for the taxa concepts and then over-ridden during query execution.  
  * _Locational uncertainty of recorded individuals_ -- As mentioned in *detection distance* and *coordinate uncertainty*, although records are recorded as x,y coordinates, there are varying degrees of uncertainty regarding the exact locations of individuals recorded ("locational uncertainty").  Here, locational uncertainty is defined as the sum of gps precision and the maximum possible detection distance of the species during the survey or sampling event, as well as the distance traveled if the sampling protocol was a transect.  Within GBIF, this information is rarely attributed to records, so researchers must make assumptions or best guesses.  The framework provides a means of documenting and explaining those choices.  Additionally, points are buffered according to locational uncertainty of each record.  The buffer radii are the reported or assumed coordinate uncertainty plus a user-provided value for maximum detection distance.  Records from transects or traveling counts may need their length added to given or assumed coordinate uncertainty from gps or other source of geolocation.   
  * _Dynamism in occurrence data sets_ -- Some datasets enable data contributors to go back and edit attributes of occurrence records.  In addition, historic records may be added that change the set of records associated with a past time period.  That is to say that a query of years past that was run today may be different than the same query run tomorrow.  This represents a challenge for provenance.  The method to handle this here is to document data request parameters and post-request filter sets as uniquely identifiable objects that are stored and documented in a database ('wildlife-wrangler.sqlite').  Records included in the wildlife-wrangler output are linked to the filter sets used to acquire them.
  * _Taxonomic issues_ -- Failure to account for taxonomic issues, such as species name changes, synonyms, homonyms, and taxon concept changes, can create problems for studies of species' geographic distributions that use species occurrence records.  The potential consequences of these errors include commission errors, inflated omission rates, and missed opportunities for model validation.  The wildlife-wrangler.sqlite database can facilitate efforts to avoid those errors.  This topic is discussed in detail within the user's guide.


## Recent changes (July 1, 2020)
* Added ability to limit requests to within geometries.
* Added ability to specify a limiting polygon for a species.
* Occurrence record database (output) now includes column for weight that users can use to omit or devalue undesirable records.
* Added config file template.  Config file is necessary to avoid sharing your email password when performing large queries of GBIF.
* Improved handling of duplicates (see Duplicate handling sections above for more information).

## Coming soon
* Ability to incorporate bird records directly from a copy of the eBird Basic Dataset ([EBD](https://ebird.org/about/download-ebird-data-products)) that the user has downloaded.  The eBird Observational Dataset ([EOD](https://ebird.org/about/download-ebird-data-products)) is available through GBIF, but that dataset does not include some valuable information that is available in the EBD.
* Making species level geometry filtering optional if polygon is present in taxa concepts table.
* Incorporating GBIF fields "dataGeneralizations", "georeferenceRemarks", and "informationWitheld".
* Overriding polygon geometry columns in output database if a "footprintWKT" value was provided.
* Incorporating tools for navigating taxa concept matching and assessment.

## Inputs
Data is gathered from databases via GBIF, so there are few inputs.  However, the 'wildlife-wrangler.sqlite' database is needed, which includes tables for taxa definitions, data request parameters, and post-request filtering criteria.

## Outputs
On a per-species, per-query basis
* A database of filtered species occurrence records with documentation.  The format supports display in various Geographic Information Systems (GIS).
* Notebook documents that describe decisions made and the data acquired.

## Constraints
* Currently only works with Windows operating systems.  
* Queries returning > 5,000,000 records may fail.
* Setup of spatialite can be difficult.
* Processing speed is limited in some cases by lack of spatial indexing, because setup of spatialite with spatial indexing enabled is very difficult.

## Dependencies
Python 3 and numerous packages including sqlite3 with the spatialite extension are needed.  Running the following code in a conda shell should create a suitable environment named "wrangler":
1. "conda create -n wrangler python=3.6 pandas jupyter basemap-data-hires notebook numpy shapely"
2. "conda activate wrangler"
3. "pip install pygbif python-dwca-reader"

See the included spatialite install notes for how to install spatialite on Windows 10.  

## Code
All code is included in this repository.  Runtimes of discrete tasks made grouping code into separate functions preferable.  

## Suggested Citation
Tarr, N. M., McKerrow, A. J., and M. J. Rubino. 2021. The Wildlife Wrangler. U.S. Geological Survey software release, accessed January 6, 2021, at https://doi.org/10.5066/P98K7E93

## Copyright and License
Unless otherwise noted, This project is in the public domain in the United States because it contains materials that originally came from the United States Geological Survey, an agency of the United States Department of Interior. For more information, see the official USGS copyright policy at https://www.usgs.gov/information-policies-and-instructions/copyrights-and-credits

Additionally, we waive copyright and related rights in the work worldwide through the CC0 1.0 Universal public domain dedication.

This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.
