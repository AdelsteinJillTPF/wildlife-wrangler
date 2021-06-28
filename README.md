# The Wildlife Wrangler
__Maintainer__

Nathan M. Tarr (nmtarr@ncsu.edu, [orcid:0000-0003-2925-8948](https://orcid.org/0000-0003-2925-8948), github.com/nmtarr)

__Contributors__

Alexa McKerrow (amckerrow@usgs.gov, [orcid:0000-0002-8312-2905](https://orcid.org/0000-0002-8312-2905))

Matthew Rubino (mjrubino@ncsu.edu, [orcid:0000-0003-0651-3053](https://orcid.org/0000-0003-0651-3053))

## USGS Software Release Information
IP-120038
https://doi.org/10.5066/P98K7E93

## Purpose
The abundance of wildlife occurrence data sets that are currently accessible can be valuable for efforts such as species distribution modeling and range delineation.  However, the task of downloading and filtering occurrence records is often complex due to errors and uncertainties that are present in data sets (Tessarolo et al. 2017).  This repository provides a framework for acquiring and filtering occurrence data that are freely available through the Global Biodiversity Information Facility's ([GBIF](https://gbif.org)) API and eBird Basic Dataset ([EBD](https://ebird.org/science/use-ebird-data/)).     

## Framework
Records are requested from occurrence data sets and filtered according to species- and query-specific parameters.  Filtered occurrence records are saved in a database.  The details of taxa concepts and filter parameters can be stored as JSON files for reuse and reference.  Additionally, Jupyter Notebook documents are created that describe the filtered data sets for the sake of documentation and filter set refinement.

## Features
This framework has certain features that support transparency and build confidence for analyses and evaluations of species distributions.

* __Automation__ -- The potential volume of data involved necessitates the processes be automated. Automation also reduces subjectivity in decision making, enables thorough documentation, and ensures repeatability.  However, some aspects of data wrangling and filtering are unavoidably analog.

* __Open source components__ -- All processes are run with Python, but are coded in Python, R, and SQL and use sqlite3, a built-in Python package.

* __Detailed parameterization__ -- Queries and filters can be parameterized on a per-species and per-event basis. Rules do not have to be applied generally across large numbers of species or evaluations.

* __Transparency__ -- Summaries of occurrence data and models using empirical data include subjectivity in the form of choices regarding parameters and rules for curating input data.  This framework is meant to provide a way to document those choices.  

* __Geospatial processing__ -- Some geospatial operations are performed including re-projecting and buffering points.  A helper function is provided that creates a shapefile of buffered points from an output database.

* __Data summarization__ -- The attributes of occurrence records are summarized in a Jupyter Notebook document.

* __Spatial filtering__ -- Queries can be limited to within geometries (polygons), and spatial restrictions can be assigned to taxa circumscriptions (i.e., species' extent of occurrence polygons).  The user can also specify a continent and/or country within which to return records.

* __Source filtering__ -- Global Biodiversity Information Facility ([GBIF](https://gbif.org)) aggregates records from many datasets, collections, and institutions.  The user can specify collections and institutions to omit.  

* __Duplicate handling__ -- Queries commonly include duplicates based on the latitude, longitude, and date fields.  The user can opt to keep or exclude duplicates.  See the User's Guide for more details.

* __Weighting__ -- Output databases include fields for weighting records based on their attributes.  The wrangler does not currently offer the capability to adjust weights, but the user can do so after the output is generated.  Weighting can be used to quantify data degradation, and thus improve models and other applications (Tessarolo et al. 2017).  

* __Wildlife-centric design__ -- The framework addresses several filter criteria that are especially relevant to species occurrence records for studies of wildlife distributions.
  * _Occurrence year_ -- Species' true distributions and taxonomic concepts can change over time, which warrants a way to select records from within user-defined time periods.
  * _Occurrence month_ -- Relevant for investigations of migratory species and individual seasons.
  * _Locational uncertainty_ -- Data providers associate occurrence records with geographic locations through a process called georeferencing.  However, the spatial accuracy, precision, and detail of georeferences can vary immensely among records (Chapman and Wieczorek 2020).  Thus, there are varying degrees of uncertainty regarding the locations of individuals recorded and/or the areas sampled by observers ("locational uncertainty").  Within GBIF and eBird, georeferences are often limited in detail, so researchers must make assumptions or best guesses.  The framework provides a means of documenting and explaining those choices and includes processes that interpret and approximate georeferences. It is also possible to filter on locational uncertainty.  See the user's guide for more information on this topic.
  * _Detection distance_ -- The distance between an observer and an individual animal that is recorded can affect the usability of data and limit the scales of analyses that are possible.  That is because detection distance adds to the locational uncertainty associated with records from other factors such as the precision of the Global Positioning System (GPS) reading and observer movement during transects and traveling counts.  Different taxa may be sampled with different methods, and the uncertainty surrounding the given locations of individuals recorded can vary among methods.  For example, small mammals that are captured in traps can confidently be assigned to the trap's location, but loud-singing birds detected in an auditory survey could be hundreds of meters away from the observer and thus, the coordinate associated with the record.  Some researchers choose not to address this, while others structure analyses around it.  The framework allows for either strategy and anything in between, and the decision can be documented. A default maximum detection distance can be specified for the taxa concepts or during query execution. See section 2.3 in Chapman and Wieczorek (2020) for more on this topic.
  * _Dynamism in occurrence data sets_ -- Some datasets enable data contributors to revise attributes of occurrence records after they are added to data sets.  In addition, historic records may be added that change the set of records associated with a past time period.  That is to say that a query of years past that was run today may be different than the same query run tomorrow.  This represents a challenge for provenance.  The method used to handle the issue is to document taxon concepts and parameters for filtering data as JSON objects that can be saved, documented, and reused.  Additionally, records included in the wildlife-wrangler output are linked to the filter sets used to acquire them.
  * _Taxonomic issues_ -- Failure to account for taxonomic issues, such as species name changes, synonyms, homonyms, and taxon concept changes, can create problems for studies of species' geographic distributions that use species occurrence records (Tessarolo et al. 2017).  The potential consequences of these errors include commission errors, inflated omission rates, and missed opportunities for model validation.  Careful specification of taxon concepts in the wildlife wrangler can help avoid such errors.
  * _Known issues_ -- Some records have known issues that limit their value.  The framework enables users to exclude records based on issues identified for records by providers.
  * _Basis of record and sampling protocols_ -- Data sets accessed through GBIF include a variety of types of records, such as preserved specimens and fossil records.  Additionally, different sampling protocols may have been employed.  The user can choose which types to filter out.


## Recent changes (June 22, 2021)
* Made species level geometry filtering optional if polygon is present in taxon concept dictionary.
* Included citations, rights, download key, and download doi for GBIF downloads.
* Increased the number of attributes retained for GBIF records.
* Ability to incorporate bird records directly from a copy of the eBird Basic Dataset ([EBD](https://ebird.org/about/download-ebird-data-products)) that the user has downloaded.  The eBird Observational Dataset ([EOD](https://ebird.org/about/download-ebird-data-products)) is available through GBIF, but that dataset does not include some valuable information that is available in the EBD.
* Abandonment of Spatialite in favor of the Python Geopandas package for spatial processes.  
* All software requirements can now be satisfied with conda. There is no need to suffer through Spatialite installation.
* Revised Jupyter Notebook document for designing and running queries.
* More sophisticated logic for interpreting and approximating the geographic locations of records in accordance with Chapman and Wieczorek (2020).
* Added option to helper function that creates shapefiles from output that enables user to generate a point, record footprint, or extract a random point from within each occurrence records' footprint.

## Coming soon
* Ability to interpret the shape georeference methods by using the Darwin Core attribute "footprintWKT"
* Incorporating tools for navigating taxa concept matching and assessment.

## Inputs
Data is gathered from databases via GBIF, and records can be queried from a copy of the EBird Basic Dataset that the user has acquired from eBird.  The user also builds or inputs taxon concepts and filter parameters.

## Outputs
On a per-species, per-query basis
* A database of filtered species occurrence records with documentation and summaries.
* Notebook documents that describe decisions made by the user and summarize the properties of acquired records.

## Constraints
* Currently may only work with Windows 10 operating systems.  
* Filtering with geometry is not functional due to pygbif bugs.
* Queries returning > 5,000,000 records may fail.
* Processing speed is limited in some cases by lack of spatial indexing and slow speed of eBird's Auk R package.

## Dependencies
Python 3 and numerous packages including sqlite3 with the spatialite extension are needed.  Running the following code in a conda shell should create a suitable environment named "wrangler":
1. "conda create -n wrangler python=3.6 ipython rpy2 pandas jupyter notebook numpy shapely matplotlib r-base geopandas descartes"
2. "conda activate wrangler"
3. "pip install pygbif python-dwca-reader"

## Code
All code is included in this repository.  Runtimes of discrete tasks made grouping code into separate functions preferable.  

## Suggested Citation
Tarr, N. M., McKerrow, A. J., and M. J. Rubino. 2021. The Wildlife Wrangler. U.S. Geological Survey software release, accessed June 6, 2021, at https://doi.org/10.5066/P98K7E93

## Copyright and License
Unless otherwise noted, This project is in the public domain in the United States because it contains materials that originally came from the United States Geological Survey, an agency of the United States Department of Interior. For more information, see the official USGS copyright policy at https://www.usgs.gov/information-policies-and-instructions/copyrights-and-credits

Additionally, we waive copyright and related rights in the work worldwide through the CC0 1.0 Universal public domain dedication.

This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.

## References
Chapman, A.D. & Wieczorek, J.R. (2020) Georeferencing Best Practices. Copenhagen: GBIF Secretariat. https://doi.org/10.15468/doc-gg7h-s853

Tessarolo, G., Ladle, R., Rangel, T., and Hortal, J. 2017. Temporal degradation of data limits biodiversity research. Ecology and Evolution.  2017;7:6863–6870. DOI: 10.1002/ece3.3259.
