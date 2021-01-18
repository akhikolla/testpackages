LocalControl 1.1.2.2
================

### Minor changes
* Added citation file and doi for vignette's JSS publication.

LocalControl 1.1.2
================

### Minor changes
* Fixed a bug causing survival confidence intervals to be incorrect.
* Updated vignette to match JSS requests.
* Changed cardSim data to match the new simulation in the vignette.

LocalControl 1.1.1
================

### Major changes
* Added vignette

### Minor changes
* Added default value (NULL) to modelForm argument (LocalControl())

LocalControl 1.1.0
================

### Major changes
* Deprecated the localControlCompetingRisks(), localControlNearestNeighbors(), and (hidden) doLocalControl() functions.
* Added a new function in their place, LocalControl().
*   The behavior of the three removed functions now exists in this function. 
*   The new “outcomeType” parameter allows users to toggle between the competingRisks and NearestNeighbors functionality.
* Deprecated the plotLocalControlCIF and plotLocalControlLTD functions, replaced with s3 (plot()) functions. 
* Added S3 functions: print() and summary() have also been added for the LocalControl classes.
* Changed the structure of the LocalControlCS and LocalControlCR objects.
* Removed the summary object. 
*   The summary is now created upon request using the summary() s3 function with the LocalControl classes. 
* Added a formula interface for LocalControl (beta).

     LocalControl(data = lindner, 
                  modelForm = formula('cardbill ~ abcix | stent + female + acutemi'))
     
     Is now a valid alternative to:
	
     LocalControl( data = lindner,
                   clusterVars = c("stent",  "female", "acutemi"),
                   treatmentColName = "abcix",
                   outcomeColName = "cardbill")




LocalControl 1.0.1
================

### Minor changes
* Changed required R version from R (>= 2.10) to R (>= 3.0.0), due to Rcpp dependency.
* Updated title field in DESCRIPTION file to remove redundancies.
* Added Apache Version 2.0 LICENSE file.

### Bug fixes
* Fixed warning: variable length arrays are a C99 feature.
