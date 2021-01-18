RWDataPlyr
=================

<!-- badges: start -->
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/RWDataPlyr)](https://cran.r-project.org/package=RWDataPlyr)

[![R build status](https://github.com/BoulderCodeHub/RWDataPlyr/workflows/R-CMD-check/badge.svg)](https://github.com/BoulderCodeHub/RWDataPlyr/actions)

[![codecov](https://codecov.io/gh/BoulderCodeHub/RWDataPlyr/branch/master/graphs/badge.svg)](https://codecov.io/gh/BoulderCodeHub/RWDataPlyr)
<!-- badges: end -->

## Overview

RWDataPlyr is a tool to read and manipulate data generated from [RiverWare<sup>TM</sup>](http://www.riverware.org) simulations in rdf, csv, and nc formats and work with those data in a dplyr pipeline. It provides functions to gather,  aggregate, and summarize data from multiple RiverWare simulations, i.e., scenarios.

## Installation

RWDataPlyr can be installed from CRAN:

```{r, eval = FALSE}
install.packages("RWDataPlyr")
```

Or the development version can be installed from GitHub:

```{r, eval=FALSE}
# install.packages("devtools")
devtools::install_github("BoulderCodeHub/RWDataPlyr")
```

## Usage

RWDataPlyr provides at least three workflows for reading and using RiverWare data:

1. Reading and manipulating a single scenario
    * Fast
    * Best for inspecting a single slot
    * If comparing scenarios, must manually repeat for each scenario
    * Relies on `read_rdf()` and `read_rw_csv()`
2. Summarizing multiple slots of data from a single scenario
    * Repeatable; allows user to process many slots at once
    * Best for producing "polished" analyses of a single scenario
    * Relies on `rdf_aggregate()` and user specified `rwd_agg` object
3. Aggregating and summarizing many scenarios
    * Repeatable; allows user to process many slots for many scenarios at once
    * Repeats summary of a single scenario on multiple scenarios and combines results together
    * Relies on `rw_scen_aggregate()` and user specified `rwd_agg` object

Check out the workflow vignette for more details:

```{r, eval = FALSE}
vignette("rwdataplyr-workflow", package = "RWDataPlyr")
```

## Log
* 2020-04-17: version 0.6.4 available
* 2020-03-03: version 0.6.3 available
* 2018-08-16: version 0.6.2 available (first version available on CRAN)
* 2018-06-07: version 0.6.1 available
* 2018-04-10: version 0.6.0 available
* 2017-05-26: version 0.5.0 available
* 2016-11-01: version 0.4.1.1 available. The package is now actually called RWDataPlyr.
* 2016-10-20: version 0.4.1 available
* Previous versions were originally available as the `RWDataPlot` package
  * 2016-07-13: version 0.4 available
  * 2016-03-22: version 0.3 available
  * 2015-07-01: version 0.2 available
  * 2014-09-16: working to create an R Package from existing code.
  
## Disclaimer

This software is in the public domain because it contains materials that originally came from the U.S. Bureau of Reclamation, an agency of the United States Department of Interior. 

Although this code has been used by Reclamation, no warranty, expressed or implied, is made by Reclamation or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by Reclamation in connection therewith.

This software is provided "AS IS."
