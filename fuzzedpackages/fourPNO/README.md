
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.org/tmsalab/fourPNO.svg)](https://travis-ci.org/tmsalab/fourPNO)
[![Package-License](http://img.shields.io/badge/license-GPL%20\(%3E=2\)-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN Version
Badge](http://www.r-pkg.org/badges/version/fourPNO)](https://cran.r-project.org/package=fourPNO)
[![CRAN
Status](https://cranchecks.info/badges/worst/fourPNO)](https://cran.r-project.org/web/checks/check_results_fourPNO.html)
[![RStudio CRAN Mirror’s Monthly
Downloads](http://cranlogs.r-pkg.org/badges/fourPNO?color=brightgreen)](http://www.r-pkg.org/pkg/fourPNO)
[![RStudio CRAN Mirror’s Total
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/fourPNO?color=brightgreen)](http://www.r-pkg.org/pkg/fourPNO)
[![Coverage
status](https://codecov.io/gh/tmsalab/fourPNO/branch/master/graph/badge.svg)](https://codecov.io/github/tmsalab/fourPNO?branch=master)

# `fourPNO` R package

Estimate Barton & Lord’s (1981)
<https://doi.org/10.1002/j.2333-8504.1981.tb01255.x> four parameter IRT
model with lower and upper asymptotes using Bayesian formulation
described by Culpepper (2016)
<https://doi.org/10.1007/s11336-015-9477-6>.

## Installation

You can install `fourPNO` from CRAN using:

``` r
install.packages("fourPNO")
```

Or, you can be on the cutting-edge development version on GitHub using:

``` r
if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("tmsalab/fourPNO")
```

## Usage

To use the `fourPNO` package, load it into *R* using:

``` r
library("fourPNO")
```

Inside the package, the estimation routines can be viewed as:

  - `Gibbs_2PNO()`
  - `Gibbs_4PNO()`

## Author

Steven Andrew Culpepper

## Citing the `fourPNO` package

To ensure future development of the package, please cite `fourPNO`
package if used during an analysis or simulation study. Citation
information for the package may be acquired by using in *R*:

``` r
citation("fourPNO")
```

## License

GPL (\>= 2)
