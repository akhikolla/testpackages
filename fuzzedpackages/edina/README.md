
<!-- README.md is generated from README.Rmd. Please edit that file -->

# edina

<!-- badges: start -->

[![R build
status](https://github.com/tmsalab/edina/workflows/R-CMD-check/badge.svg)](https://github.com/tmsalab/edina/actions)
[![Package-License](http://img.shields.io/badge/license-GPL%20\(%3E=2\)-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN Version
Badge](http://www.r-pkg.org/badges/version/edina)](https://cran.r-project.org/package=edina)
[![CRAN
Status](https://cranchecks.info/badges/worst/edina)](https://cran.r-project.org/web/checks/check_results_edina.html)
[![RStudio CRAN Mirror’s Monthly
Downloads](http://cranlogs.r-pkg.org/badges/edina?color=brightgreen)](http://www.r-pkg.org/pkg/edina)
[![RStudio CRAN Mirror’s Total
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/edina?color=brightgreen)](http://www.r-pkg.org/pkg/edina)
<!-- badges: end -->

Perform a Bayesian estimation of the Exploratory Deterministic Input,
Noisy “And” Gate (EDINA) cognitive diagnostic model described by Chen et
al. (2018).

## Installation

You can install `edina` from CRAN using:

``` r
install.packages("edina")
```

Or, you can be on the cutting-edge development version on GitHub using:

``` r
if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("tmsalab/edina")
```

## Usage

To use the `edina` package, load it into *R* using:

``` r
library("edina")
```

From there, the EDINA model can be estimated using:

``` r
edina_model = edina(<data>, chain_length = 10000)
```

To compute a model underneath different *K* attribute configured *Q*
matrices, use:

``` r
edina_model = auto_edina(<data>, k = 2:4, chain_length = 10000)
```

**Note:** Higher *K* configured *Q* matrices take longer to estimate.

## Authors

James Joseph Balamuta, Steven Andrew Culpepper, and Jeffrey A. Douglas

## Citing the `edina` package

To ensure future development of the package, please cite `edina` package
if used during an analysis or simulation studies. Citation information
for the package may be acquired by using in *R*:

``` r
citation("edina")
```

## License

GPL (\>= 2)
