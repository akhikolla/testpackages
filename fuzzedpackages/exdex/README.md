
<!-- README.md is generated from README.Rmd. Please edit that file -->
exdex
=====

[![Travis-CI Build Status](https://travis-ci.org/paulnorthrop/exdex.svg?branch=master)](https://travis-ci.org/paulnorthrop/exdex) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/paulnorthrop/exdex?branch=master&svg=true)](https://ci.appveyor.com/project/paulnorthrop/exdex) [![Coverage Status](https://codecov.io/github/paulnorthrop/exdex/coverage.svg?branch=master)](https://codecov.io/github/paulnorthrop/exdex?branch=master) [![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/exdex)](https://cran.r-project.org/package=exdex)

Estimation of the Extremal Index
--------------------------------

### What does exdex do?

The extremal index *θ* is a measure of the degree of local dependence in the extremes of a stationary process. The `exdex` package performs frequentist inference about *θ* using two types of methodology.

One type ([Northrop, 2015](https://doi.org/10.1007/s10687-015-0221-5)) is based on a model that relates the distribution of block maxima to the marginal distribution of the data, leading to a semiparametric maxima estimator. Two versions of this type of estimator are provided, following [Northrop, 2015](https://doi.org/10.1007/s10687-015-0221-5) and [Berghaus and Bücher, 2018](https://doi.org/10.1214/17-AOS1621)). A slightly modified version of the latter is also provided. Estimates are produced using both disjoint and sliding block maxima, that latter providing greater precision of estimation.

The other type of methodology uses a model for the distribution of threshold inter-exceedance times ([Ferro and Segers, 2003](https://doi.org/10.1111/1467-9868.00401)). Two versions of this type of approach are provided, following [Süveges (2007)](https://doi.org/10.1007/s10687-007-0034-2) and [Süveges and Davison (2010)](https://doi.org/10.1214/09-AOAS292).

### A simple example

The following code estimates the extremal index using the semiparametric maxima estimators, for an example dataset containing a time series of sea surges measured at Newlyn, Cornwall, UK over the period 1971-1976.

``` r
library(exdex)
theta <- spm(newlyn, 20)
theta
#> 
#> Call:
#> spm(data = newlyn, b = 20)
#> 
#> Estimates of the extremal index theta:
#>           N2015   BB2018  BB2018b
#> sliding   0.2392  0.3078  0.2578 
#> disjoint  0.2350  0.3042  0.2542
summary(theta)
#> 
#> Call:
#> spm(data = newlyn, b = 20)
#> 
#>                   Estimate Std. Error Bias adj.
#> N2015, sliding      0.2392    0.01990  0.003317
#> BB2018, sliding     0.3078    0.01642  0.003026
#> BB2018b, sliding    0.2578    0.01642  0.053030
#> N2015, disjoint     0.2350    0.02222  0.003726
#> BB2018, disjoint    0.3042    0.02101  0.003571
#> BB2018b, disjoint   0.2542    0.02101  0.053570
```

### Installation

To get the current released version from CRAN:

``` r
install.packages("exdex")
```

### Vignette

See `vignette("exdex-vignette", package = "exdex")` for an overview of the package.
