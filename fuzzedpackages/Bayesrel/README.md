
# Bayesrel

<!-- badges: start -->
[![R build status](https://github.com/juliuspf/Bayesrel/workflows/R-CMD-check/badge.svg)](https://github.com/juliuspf/Bayesrel/actions)
<!-- badges: end -->

The goal of Bayesrel is to provide popular single-test reliability coefficients in the Bayesian statistics framework.
These coefficients include alpha, lambda2, lambda6, the glb, and omega.

## Installation

You can install the released version of Bayesrel from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("Bayesrel")
```
or install the latest version of Bayesrel from [github] (https:://github.com) with the help of the devtools-package:

```r
devtools::install_github("juliuspf/Bayesrel")
```

## Example

This is a basic example which shows you how to compute alpha, lambda2, the glb, and omega for an example real data set:

``` r
library(Bayesrel)
## basic example code
## load example data set from the package
data <- asrm
## run the main reliability function
res <- strel(asrm)
## get a full result output
summary(strel)
## return the probability that coefficient alpha is larger than .70
p_strel(res, estimate = "alpha", low.bound = .70)
```
