
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gsisdecoder

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/gsisdecoder)](https://CRAN.R-project.org/package=gsisdecoder)
[![R build
status](https://github.com/mrcaseb/gsisdecoder/workflows/R-CMD-check/badge.svg)](https://github.com/mrcaseb/gsisdecoder/actions)
[![Travis build
status](https://travis-ci.com/mrcaseb/gsisdecoder.svg?branch=master)](https://travis-ci.com/mrcaseb/gsisdecoder)
[![Codecov test
coverage](https://codecov.io/gh/mrcaseb/gsisdecoder/branch/master/graph/badge.svg)](https://codecov.io/gh/mrcaseb/gsisdecoder?branch=master)
<!-- badges: end -->

The goal of gsisdecoder is to provide a function wrapper for the high
efficient GSIS ID decoder written in c++. The function will be used
within the package
[nflfastR](https://mrcaseb.github.io/nflfastR/index.html) but since it
needs to be precompiled it was outsourced from the mentioned package
into this standalone package.

## Installation

You can install the released version of gsisdecoder from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("gsisdecoder")
```
