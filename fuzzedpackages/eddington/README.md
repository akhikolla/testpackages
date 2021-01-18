
<!-- README.md is generated from README.Rmd. Please edit that file -->
eddington
=========

<!-- badges: start -->
[![Travis-CI Build Status](https://travis-ci.org/pegeler/eddington2.svg?branch=master)](https://travis-ci.org/pegeler/eddington2) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/pegeler/eddington2?branch=master&svg=true)](https://ci.appveyor.com/project/pegeler/eddington2) <!-- badges: end -->

One statistic that cyclists are known to track is their Eddington number. The Eddington number for cycling, *E*, is the maximum number where a cyclist has ridden *E* miles in *E* days. So to get a number of 30, you need to have ridden 30 miles or more on 30 separate days.

This package provides functions to compute a cyclist's Eddington number, including efficiently computing cumulative *E* over a vector. These functions may also be used for *h*-indices for authors. Both are specific applications of computing the side length of a Durfee square.

Installation
------------

<!--
You can install the released version of eddington from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("eddington")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pegeler/eddington2")
```
-->
``` r
# install.packages("devtools")
devtools::install_github("pegeler/eddington2/R/package")
```

Example
-------

Here is a basic example showing how to get the summary Eddington number of the included **simulated** `rides` dataset. Note that we first have to aggregate ride mileage by date.

``` r
library(eddington)
library(dplyr)

rides %>%
  group_by(ride_date) %>%
  summarize(n = n(), total = sum(ride_length)) %>% 
  summarize(E = E_num(total)) %>% 
  pull
#> [1] 29
```

See the package vignette for detailed usage.
