
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/kota7/combiter.svg?branch=master)](https://travis-ci.org/kota7/combiter) [![CRAN Status](https://www.r-pkg.org/badges/version/combiter)](https://www.r-pkg.org/badges/version/combiter) [![](https://cranlogs.r-pkg.org/badges/combiter)](https://cran.r-project.org/package=combiter) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/kota7/combiter?branch=master&svg=true)](https://ci.appveyor.com/project/kota7/combiter)

combiter: Combinatorics Iterators for R
=======================================

This package provides iterators for combinations, permutations, subsets, and cartesian product, with which one can go through the elements without creating a huge set of all possible values.

Since `v1.0.2`, the iterator objects of the package inherit the `iter` class from [iterators](https://cran.r-project.org/package=iterators) package. As a result, the objects are also compatible with [foreach](https://cran.r-project.org/package=foreach) library.

Installation and import
-----------------------

Install from CRAN by:

``` r
install.packages("combiter")
```

Or install development version by:

``` r
devtools::install_github("kota7/combiter")
```

Use by

``` r
library(combiter)
```

Usage
-----

Please visit the [Documentation Page](https://kota7.github.io/combiter) for the package introduction.
