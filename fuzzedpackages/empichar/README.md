
<!-- README.md is generated from README.Rmd. Please edit that file -->
empichar
========

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/gbasulto/empichar.svg?branch=master)](https://travis-ci.org/gbasulto/empichar) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/gbasulto/empichar?branch=master&svg=true)](https://ci.appveyor.com/project/gbasulto/empichar) [![Codecov test coverage](https://codecov.io/gh/gbasulto/empichar/branch/master/graph/badge.svg)](https://codecov.io/gh/gbasulto/empichar?branch=master) [![CRAN status](https://www.r-pkg.org/badges/version/empichar)](https://cran.r-project.org/package=empichar) <!-- badges: end --> <!-- -->

The goal of empichar is to evaluate the empirical characteristic function meeting the following criteria:

1.  Allow the `C++` functions to be imported in other `R` packages.
2.  Make a fast evaluation (using `RcppArmadillo`).
3.  Check the dimensions of the input arguments.
4.  If desired, compute only the real part, imaginary part or modulus of the empirical characteristic function (faster than taking the real part of the empirical characteristic function).

Installation
------------

You can install the released version of empichar from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("empichar")
```

Available functions
-------------------

| Function   | Description                                             |
|:-----------|:--------------------------------------------------------|
| `ecf`      | Empirical characteristic function of a given sample     |
| `ecf_real` | Real part of the empirical characteristic function      |
| `ecf_imag` | Imaginary part of the empirical characteristic function |
| `ecf_mod`  | Modulus of the empirical characteristic function        |

Documentation
-------------

I will add a vignette later on. Meanwhile, consult the documentation with `help("ecf", "empichar")`, `help("ecf_real", "empichar")`, `help("ecf_imag", "empichar")`, `help("ecf_mod", "empichar")`.
