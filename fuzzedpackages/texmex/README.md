The texmex package for R
========================

[![Build Status](https://img.shields.io/travis/harrysouthworth/texmex/master.svg)](https://travis-ci.org/harrysouthworth/texmex)
[![Code Coverage](https://codecov.io/github/harrysouthworth/texmex/branch/master/graphs/badge.svg)](https://codecov.io/github/harrysouthworth/texmex)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/texmex)](https://CRAN.R-project.org/package=texmex)

Extreme value modelling with R. Includes univariate
modelling using generalized Pareto,
generalized extreme value, Weibull and Gumbel
distributions, and the
multivariate conditional approach of Heffernan and
Tawn.

The package contains a test suite that depends on the
testthat package. To use the test suite, install
using 'R CMD INSTALL --install-tests' and then running
'devtools::test("texmex")' within R, where "texmex" points
to the package location.

This work was partially funded by AstraZeneca.

Currently (2019-11-01) for me (Harry) to get roxygen to work, it seems necessary
to do

remotes::install_version("roxygen2", "6.0.1")

I'm sure there are other workarounds and the issue may be system
dependent. There's a whole lovefest about it under Issue 822
in the roxygen2 GitHub.
