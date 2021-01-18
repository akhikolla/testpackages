[![](https://cranlogs.r-pkg.org/badges/beam)](https://cran.r-project.org/package=beam)
[![](https://cranlogs.r-pkg.org/badges/grand-total/beam)](https://cran.r-project.org/package=beam)
[![Travis build status](https://travis-ci.org/gleday/beam.svg?branch=master)](https://travis-ci.org/gleday/beam)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/gleday/beam?branch=master&svg=true)](https://ci.appveyor.com/project/gleday/beam)
[![Coverage status](https://codecov.io/gh/gleday/beam/branch/master/graph/badge.svg)](https://codecov.io/github/gleday/beam?branch=master)

# beam

This R package implements the method described in

Leday, G.G.R. and Richardson, S. (2019). [Fast Bayesian inference in large Gaussian graphical models](https://doi.org/10.1111/biom.13064). *Biometrics.* 75(4), 1288--1298.

## Description

Fast Bayesian inference of marginal and conditional independence structures from high-dimensional data.

## Installation

If you wish to install **beam** from R:

```R
# Install/load R package devtools
install.packages("devtools")
library(devtools)

# Install/load R package beam from github
install_github("gleday/beam")
library(beam)
```

Note that **beam** is maintained on github and only updated on CRAN every so often.
Therefore, software versions may differ.
