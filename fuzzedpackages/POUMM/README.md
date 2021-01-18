
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis build
status](https://travis-ci.com/venelin/POUMM.svg?branch=master)](https://travis-ci.com/venelin/POUMM)
[![Coverage
status](https://codecov.io/gh/venelin/POUMM/branch/master/graph/badge.svg)](https://codecov.io/github/venelin/POUMM?branch=master)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/POUMM?color=blue)](https://cran.r-project.org/package=POUMM)
[![Downloads](http://cranlogs.r-pkg.org/badges/POUMM?color=blue)](https://cran.r-project.org/package=POUMM)
[![DOI](https://zenodo.org/badge/115860927.svg)](https://zenodo.org/badge/latestdoi/115860927)

# The Phylogenetic Ornstein-Uhlenbeck Mixed Model

The Phylogenetic Ornstein-Uhlenbeck Mixed Model (POUMM) allows to
estimate the phylogenetic heritability of a continuous trait, to test
hypotheses of neutral evolution versus stabilizing selection, to
quantify the strength of stabilizing selection, to estimate measurement
error and to make predictions about the evolution of a phenotype and
phenotypic variation in a population. The POUMM package provides an easy
and efficient way to perform this variety of analyses on large
macro-evolutionary or epidemic trees. It implements a fast-likelihood
calculation algorithm enabling MCMC-sampling with millions of iterations
within minutes on contemporary multiple core processors. The package
provides functions for configuring the fit of the model and a number of
standard generic functions such as logLik, plot, summary, allowing a
visual and a statistical assessment of the goodness of fit. This is an
important step before using the model fit to answer relevant biological
questions.

# Using the R-package

Here is a quick example on how to use the package on a simulated tree
and trait data:

``` r
# number of tips
N <- 500 

# phylogeny
tr <- ape::rtree(N)

# for the example, simulate trait values on the tree according to a POUMM model.
z <- rVNodesGivenTreePOUMM(
  tree = tr,   
  z0 = 0,      # fixed value at the root
  alpha = 2,   # selection strength of the OU process
  theta = 3,   # long term mean of the OU process
  sigma = 1,   # unit-time standard deviation of the OU process
  sigmae = 1   # standard deviation of the non-heritable component
)[1:N]         # only the values at the N tips will be available in reality

# A combined ML and MCMC fit of the model with default parameter settings.
fit <- POUMM(z, tr)

plot(fit)
summary(fit)
AIC(fit)
BIC(fit)
coef(fit)
logLik(fit)
fitted(fit)
plot(resid(fit))
abline(h=0)

# fit PMM to the same data and do a likelihood ratio test
fitPMM <- POUMM(z, tr)
lmtest::lrtest(fitPMM, fit)
```

For an introduction to the model parameters and the package, read the
[User guide](https://venelin.github.io/POUMM/articles/UserGuide.html).
More advanced topics, such as parametrizations and interpretations of
the model fit are covered in the other package
[vignettes](https://venelin.github.io/POUMM/articles/index.html) and in
the package help-pages, e.g. `?POUMM`, `?specifyPOUMM`,
`?summary.POUMM`, `?plot.POUMM`.

# Installing the R-package

Read the section [Installing the POUMM
R-package](https://venelin.github.io/POUMM/articles/UserGuide.html#Installing)
in [Get
started](https://venelin.github.io/POUMM/articles/UserGuide.html).

# Package source-code

The package source-code is available on
[github](https://github.com/venelin/POUMM).

# Licence and copyright

Copyright 2015-2020 Venelin Mitov

Source code to POUMM is made available under the terms of the GNU
General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.
POUMM is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
details.

# Package web-page

Check-out the package
[web-page](https://venelin.github.io/POUMM/index.html) for the latest
news and further documentation.
