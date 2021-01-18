
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis-CI Build
Status](https://travis-ci.org/SMAC-Group/avar.svg?branch=master)](https://travis-ci.org/SMAC-Group/avar)
[![Project Status:
Active](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Licence](https://img.shields.io/badge/licence-AGPL--3.0-blue.svg)](https://opensource.org/licenses/AGPL-3.0)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.4.0-6666ff.svg)](https://cran.r-project.org/)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/avar)](http://www.r-pkg.org/pkg/avar)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/avar)](http://www.r-pkg.org/pkg/avar)
[![Last-changedate](https://img.shields.io/badge/last%20change-2020--01--15-yellowgreen.svg)](https://github.com/SMAC-Group/avar)

# Welcome to the `avar` package <a href="https://smac-group.com/"><img src="man/figures/logo.png" align="right" style="width: 16%; height: 16%"/></a>

This package provides the tools necessary to compute the empirical Allan
Variance (AV) and use it to estimate the parameters of (latent) time
series models. The estimation of the Allan Variance is performed through
the estimator proposed by Allan (1966) and, based on this quantity, the
Allan Variance Linear Regression (AVLR) approach (or Allan Variance
Slope Method) is often used by engineers to retrieve the parameters of
time series models which are assumed to underlie the observed signals
(see for example Guerrier, Molinari, and Stebler 2016). These estimators
are implemented in this package along with the relevant plotting and
summary functions.

## Install Instructions

The `avar` package is available on both CRAN and GitHub. The CRAN
version is considered stable while the GitHub version is subject to
modifications/updates which may lead to installation problems or broken
functions. You can install the stable version of the `avar` package
with:

``` r
install.packages("avar")
```

For users who are interested in having the latest developments, the
[GitHub](https://github.com/SMAC-Group/avar) version is ideal although
more dependencies are required to run a stable version of the package.
Most importantly, users **must** have a (C++) compiler installed on
their machine that is compatible with R (e.g. Clang). Once you’ve made
sure that you have a compatible C++ compiler installed on your computer,
run the following code in an R session and you will be ready to use the
devlopment version of `avar`.

``` r
# Install dependencies
install.packages(c("devtools"))

# Install/Update the package from GitHub
devtools::install_github("SMAC-Group/avar")

# Install the package with Vignettes/User Guides 
devtools::install_github("SMAC-Group/avar", build_vignettes = TRUE)
```

## References

<div id="refs" class="references">

<div id="ref-allan1966statistics">

Allan, David W. 1966. “Statistics of Atomic Frequency Standards.”
*Proceedings of the IEEE* 54 (2). IEEE: 221–30.

</div>

<div id="ref-guerrier2016theoretical">

Guerrier, Stéphane, Roberto Molinari, and Yannick Stebler. 2016.
“Theoretical Limitations of Allan Variance-Based Regression for Time
Series Model Estimation.” *IEEE Signal Processing Letters* 23 (5). IEEE:
597–601.

</div>

</div>
