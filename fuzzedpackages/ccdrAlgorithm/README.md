ccdrAlgorithm
=============

[![Project Status: Active The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![Travis-CI Build Status](https://travis-ci.org/itsrainingdata/ccdrAlgorithm.svg?branch=master)](https://travis-ci.org/itsrainingdata/ccdrAlgorithm) [![](http://www.r-pkg.org/badges/version/ccdrAlgorithm)](http://www.r-pkg.org/pkg/ccdrAlgorithm) [![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/ccdrAlgorithm)](http://www.r-pkg.org/pkg/ccdrAlgorithm)

`ccdrAlgorithm` implements the CCDr structure learning algorithm described in \[[1-2](#references)\]. This algorithm estimates the structure of a Bayesian network from mixed observational and experimental data using penalized maximum likelihood based on L1 or concave (MCP) regularization.

Presently, this package implements the main algorithm and provides a method to simulate data from a Gaussian Bayesian network. To simulate random networks, it is recommended to use the [`sparsebnUtils`](https://cran.r-project.org/package=sparsebnUtils) package. Other packages for simulating DAGs and observational data include [`bnlearn`](https://cran.r-project.org/package=bnlearn), [`pcalg`](https://cran.r-project.org/package=pcalg), and [`igraph`](https://cran.r-project.org/package=igraph).

Overview
--------

The main method is `ccdr.run`, which runs the CCDr structure learning algorithm as described in \[[1-2](#references)\]. For simulating data from a Gaussian Bayesian network, the package provides the method `generate_mvn_data`. This method can simulate observational data or experimental data with interventions (or combinations of both).

Installation
------------

You can install:

-   the latest CRAN version with

    ``` r
    install.packages("ccdrAlgorithm")
    ```

-   the latest development version from GitHub with

    ``` r
    devtools::install_github(c("itsrainingdata/sparsebnUtils/dev", "itsrainingdata/ccdrAlgorithm/dev"))
    ```

References
----------

\[1\] Aragam, B. and Zhou, Q. (2015). [Concave penalized estimation of sparse Gaussian Bayesian networks.](http://jmlr.org/papers/v16/aragam15a.html) *The Journal of Machine Learning Research*. 16(Nov):2273−2328.

\[2\] Zhang, D. (2016). Concave Penalized Estimation of Causal Gaussian Networks with Intervention. Master’s thesis, UCLA.

\[3\] Fu, F. and Zhou, Q. (2013). [Learning sparse causal Gaussian networks with experimental intervention: Regularization and coordinate descent.](http://amstat.tandfonline.com/doi/abs/10.1080/01621459.2012.754359) Journal of the American Statistical Association, 108: 288-300.
