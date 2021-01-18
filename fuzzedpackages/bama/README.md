<!-- badges: start -->
[![CRAN
Version](https://img.shields.io/cran/v/bama?style=flat-square&color=blue&label=CRAN)](https://cran.r-project.org/package=bama)
[![GitHub
Release](https://img.shields.io/github/v/release/umich-cphds/bama?include_prereleases&label=Github&style=flat-square)](https://github.com/umich-cphds/bama)
[![Travis
CI](https://img.shields.io/travis/umich-cphds/bama?style=flat-square)](https://travis-ci.org/umich-cphds/bama)

Bayesian Mediation Analysis
===========================

Perform mediation analysis in the presence of high-dimensional mediators
based on the potential outcome framework. Bayesian Mediation Analysis
(BAMA), developed by Song et al (2019), relies on two Bayesian sparse
linear mixed models to simultaneously analyze a relatively large number
of mediators for a continuous exposure and outcome assuming a small
number of mediators are truly active. This sparsity assumption also
allows the extension of univariate mediator analysis by casting the
identification of active mediators as a variable selection problem and
applying Bayesian methods with continuous shrinkage priors on the
effects.

Installation
------------

You can install `bama` via CRAN

    install.packages("bama")

Or devtools

    devtools::install_github("umich-cphds/bama", build_opts = c())

The Github version may contain new features or bug fixes not yet present
on CRAN, so if you are experiencing issues, you may want to try the
Github version of the package.

If you wish to install the package via `devtools`, you will need a C++
compiler installed. This can be accomplished by installing Rtools on
Windows and Xcode on MacOS.

Example
-------

This example is taken from the `bama` help file to help you get started
using the method. Please check the documentation of the function by
typing `?bama::bama`, and the vignette by typing `vingette("bama")` in
R.

`bama` includes an example dataset, `bama.data`. It is a `data.frame`
with a numeric response `y`, numeric exposure `a` and 100 numeric
mediators named `m1, m2, ..., m100`.

We recommend using much larger numbers for `burnin` and `ndraws`, for
example (30000, 1000).

    library(bama)

    Y <- bama.data$y
    A <- bama.data$a

    # grab the mediators from the example data.frame
    M <- as.matrix(bama.data[, paste0("m", 1:100)], nrow(bama.data))

    # We just include the intercept term in this example.
    C1 <- matrix(1, nrow(M), 1)
    C2 <- matrix(1, nrow(M), 1)
    beta.m  <- rep(0, 100)
    alpha.a <- rep(0, 100)

    set.seed(1245)
    bama.out <- bama(Y, A, M, beta.m, alpha.a, C1 = C2, C2 = C2, burnin = 1000,
                     ndraws = 100)

    # Rank mediators and see summary information
    head(summary(bama.out, rank = T))
    #>        estimate    ci.lower     ci.upper  pip
    #> m12  0.18576755  0.12218219  0.253780289 0.99
    #> m65 -0.23979113 -0.30642142 -0.165385815 0.99
    #> m89 -0.14252621 -0.21202100 -0.065313879 0.91
    #> m97 -0.04289325 -0.10325313  0.007526357 0.31
    #> m90 -0.03041407 -0.08301040  0.012826983 0.17
    #> m93  0.04051381 -0.00547444  0.118579967 0.17

Reference
=========

Song, Y, Zhou, X, Zhang, M, et al. Bayesian shrinkage estimation of high
dimensional causal mediation effects in omics studies. Biometrics. 2019;
1-11. [doi:10.1111/biom.13189](https://doi.org/10.1111/biom.13189)
