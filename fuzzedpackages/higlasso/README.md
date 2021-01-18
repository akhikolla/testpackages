<!-- Badges -->
[![CRAN
Version](https://img.shields.io/cran/v/higlaso?style=flat-square&color=blue&label=CRAN)](https://cran.r-project.org/package=higlasso)
[![GitHub
Release](https://img.shields.io/github/v/release/umich-cphds/higlasso?include_prereleases&label=Github&style=flat-square&color=blue)](https://github.com/umich-cphds/higlasso)
[![Travis
CI](https://img.shields.io/travis/umich-cphds/higlasso?style=flat-square)](https://travis-ci.org/umich-cphds/higlasso)

Hierarchical Integrative Group LASSO
====================================

Environmental health studies are increasingly measuring multiple
pollutants to characterize the joint health effects attributable to
exposure mixtures. However, the underlying dose-response relationship
between toxicants and health outcomes of interest may be highly
nonlinear, with possible nonlinear interaction effects. Hierarchical
integrative group least absolute shrinkage and selection operator
(HiGLASSO) is a general framework to identify noteworthy nonlinear main
and interaction effects in the presence of group structures among a set
of exposures.

Installation
------------

`higlasso` can be installed via Github using `devtools`

    # install.packages("devtools")
    devtools::install_github("umich-cphds/higlasso")

You'll need a working C++11 compiler, which can obtained by installing
Xcode on MacOS, and RTools on Windows.

Example
-------

`higlasso` can be slow, so it may may be beneficial to tweak some its
settings (for example, `nlambda1` and `nlambda2`) to get a handle on how
long the method will take before running the full model.

    library(higlasso)

    set.seed(48109)

    X <- as.matrix(higlasso.df[, paste0("V", 1:10)])
    Y <- higlasso.df$Y
    Z <- matrix(1, nrow(X))

    # This can take a bit of time
    fit <- cv.higlasso(Y, X, Z)

    print(fit)

    ## 'cv.higlaso' fit:
    ## Average cross validation error for each (lambda1, lambda2)
    ##           l2.1     l2.2     l2.3     l2.4     l2.5     l2.6     l2.7     l2.8
    ## l1.1  31.27646 30.16180 26.95437 27.53877 29.02558 30.07576 26.47052 31.89449
    ## l1.2  29.41992 29.10689 27.10846 26.87978 27.47474 27.63711 26.47278 26.85553
    ## l1.3  27.54288 29.05916 26.93971 27.08178 28.45151 29.09442 24.79442 21.60319
    ## l1.4  26.54363 28.53124 26.53870 25.96099 26.54686 28.31096 30.42201 28.16428
    ## l1.5  26.81162 27.25271 24.09822 25.15356 22.58817 23.14864 27.22033 30.09574
    ## l1.6  27.16120 27.58498 23.66292 22.17489 22.00758 24.25419 26.21524 29.39816
    ## l1.7  27.22267 27.51427 23.86357 22.68843 22.47261 25.06126 29.36773 29.38315
    ## l1.8  27.38411 27.67919 24.78215 22.76672 21.47550 24.58022 31.35179 29.82478
    ## l1.9  27.50484 27.74833 24.55036 23.56705 22.28636 24.74533 31.95273 36.06778
    ## l1.10 27.55453 27.81106 25.47948 24.48003 22.49150 23.45593 30.46539 30.89732
    ##           l2.9    l2.10
    ## l1.1  30.40880 29.62060
    ## l1.2  26.55717 30.37125
    ## l1.3  25.07946 30.30343
    ## l1.4  27.52339 27.48728
    ## l1.5  30.18009 30.65290
    ## l1.6  28.36495 28.98296
    ## l1.7  29.70256 26.97493
    ## l1.8  31.19091 36.82322
    ## l1.9  28.48208 25.13074
    ## l1.10 29.82692 33.11209
    ## Lambda min:
    ## 0.292944 0.9037191 
    ## Lambda 1 SE:
    ## 0.7951724 0.9037191

Bugs
----

If you encounter a bug, please open an issue on the
[Issues](https://github.com/umich-cphds/higlasso/issues) tab on Github
or send us an email.

Contact
-------

For questions or feedback, please email Jonathan Boss at
<bossjona@umich.edu> or Alexander Rix <alexrix@umich.edu>.

References
----------

A Hierarchical Integrative Group LASSO (HiGLASSO) Framework for
Analyzing Environmental Mixtures. Jonathan Boss, Alexander Rix, Yin-Hsiu
Chen, Naveen N. Narisetty, Zhenke Wu, Kelly K. Ferguson, Thomas F.
McElrath, John D. Meeker, Bhramar Mukherjee. 2020. arXiv:2003.12844
