This file lists changes in the R package for each release.

## beam 2.0.2

* changed contact email address
* seperate centering and scaling c++ functions to compute sample variances

## beam 2.0.1

* Small bugs fixed in c++ code

## beam 2.0.0

### Major changes

* the main code of `beam()` is now written in C++

* new `D` argument for `beam()` to specify the shrinkage target of the linear shrinkage estimator

* new function `lightbeam()` that is a wrapper for `beam()` and `beam.select()`. It is faster and outputs a lighter object (sparse matrix).

* new methods `postExpSigma` and `postExpOmega` for `beam-class` objects that provide the posterior expectation of covariance and inverse-covariance matrices, respectively. These methods also allow the shrinkage estimation of the variances (or not).

### Minor changes

* Dependency on Rcpp has been added

* date and description have been added in ``DESCRIPTION`

* a `NEWS` file has been added

* paper citations have been updated

