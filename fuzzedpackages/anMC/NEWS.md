# anMC 0.2.2

## Changes with respect to anMC 0.2.1

* Fixed an bug in `ProbaMax`, `ProbaMin` that previously resulted in an error in the Cholesky decomposition in rare instances.

* updated DESCRIPTION file with URL and author.

* Fixed bug in `selectQdims` that could cause crashes for some user-defined samplers.


## Changes with respect to anMC 0.2.0

* Fixed an issue in `trmvrnorm_rej_cpp` that could prevent a call to the function `conservativeEstimate` from inside another package. 

* updated references.

## Changes with respect to anMC 0.1.0 

* The package requires a C++11 compiler.

* The time measurements now are taken with the `chrono` library from C++11.

* Dropped dependency on `microbenchmark` package.

## Major changes with respect to ConservativeEstimates 0.2.0 

* The function `selectEq` is now called `selectActiveDims`

* The variable `Thresh` is now called `threshold`

* The functions `ProbaMax` and `ProbaMin` now have a user option to choose the function `trmvrnorm` to generate truncated multivariate normal samples. The default option is `trmvrnorm_rej_cpp`. 
