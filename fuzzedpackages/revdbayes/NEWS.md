# revdbayes 1.3.9

## Bug fixes and minor improvements

* Tests in `test-gp.R`, `test-gev.R` and `test-bingp.R` have been modified to avoid errors in the upcoming new release of the `testthat` package.

# revdbayes 1.3.8

## Bug fixes and minor improvements

* The functions `grimshaw_gp_mle()`, `gp_pwm()` and `gp_lrs()` are now exported, so that the rust package can access them using :: not :::.

* The hyperlinks to the Grimshaw (1993) paper in the documentation to `grimshaw_gp_mle()` and `set_prior()` have been corrected.

# revdbayes 1.3.7

## Bug fixes and minor improvements

* Fixed a bug in `dgp()` that produced an incorrect value for the log-density (`log = TRUE`) when `shape` is negative and very close to zero and `x = -1/shape`.

# revdbayes 1.3.6

## Bug fixes and minor improvements

* Use `inherits()` to check the class of objects returned from `try()`, rather than `class()`.

* pkgdown documentation at [https://paulnorthrop.github.io/revdbayes/](https://paulnorthrop.github.io/revdbayes/)

# revdbayes 1.3.5

## Bug fixes and minor improvements

* The d/p/q function for the GEV and GP distributions now handle correctly cases where the input has length 0 and/or is `NA` and inputs `Inf` and `-Inf`. 

# revdbayes 1.3.4

## New features

* In `set_bin_prior()` the user can specify their own prior for the binomial probability, by providing an R function.

## Bug fixes and minor improvements

* In `rpost()` and `rpost_rcpp()` an error is thrown if the prior and the model are not compatible.  Previously a warning was given.

* The penultimate example in the documentation for `set_prior()` has been corrected by adding `model = "gp".  The default `model = "gev"` is not appropriate here because the prior is set up for the GP model.

* (This is an amendment to the third minor improvement in the NEWS for v1.3.3.) In `rpost()` and `rpost_rcpp()` an error is thrown if the input threshold `thresh` is lower than the smallest observation in `data`.  This is only checked when `model = "bingp"` or `model = "pp"`.  This not checked when `model = "gp"` because the user may legitimately supply only threshold excesses. (Many thanks to Leo Belzile for spotting this.)

# revdbayes 1.3.3

## Bug fixes and minor improvements

* LF line endings used in inst/include/revdbayes.h and inst/include/revdbayes_RcppExports.h to avoid CRAN NOTE.

* The format of the `data` supplied to `rpost()` and `rpost_rcpp()` is checked and an error is thrown if it is not appropriate.

* In `rpost()` and `rpost_rcpp()` an error is thrown if the input threshold `thresh` is lower than the smallest observation in `data`.  This is only relevant when `model = "gp"`, `model = "bingp"` or `model = "pp"`.

* The summary method for class "evpost" is now set up according to Section 8.1 of the R FAQ at (https://cran.r-project.org/doc/FAQ/R-FAQ.html).

* A bug in `grimshaw_gp_mle` has been fixed, so that now solutions with K greater than 1 are discarded.  (Many thanks to Leo Belzile.)

* In `grimshaw_gp_mle` using the starting value equal to the upper bound can result in early termination of the Newton-Raphson search.  A starting value away from the upper bound is now used (lines 282 and 519 of frequentist.R). (Many thanks to Jeremy Rohmer for sending me a dataset that triggered this problem.)

* In `set_prior()` if `prior = "norm"` or `prior = "loglognorm"` then an explicit error is thrown if `cov` is not supplied. (Many thanks to Leo Belzile.)

* The mathematics in the reference manual has been tidied.

# revdbayes 1.3.2

## Bug fixes and minor improvements

* The arguments to `d/p/q/rgev` and `d/p/q/rgp` now obey the usual conventions for R's dpqr probability distribution functions.

* In `pp_check.evpost` the argument `subtype` is now documented properly.

* The `conf` argument to `kgaps_mle` didn't work properly: `conf = 95` was always used.  This has been corrected.

# revdbayes 1.3.1

## New features

* Bayesian and maximum likelihood inference for the K-gaps model for inferring the extremal index using threshold inter-exceedances times. [Suveges, M. and Davison, A. C. (2010), Model misspecification in peaks over threshold analysis, The Annals of Applied Statistics, 4(1), 203-221. doi:10.1214/09-AOAS292.]

* New vignette: "Inference for the extremal index using the K-gaps model".

## Bug fixes and minor improvements

* Added the attribute `attr(gom, "npy")` (with value 3) to the `gom` dataset.  This is for compatability with the **threshr** package.

* Give an explicit error message if `plot.evpost` is called with the logically incompatible arguments `add_pu = TRUE` and `pu_only = TRUE`.

* The documentation for `set_bin_prior` has been corrected: only in-built priors are available, i.e. it is not possible for the user to supply their own prior.

# revdbayes 1.2.1

## Bug fixes and minor improvements

* In some extreme cases (datasets with very small numbers of threshold excesses) calling `predict.evpost` with `type = "q"` and `x` close to 1 returns an imprecise value for the requested predictive quantiles.  This has been corrected by using `stats::uniroot` rather than `stats::nlminb`.

* A bug (missing `drop = FALSE` in subsetting a matrix) in `plot.evpred` produced an error message if `n_years` was scalar in the prior call to `predict.evpost`.  This bug has been corrected.

* The placing of ... in the function definitions of `rpost` and `rpost_rcpp` meant that it was not possible to supply the argument `r` to be passed to `rust::ru` or `rust::ru_rcpp` to change the ratio-of-uniforms tuning parameter `r`.  Furthermore, if `model = "os"` then trying to do this sets `ros` in error. This has been corrected.

* A bug meant that the values returned by `predict(evpost_object, type = "d")` being incorrect if `evpost_object` was returned from a call to `rpost` using `model = bingp`.  The values returned were too small: they differ from the correct values by a factor approximately equal to the proportion of observations that lie above the threshold.  This bug has been corrected.

# revdbayes 1.2.0

## New features

* Faster computation, owing to the use of packages Rcpp and RcppArmadillo in package rust (https://CRAN.R-project.org/package=rust).

* New function: `rpost_rcpp`.  

* New vignette. "Faster simulation using revdbayes".

* `set_prior` has been extended so that informative priors for GEV parameters can be specified using the arguments `prior = "prob"` or `prior = "quant"`.  It is no longer necessary to use the functions `prior.prob` and `prior.quant` from the evdbayes package to set these priors.

## Bug fixes and minor improvements

* The list returned from `set_prior` now contains default values for all the required arguments of a given in-built prior, if these haven't been specified by the user.  This simplifies the evaluation of prior densities using C++.

* The GEV functions `dgev`, `pgev`, `qgev`, `rgev` and the GP functions `dgp`, `pgp`, `qgp`, `rgp` have been rewritten to conform with the vectorised style of the standard functions for distributions, e.g. those found at `?Normal`.  This makes these functions more flexible, but also means that the user take care when calling them with vectors arguments or different lengths.

* The documentation for `rpost` has been corrected: previously it stated that the default for `use_noy` is `use_noy = FALSE`, when in fact it is `use_noy = TRUE`.

* Bug fixed in `plot.evpost` : previously, in the `d = 2` case, providing the graphical parameter `col` produced an error because `col = 8` was hard-coded in a call to `points`. Now the extra argument `points_par` enables the user to provide a list of arguments to `points`.

* All the (R, not C++) prior functions described in the documentation of `set_prior` are now exported.  This means that they can now be used in the function `posterior` in the `evdbayes` package.

* Unnecessary dependence on package `devtools` via Suggests is removed.

* Bugs fixed in the (R) prior functions `gp_norm`, `gev_norm` and `gev_loglognorm`.  The effect of the bug was negligible unless the prior variances are not chosen to be large.

* In a call to `rpost` or `rpost_rcpp` with `model = "os"` the user may provide `data` in the form of a vector of block maxima.  In this instance the output is equivalent to a call to these functions with `model = "gev"` with the same data. 

# revdbayes 1.1.0

## New features

* A new vignette (Posterior Predictive Extreme Value Inference using the revdbayes Package) provides an overview of most of the new features. Run browseVignettes("revdbayes") to access.

* S3 `predict()` method for class 'evpost' performs predictive inference about the largest observation observed in N years, returning an object of class `evpred`.
  
* S3 `plot()` for the `evpred` object returned by `predict.evpost`.

* S3 `pp_check()` method for class 'evpost' performs posterior predictive checks using the bayesplot package.

* Interface to the bayesplot package added in the S3 `plot.evpost` method.

* `model = bingp` can now be supplied to `rpost()` to add inferences about the probability of threshold exceedance to inferences about threshold excesses based on the Generalised Pareto (GP) model.  `set_bin_prior()` can be used to set a prior for this probability.

* `rprior_quant()`: to simulate from the prior distribution for GEV parameters proposed in Coles and Tawn (1996) [A Bayesian analysis of extreme rainfall data. Appl. Statist., 45, 463-478], based on independent gamma priors for differences between quantiles.  
   
* `prior_prob()`: to simulate from the prior distribution for GEV parameters based on Crowder (1992), in which independent beta priors are specified for ratios of probabilities (which is equivalent to a Dirichlet prior on differences between these probabilities).

## Bug fixes and minor improvements

* The spurious warning messages relating to checking that the model argument to `rpost()` is consistent with the prior set using `set-prior()` have been corrected.  These occurred when `model = "pp"` or `model = "os"`.
  
* The hyperparameter in the MDI prior was `a` in the documentation and `a_mdi` in the code.  Now it is `a` everywhere.
  
* In `set_prior` with `prior = "beta"` parameter vector `ab` has been corrected to `pq`.
  
* In the documentation of `rpost()` the description of the argument `noy` has been corrected.
  
* Package spatstat removed from the Imports field in description to avoid NOTE in CRAN checks.  

