# rust 1.3.10

## Bug fixes and minor improvements

* This is a patch that avoids the CRAN package check NOTEs resulting from the use of revdbayes::: instead of revdbayes::.

# rust 1.3.9

* This is a patch aimed at avoiding false positive CRAN package check warnings on some platforms.

# rust 1.3.8

## Bug fixes and minor improvements

* Dependence on the spatstat package has been removed.

* In the S3 plot method `plot.ru()` warnings, stemming from stats::regularize.values(), which is called by stats::approx(), have been suppressed.

* A very basic S3 print method `print.ru()` has been added, which simply prints the call to `ru()` or `ru_rcpp()`.

# rust 1.3.7

## Bug fixes and minor improvements

* An extra argument `shoof` has been added to `ru()` and `ru_rcpp()` to enable the user to have control over the initial value used to check that a non-zero convergence indicator is indeed spurious.

* Functions that are not intended to be called by the user have been moved to rust-internal.R

* Use `inherits()` to check the class of objects returned from `try()`, rather than `class()`.

* Hyperlinks corrected in vignettes.

* pkgdown documentation at [https://paulnorthrop.github.io/rust/](https://paulnorthrop.github.io/rust/)

# rust 1.3.6

## Bug fixes and minor improvements

* Further steps are taken in the opimisations used to set up the ratio-of-uniforms bounding box (specifically in find_a() and in cpp_find_a()) to avoid erroneous (lack of) convergence warnings. These can occur if we start the optimisation algorithm too close to the solution.

* LF line endings used in inst/include/rust.h and inst/include/rust_RcppExports.h to avoid CRAN NOTE

# rust 1.3.5

## Bug fixes and minor improvements

* The summary method for class "ru" is now set up according to Section 8.1 of the R FAQ at (https://cran.r-project.org/doc/FAQ/R-FAQ.html).

* In the extra checks used to try to avoid erroneous convergence warnings lower and upper bounds on the parameters are now used.  

* Corrected minor typo in the penultimate sentence of the "Cauchy density" subsection in the "When can rust be used vignette": b1* should be b1+.

# rust 1.3.4

## Bug fixes and minor improvements

* Removed unnecessary dependence on the current release of R that caused an installation error on r-oldrel.

# rust 1.3.3

## New features

* New vignette. "When can rust be used?".

# rust 1.2.3

## Bug fixes and minor improvements

* In `ru()` and `ru_rcpp()` the constant `hscale` is used to shift `logf` (and hence scale the target density f) in function `logf` in the returned object.  This helps to avoid over/under-flow when contouring f in `plot.ru` when `d = 2`.

* The `var_names` argument to `ru_rcpp` didn't work.  This has been corrected.

* The arguments `lower, upper` to `ru` and `ru_rcpp` are now used, at least partly even when `rotate = TRUE`.  See the updated description of `lower, upper` in the documentation of `ru` and `ru_rcpp`.

* That the function `logf` supplied to `ru` or `ru_rcpp` should return `-Inf` when the density f is zero is stated explicitly in the help files.

* ru() did not work when `trans = "user"` and `d` > 1. This has been corrected.

* Extra checks are used to try to avoid erroneous convergence warnings.

* Incorrectly formatted link to the Rcpp Gallery webpage corrected in the vignette "Rusting Faster: Simulation using Rcpp".

* Extra examples provided for `ru` and for `ru_rcpp`: (a) Cauchy, (b) half-Cauchy and (c) bivariate normal x bivariate student-t.

# rust 1.2.2

## Bug fixes and minor improvements

* An overloading ambiguity has been corrected to ensure installation on Solaris.

# rust 1.2.1

## Bug fixes and minor improvements

* Corrected C++ function `vecpow` to avoid compilation errors on some platforms.

* Unnecessary dependence on packages `devtools` and `roxygen2` via Suggests is removed.

* Minor edit to vignette: provide link directly to example C++ file `user_fns.cpp` in `src` directory, rather than the the (identical) `example_user_fns.cpp` file in the `vignettes` directory.

# rust 1.2.0

## New features

* Packages Rcpp (https://CRAN.R-project.org/package=Rcpp) and RcppArmadillo (https://CRAN.R-project.org/package=RcppArmadillo) are used to speed up the computations if the user provides a C++ function to evaluate their target log-density. 

* New functions: `ru_rcpp`, `find_lambda_rcpp` and `find_lambda_one_d_rcpp`.  

* New vignette. "Rusting faster: Simulation using Rcpp".

## Bug fixes and minor improvements

* Bug fixed in `plot.ru()`: previously `plot.ru()` failed when `d > 2`  and no axis label names were provided.

* Bug fixed in `plot.ru` : previously, in the `d = 2` case, providing the graphical parameter `col` produced an error because `col = 8` was hard-coded in a call to `points`. Now the extra argument `points_par` enables the user to provide a list of arguments to `points`.

* "using `pairs()`" removed from the last sentence Description of `plot.ru()` because `pairs()` is not used when `d > 2`, rather a single plot is produced for each pair of variables.

* Obsolete function `rho_to_theta()` removed from function `ru` in `ru_sampling.R`.

* If the user calls `ru` (or `ru_rcpp`) with `trans = "user"` but doesn't supply `phi_to_theta` then an error is returned.

* `plot.ru` edited to avoid warning message that occurred in the `d=1` case when `breaks` was supplied as an argument.

* The functions `rgpd`, `gpd_sum_stats`, `gpd_init` and `gpd_logpost` are now exported.

# rust 1.1.0

## New features

* `plot.ru()` can now be used when `d > 2`: it produces pairwise plots of the simulated values.

* `find_lamba()`: argument `init_lambda` added to enable the user to supply an initial value for the Box-Cox transformation vector `lambda`.

## Bug fixes and minor improvements

* Unnecessary print statement `print(dim(phi))` removed from function `find_lambda()`.

* Unnecessary print statement `print(a_algor)` removed from function `ru()`.

* Correct `lambda$init` to `lambda$init_psi` in `ru()` when extracting Box-Cox information.
   
* Documentation of `ru()` updated to include a description of the returned function `logf_rho()` and simulated values `sim_vals_rho` and to clarify the meaning of the returned value of `f_mode`.

* `ru()`: the expression for the inverse Box-Cox transformation in the case where lambda is exactly 0 has been corrected. 

* `find_lambda()`: carry out calculation of the target on a shifted log scale to avoid underflow.

* Set up `plot.ru()` so that if the user supplies axis labels then they are used and otherwise the column name(s) of `ru_object$sim_vals` are used. Also enable plotmath symbols to be rendered in the axis labels.
