# Countr 3.5.4 (CRAN)

- fixed wrong links in inst/doc/index.html.

# Countr 3.5.3

- consolidation after the refactoring of the reposotory.

- created pkgdown website.

- dealt with 'xtimes' problem on Windows (issue#1: the data frame
  returned by `optimx_2018-7.10` on Windows has a column `xtime` not
  `xtimes`).

- a number of updates in connection of the forthcoming release of the
  JSS paper.

# Countr 3.5.2 (CRAN)

- small adjustments before submission.

- renamed 'Countr.R' to 'Countr-package.R'.

- updated the main vignette with the accepted version of the JSS paper.

# Countr 3.5.0

- now covariates for ancilliary parameters for the renewal regression model can
  be specified as part of argument 'formula'. This uses the extended formulas
  provided by package `Formula`. The syntax with argument `anc` remains
  available.

# Countr 3.4.3
- dealt with warnigs from the compiler about unused variables (and a few
  others).
- When there is no precomputed bootstrap sample in the object,
  `se.coef.renewal()` now prints a slightly more informative message without
  raising a warning.

# Countr 3.4.2
- fixed the url's in the references for the vignettes in REFERENCES.bib.

# Countr 3.4.1 (2017-11-20)
- bugs fixed in convolution methods.
- option added to rescale covariates (standardise, standardise_scale).
- new vignettes added + many examples with different datasets fitted.
- Bibtex entries for citation of the vignettes added to REFERENCES.bib.
- numerous small improvements throughout the package.

# Countr 3.3.1
- new `type = "prob"` for residuals - changes in `residuals.renewal`,
     `.objectiveFunction` and `renewal()`. Comments and asome examples are
     temporarily put in the documentation of `renewal()`.
- new argument 'log' for `.objectiveFunction()`. It is now passed on to
  probability computing functions to specify whether logged values are required
  when `summa = FALSE`.

# Countr 3.3.0
- first changes on bitbucket.

# Countr 3.2.8  (2016-12-20)
- the call in the printout of renewal objects (and `summary()`) could be very
  wide (annoying in Sweave output), now `print.renewal()` and
  `print.summary.renewal()` obey 'width' when possible (see helper function
  `.deparseCall()`).
- bugs corrected in `print.summary.renewal()` when failure to compute variance
  or residuals.

# Countr 3.2.7  (2016-10-06)
- fixed some bugs in printing methods.
- added a jss_paper vignette.

# Countr 3.2.6 (2016-09-14)

# Countr 3.2.4
- merged changes done separately by T. and G. to Version 3.2.2.
- fixed some bugs introduced in 3.2.3.

# Countr 3.2.3
- all count functions vectorized in c++.
- likelihood functions improved.

# Countr 3.2.2
- minor bug fixed in `.checkInitialValues` in the anc parameters.
- all bivariate functionalities removed.
- fertility dataset now exported with data().

# Countr 3.2.1
- minor bugs fixed in `.objectiveFunction` and predict with std.error.
- last location "/..." vectorial parameters returned.
- sign in the Frank Copula package changed to be more consistent with
  the intuitive results: theta <0 -> negative dependence.

# Countr 3.2.0 (2016-03-23)
- consolidated the documentation.
- to be first CRAN version.

# Countr 3.1.x
- clean-up for CRAN submission.

# Countr 3.0.0
- renewal constructor + formula + methods.
- dePril convolution added.
- changed the way built in distributions are called.

# Countr 2.0.3
- inclusion of covariates using formula.
- Burr (weibull-gamma) distribution added.
- delayed first arrival included.

# Countr 2.0.1
- dePril convolution option added.

# Countr 2.0.0
- header files added taking advantage of Rcpp attributes.
- `conv_utils.cpp` added (contains the convolution general methods).
- weibull, gamma, gengamma counting processes added.
- option to add user based inter-arrival times distribution using routines
  `getProbs_directConv`, `getProbs_minConv`; see tests for examples.

# Countr 1.1.9
- incorporated the convolution computations from Tarak's paper.

# Countr 1.1.8
- introduced file `.Rbuildignore` in the root directory of the package.
- shortened the name of subdirectory 'Euler-van ...' to remove the complaint
  from 'R CMD build' about long paths. Also, to remove the warning by 'R CMD
  check' about 'not fully portable file names', replaced space with undescore in
  the directory name and the pdf file in it.
- in file 'DESCRIPTION' - moved Rcpp from DEPENDS to IMPORTS:.
- completed the roxigen argument descriptions (in the sense that 'R CMD check'
  doesn't complain).

# Countr 1.1.7
- `gam_weiH` and `gam_weiA` renamed `ccX` and `ccY`.
- vectorial version `dWeibullInterArrivalCountInd` and
  `dWeibullInterArrivalCountFrankCopula` improved.
- vectorial version with gamma heterogeinity added: Note that you only
  pass X %*% beta and not the exponential.
- tests added based on Examples from McShane(2008) paper.
- tests added for bivariate vectorial functions.
- `cRes` object removed from "R/" to free some space.

# Countr 1.1.5
- pkg renamed from BivCount to Countr.
- `alphagen2` changed to `alphagen`.
- oldCpp file removed.
- gamma-heterogeneity functions added to the original and fast version.
- covariates-heterogeneity functions added to the original and fast version.

# Countr 1.1.0
- renamed `alphagen` to `alphagenOrig`.
- the signatures of the cpp functions now have namespace prefixes for namespaces
  other than Rcpp. In this way RcppExports.cpp does not need `using namespace
  arma` which is not included by `Rcpp::compileAttributes` and
  `devtools::test()` and so needed to be added manually after these calls.

  I actually used search and replace to replace occurences in the body of the
  functions in our cpp files, but since we can put the 'using namespace'
  directives in them, only the signatures of exported functions are important
  for the abovde purpose.

- replaced std::cout with calls to Rprintf() in C++ code, since direct printing
  from C++ inside R can cause problems (and is caught by R CMD check).  Also
  removed <iostream.
