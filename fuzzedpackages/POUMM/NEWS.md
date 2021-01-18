---
title: "NEWS about the POUMM R-package"
author: "Venelin Mitov"
date: "20 October, 2020"
output: html_document
---

# POUMM 2.1.7 
* Updated dependency version of data.table (>=1.13.2). This is to always the fix https://github.com/Rdatatable/data.table/issues/4746.

# POUMM 2.1.6
* Fixed an S3-class related problem on r-devel.
* Disabled OPENMP compilation flags due to potential test failures on some platforms, e.g. Solaris.

# POUMM 2.1.5
* Improved documentation for some functions (specifyPOUMM). 
* Fixed a bug: duplicated rows for "g0" in the summary of a POUMM object.
* Compatibility with v1.12.2 of package data.table - using rep instead of recycling in columns with more than 1 but fewer than the max-number of rows (see https://github.com/venelin/POUMM/pull/2).

# POUMM 2.1.3
* Improved man page for gPOUMM function.
* Testing minor changes in the source package - store RData files in lfs. 
* Testing github release and zanodo doi-issuing. 

# POUMM 2.1.2

* Temporarily removed bivariate posterior density plotting functionality; univariate posterior density plotting still 
available. 
* Removed dependency to RcppArmadillo 
* Fixed a bug in the C++ likelihood calculation for the case when both sigmae = 0 and se = 0.
* Refactored unit tests. 
* Enabled test-coverage report with inserted badge in the README.Rmd.
* Enabled travis continuous integration with a status badge in the README.Rmd.
* Some bugs were fixed. 


# POUMM 2.1.1

* Release available from github only

# POUMM 2.1.0

* SPLiTTree dependency is renamed to SPLITT
* Release available from github only


# POUMM 2.0.0
* Using SPLiTTree for the pruning likelihood calculation.
* Automatic web-site generation using the R-package pkgdown.
* Project source-code published on github: https://github.com/venelin/POUMM.git.
* Written a new vignette "Interpreting the POUMM".

# POUMM 1.3.2
* Fixing test-errors on r-patched-solaris-x86 and notes on r-devel-linux-x86_64-fedora-clang.
* Fixed memory errors reported by ASAN and valgrind.
* Improved default values for the parameter sigmae.

# POUMM 1.3
* Parallel likelihood calculation: (tested on Linux using Intel complier v16.0.0).
* Added possibility to specify known measurement standard error for each tip of the tree. 
(see `?POUMM`).
* Added a bivariate densplot of the posterior samples displaying a pairs-plot, i.e. a 
matrix of plots with marginal densities at the diagonal, scatter plots on the
lower triangle and correlation on the upper (see `?plot.summary.POUMM`);
possibility to change the color-palette used in the plot with color-blindness
friendly default colors.
* New default parameter settings: now the ML-limits (parUpper and parLower) as well as the
priors for the parameters alpha and sigma are scaled according to the length of the tree (tMean)
and the empirical variance of the trait (zVar) (see ?specifyPOUMM). 
Please note that the default parameter settings may change in future releases. 
* Improved ML-fit: now the optimization is run 50 times starting from different points in the 
search region. The ML-fit is still kept deterministic, i.e. consecutive runs should find the same
optimum.
* Fixed a dependency on the adaptMCMC package (now added as an import).
* Update references.

# POUMM 1.2.1
* First version of the package released on CRAN.

