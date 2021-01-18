
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis build
status](https://travis-ci.com/venelin/PCMBaseCpp.svg?branch=master)](https://travis-ci.com/venelin/PCMBaseCpp)
[![Coverage
status](https://codecov.io/gh/venelin/PCMBaseCpp/branch/master/graph/badge.svg)](https://codecov.io/github/venelin/PCMBaseCpp?branch=master)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/PCMBaseCpp?color=blue)](https://cran.r-project.org/package=PCMBaseCpp)
[![Downloads](http://cranlogs.r-pkg.org/badges/PCMBaseCpp?color=blue)](https://cran.r-project.org/package=PCMBaseCpp)

# PCMBaseCpp

This is a fast C++ backend for the
[PCMBase](https://venelin.github.io/PCMBase) R-package.

# Installation

The package needs a C++ 11 compiler and Rcpp to be installed in you
R-environment. Once this is done, you can install the most recent
version of the package from github:

``` r
devtools::install_github("venelin/PCMBaseCpp")
```

If you experience problems installing the package from github, you may
try installing a possibly older version from CRAN:

``` r
install.packages("PCMBaseCpp")
```

Once the package is installed, use the function `BenchmarkRvsCpp` to
evaluate the gain in speed of the likelihood calculation on your
machine, relative to the R implementation:

``` r
library(PCMBaseCpp)
library(data.table)
options(digits = 4)

# Depending on your use case, you can change the number of traits, as well as the 
# other arguments:
benchRes <- BenchmarkRvsCpp(ks = 2, includeParallelMode = FALSE, verbose = TRUE)
# Example output:
# Performing benchmark for k:  2 ; optionSet:  serial / 1D-multiv. ...
#     k  modelType     N  R mode     logLik  logLikCpp  timeR timeCpp
#  1: 2 MGPM (A-F)    10  2   11 -7.416e+02 -7.416e+02  0.010  0.0007
#  2: 2 MGPM (A-F)   100  4   11 -4.294e+03 -4.294e+03  0.107  0.0016
#  3: 2 MGPM (A-F)  1000 11   11 -1.700e+05 -1.700e+05  1.221  0.0095
#  4: 2 MGPM (A-F) 10000 11   11 -1.210e+06 -1.210e+06 12.443  0.0795
#  5: 2     BM (B)    10  2   11 -4.451e+03 -4.451e+03  0.010  0.0003
#  6: 2     BM (B)   100  4   11 -8.427e+03 -8.427e+03  0.082  0.0008
#  7: 2     BM (B)  1000 11   11 -1.830e+04 -1.830e+04  0.847  0.0064
#  8: 2     BM (B) 10000 11   11 -6.574e+05 -6.574e+05  8.414  0.0663
#  9: 2     OU (E)    10  2   11 -1.126e+04 -1.126e+04  0.016  0.0006
# 10: 2     OU (E)   100  4   11 -8.486e+05 -8.486e+05  0.147  0.0015
# 11: 2     OU (E)  1000 11   11 -1.234e+06 -1.234e+06  1.505  0.0096
# 12: 2     OU (E) 10000 11   11 -1.058e+07 -1.058e+07 15.062  0.0854
```

For further examples, read the [Getting
started](https://venelin.github.io/PCMBaseCpp/articles/PCMBaseCpp.html)
guide and the reference available on the package
[homepage](https://venelin.github.io/PCMBaseCpp).

# Citing PCMBase

To give credit to the PCMBase package in a publication, please cite the
following articles:

Mitov, V., & Stadler, T. (2018). Parallel likelihood calculation for
phylogenetic comparative models: The SPLITT C++ library. Methods in
Ecology and Evolution, 2041–210X.13136.
<http://doi.org/10.1111/2041-210X.13136>

Mitov, V., Bartoszek, K., Asimomitis, G., & Stadler, T. (2019). Fast
likelihood calculation for multivariate Gaussian phylogenetic models
with shifts. Theor. Popul. Biol.
<https://doi.org/10.1016/j.tpb.2019.11.005>

# Used 3rd party libraries

The PCMBaseCpp R-package uses the following R-packages and C++
libraries:

  - For tree processing in C++: The [SPLITT
    library](https://venelin.github.io/SPLITT) (Mitov and Stadler 2018);
  - For data processing in R: data.table v1.12.8 (Dowle and Srinivasan
    2019);
  - For algebraic manipulation: The [Armadillo C++ template
    library](http://arma.sourceforge.net/) (Sanderson and Curtin 2016)
    and its port to R RcppArmadillo v0.9.700.2.0 (Eddelbuettel et al.
    2019);
  - For unit-testing: testthat v2.1.1 (Wickham 2019), covr v3.2.1
    (Hester 2018);
  - For documentation and web-site generation: roxygen2 v6.1.1 (Wickham,
    Danenberg, and Eugster 2018), pkgdown v1.3.0 (Wickham and
    Hesselberth 2018);

# Licence and copyright

Copyright 2016-2020 Venelin Mitov

Source code to PCMBaseCpp is made available under the terms of the GNU
General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.
PCMBaseCpp is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

# References

<div id="refs" class="references">

<div id="ref-R-data.table">

Dowle, Matt, and Arun Srinivasan. 2019. *Data.table: Extension of
‘Data.frame‘*. <https://CRAN.R-project.org/package=data.table>.

</div>

<div id="ref-R-RcppArmadillo">

Eddelbuettel, Dirk, Romain Francois, Doug Bates, and Binxiang Ni. 2019.
*RcppArmadillo: ’Rcpp’ Integration for the ’Armadillo’ Templated Linear
Algebra Library*. <https://CRAN.R-project.org/package=RcppArmadillo>.

</div>

<div id="ref-R-covr">

Hester, Jim. 2018. *Covr: Test Coverage for Packages*.
<https://CRAN.R-project.org/package=covr>.

</div>

<div id="ref-Mitov:2018dqa">

Mitov, Venelin, and Tanja Stadler. 2018. “Parallel likelihood
calculation for phylogenetic comparative models: The SPLITT C++
library.” *Methods in Ecology and Evolution*, December,
2041–210X.13136.

</div>

<div id="ref-Sanderson:2016cs">

Sanderson, Conrad, and Ryan Curtin. 2016. “Armadillo: a template-based
C++ library for linear algebra.” *Journal of Open Source Software* 1
(2).

</div>

<div id="ref-R-testthat">

Wickham, Hadley. 2019. *Testthat: Unit Testing for R*.
<https://CRAN.R-project.org/package=testthat>.

</div>

<div id="ref-R-roxygen2">

Wickham, Hadley, Peter Danenberg, and Manuel Eugster. 2018. *Roxygen2:
In-Line Documentation for R*.
<https://CRAN.R-project.org/package=roxygen2>.

</div>

<div id="ref-R-pkgdown">

Wickham, Hadley, and Jay Hesselberth. 2018. *Pkgdown: Make Static Html
Documentation for a Package*.
<https://CRAN.R-project.org/package=pkgdown>.

</div>

</div>
