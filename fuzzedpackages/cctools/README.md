# cctools
[![Build status Linux](https://travis-ci.org/tnagler/cctools.svg?branch=master)](https://travis-ci.org/tnagler/cctools)
[![Build status Windows](https://ci.appveyor.com/api/projects/status/github/tnagler/cctools?branch=master&svg=true)](https://ci.appveyor.com/project/tnagler/cctools)
[![codecov.io](https://codecov.io/github/tnagler/cctools/coverage.svg?branch=master)](https://codecov.io/github/tnagler/cctools?branch=master)
[![CRAN version](http://www.r-pkg.org/badges/version/cctools)](https://cran.r-project.org/package=cctools)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

cctools is an R package implementing the uniform scaled beta distribution, a 
generic function for continuous convolution, and the continuous convolution
kernel density estimator, see 
[Nagler (2017)](https://arxiv.org/abs/1704.07457).

### How to install

* the stable release on CRAN:
``` r
install.packages("cctools")
```

* the latest development version:
``` r
devtools::install_github("tnagler/cctools")
```
    
### Functions

For a detailed description, see the 
[API documentation](https://tnagler.github.io/cctools/reference/index.html).

* `dusb()`, `rusb()`: Density and simulation functions for the uniform scaled
  beta distribution.
* `cont_conv()`: Expands all factor variables in a data frame and applies the 
  continuous convolution tricks to all `ordered()` variables.
* `cckde()`, `dcckde()`, `predict.cckde()`: fit and evaluate the continuous 
  convolution kernel density estimator.


### References 

Nagler, T. (2017). *A generic approach to nonparametric function estimation
with mixed data.*
[arXiv:1704.07457](https://arxiv.org/pdf/1704.07457.pdf)

