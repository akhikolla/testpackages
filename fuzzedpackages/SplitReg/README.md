
[![Build Status](https://travis-ci.com/AnthonyChristidis/SplitReg.svg?branch=master)](https://travis-ci.com/AnthonyChristidis/SplitReg) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/SplitReg)](https://cran.r-project.org/package=SplitReg) [![Downloads](http://cranlogs.r-pkg.org/badges/SplitReg)](https://cran.r-project.org/package=SplitReg)

SplitReg
========

This package provides functions for computing the split regularized regression estimators defined in [Christidis, Lakshmanan, Smucler and Zamar (2019)](https://arxiv.org/abs/1712.03561).

------------------------------------------------------------------------

### Installation

You can install the **stable** version on [R CRAN](https://cran.r-project.org/package=SplitReg).

``` r
install.packages("SplitReg", dependencies = TRUE)
```

You can install the **development** version from [GitHub](https://github.com/AnthonyChristidis/SplitReg)

``` r
library(devtools)
devtools::install_github("AnthonyChristidis/SplitReg")
```

### Usage

``` r
# A small example
library(MASS)
library(SplitReg)
set.seed(1)
beta <- c(rep(5, 5), rep(0, 45))
Sigma <- matrix(0.5, 50, 50)
diag(Sigma) <- 1
x <- mvrnorm(50, mu = rep(0, 50), Sigma = Sigma)
y <- x %*% beta + rnorm(50)
fit <- cv.SplitReg(x, y, num_models=10) # Use 10 models
coefs <- predict(fit, type="coefficients")
```

### License

This package is free and open source software, licensed under GPL (&gt;= 2).
