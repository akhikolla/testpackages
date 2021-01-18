
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rust <img src="tools/rust_logo.png" align="right" />

[![Travis-CI Build
Status](https://travis-ci.org/paulnorthrop/rust.svg?branch=master)](https://travis-ci.org/paulnorthrop/rust)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/paulnorthrop/rust?branch=master&svg=true)](https://ci.appveyor.com/project/paulnorthrop/rust)
[![Coverage
Status](https://codecov.io/github/paulnorthrop/rust/coverage.svg?branch=master)](https://codecov.io/github/paulnorthrop/rust?branch=master)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/rust)](https://cran.r-project.org/package=rust)

## Ratio-of-uniforms simulation with transformation

### What does rust do?

The `rust` package implements the multivariate generalized
ratio-of-uniforms method of simulating random variates from a
d-dimensional continuous distribution. The user specifies (the log of) a
positive target function `f` that is proportional to the density
function of the distribution.

### A simple example

We use the function `ru` to simulate a sample of size 1000 from a
two-dimensional standard normal distribution with strong positive
correlation between the components. Of course, this particular example
is purely illustrative: there are better ways to simulate from a
multivariate normal distribution.

``` r
rho <- 0.9
covmat <- matrix(c(1, rho, rho, 1), 2, 2)
log_dmvnorm <- function(x, mean = rep(0, d), sigma = diag(d)) {
  x <- matrix(x, ncol = length(x))
  d <- ncol(x)
  - 0.5 * (x - mean) %*% solve(sigma) %*% t(x - mean)
}
x <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1000, init = c(0, 0))
```

From version 1.2.0 onwards the faster function `ru_rcpp` can be used.
See the vignette “Rusting Faster: Simulation using Rcpp” for details.

``` r
# Create an external pointer to a C++ function to evaluate the log-density.
ptr_bvn <- create_xptr("logdnorm2")
# Pass the external pointer to `ru_rcpp`.
x <- ru_rcpp(logf = ptr_bvn, rho = rho, d = 2, n = 1000, init = c(0, 0))
```

### Installation

To get the current released version from CRAN:

``` r
install.packages("rust")
```

### Vignettes

See `vignette("rust-a-vignette", package = "rust")` for an overview of
the package, `vignette("rust-b-when-to-use-vignette", package = "rust")`
for guidance on when `rust` can be used and
`vignette("rust-c-using-rcpp-vignette", package = "rust")` for
information on how to take advantage of the Rcpp package.
