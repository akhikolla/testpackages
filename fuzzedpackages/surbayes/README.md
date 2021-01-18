
<!-- README.md is generated from README.Rmd. Please edit that file -->

# surbayes

<!-- badges: start -->

<!-- badges: end -->

The goal of surbayes is to provide tools for Bayesian analysis of the
seemingly unrelated regression (SUR) model. In particular, we implement
the direct Monte Carlo (DMC) approach of Zellner and Ando (2010). We
also implement a Gibbs sampler to sample from a power prior on the SUR
model.

## Installation

You can install the released version of surbayes from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("surbayes")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ethan-alt/surbayes")
```

## Example

This is a basic example which shows you how to sample from the posterior

``` r
library(surbayes)
## Taken from bayesm package
M = 10 ## number of samples
set.seed(66)
## simulate data from SUR
beta1 = c(1,2)
beta2 = c(1,-1,-2)
nobs = 100
nreg = 2
iota = c(rep(1, nobs))
X1 = cbind(iota, runif(nobs))
X2 = cbind(iota, runif(nobs), runif(nobs))
Sigma = matrix(c(0.5, 0.2, 0.2, 0.5), ncol = 2)
U = chol(Sigma)
E = matrix( rnorm( 2 * nobs ), ncol = 2) %*% U
y1 = X1 %*% beta1 + E[,1]
y2 = X2 %*% beta2 + E[,2]
X1 = X1[, -1]
X2 = X2[, -1]
data = data.frame(y1, y2, X1, X2)
names(data) = c( paste0( 'y', 1:2 ), paste0('x', 1:(ncol(data) - 2) ))
## run DMC sampler
formula.list = list(y1 ~ x1, y2 ~ x2 + x3)

## Fit models
out_dmc = sur_sample( formula.list, data, M = M )            ## DMC used
#> Direct Monte Carlo sampling used
out_powerprior = sur_sample( formula.list, data, M, data )   ## Gibbs used
#> Gibbs sampling used for power prior model
```
