
<!-- README.md is generated from README.Rmd. Please edit that file -->
anMC
====

`anMC` is a R package to efficiently compute orthant probabilities of high-dimensional Gaussian vectors. The method is applied to compute conservative estimates of excursion sets of functions under Gaussian random field priors. This is an upgrade on the previously existent package [ConservativeEstimates](https://github.com/dazzimonti/ConservativeEstimates). See the paper [Azzimonti, D. and Ginsbourger D. (2018)](https://hal.archives-ouvertes.fr/hal-01289126) for more details.

### Features

The package main functions are:

-   `ProbaMax`: the main function for high dimensional othant probabilities. Computes *P(max X &gt; t)*, where *X* is a Gaussian vector and *t* is the selected threshold. The function computes the probability with the decomposition explained [here](https://hal.archives-ouvertes.fr/hal-01289126). It implements both the `GMC` and `GANMC` algorithms. It allows user-defined functions for the core probability estimate (defaults to `pmvnorm` of the package `mvtnorm`) and the truncated normal sampler (defaults to `trmvrnorm_rej_cpp`) required in the method.

-   `ProbaMin`: analogous of `ProbaMax` but used to compute *P(min X &lt; t)*, where *X* is a Gaussian vector and *t* is the selected threshold. This function computes the probability with the decomposition explained [here](https://hal.archives-ouvertes.fr/hal-01289126). It implements both the `GMC` and `GANMC` algorithms.

-   `conservativeEstimate` : the main function for conservative estimates computation. Requires the mean and covariance of the posterior field at a discretization design.

### Installation

To install the latest version of the package run the following code from a R console:

``` r
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("dazzimonti/anMC")
```

### References

Azzimonti, D. and Ginsbourger, D. (2018). Estimating orthant probabilities of high dimensional Gaussian vectors with an application to set estimation. Journal of Computational and Graphical Statistics, 27(2), 255-267. [DOI: 10.1080/10618600.2017.1360781](https://doi.org/10.1080/10618600.2017.1360781). Preprint at [hal-01289126](https://hal.archives-ouvertes.fr/hal-01289126)

Azzimonti, D. (2016). Contributions to Bayesian set estimation relying on random field priors. PhD thesis, University of Bern. Available at [link](http://biblio.unibe.ch/download/eldiss/16azzimonti_d.pdf)
