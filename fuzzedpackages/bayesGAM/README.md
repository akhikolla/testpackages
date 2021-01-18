# bayesGAM
Bayesian generalized additive models using Stan

The bayesGAM package is designed to provide a user friendly option to fit univariate and multivariate response Generalized Additive Models (GAM) using Hamiltonian Monte Carlo (HMC) with few technical burdens.  The R functions in this package use rstan (Stan Development Team 2020)  to call Stan routines that run the HMC simulations. The Stan code for these models is already translated to C++ and pre-compiled for the user. The programming formulation for models in bayesGAM is designed to be familiar to statisticians and analysts who fit statistical models in R.

To install, 

```r
require(remotes)
remotes::install_github("sthomas522/bayesGAM")
```


