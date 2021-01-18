
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stratEst

stratEst is a statistical software package for strategy estimation (Dal
Bo and Frechette, 2011). The goal of strategy estimation is to explain
choices of a sample of individuals by a finite mixture of discrete
choice strategies. The discrete choice strategies are deterministic
finite state automata that can be customized by the user to fit the
structure of the data. The parameters of the strategy estimation model
are the relative frequencies and the choice parameters of the
strategies. The model can be extended by adding individual level
covariates to explain the selection of strategies by individuals. The
estimation function of the package uses expectation maximization
(Dempster, Laird, and Rubin, 1977) and Newton-Raphson methods to find
the maximum likelihood estimates of the model parameters. To speed up
the estimation, the package integrates C++ and R with the help of the R
packages Rcpp (Eddelbuettel and Francois 2011) and the open source
linear algebra library for the C++ language RppArmadillo (Sanderson and
Curtin 2016). The package contains additional functions for data
processing and simulation, strategy generation, parameter tests, model
checking, and model selection.

## Installation

The most recent CRAN version of stratEst is installed by executing the
following command in the R console:

``` r
install.packages("stratEst")
```

The development version of the package can be installed from GitHub with
the help of the package devtools (Wickham, Hester, and Chang 2020):

``` r
install.packages("devtools")
devtools::install_github("fdvorak/stratEst")
```

## Example

Fit a strategy estimation model with two strategies to the
rock-paper-scissors data of Wang, Xu, and Zhou (2014). The model is a
mixture of the Nash strategy and a strategy that imitates the last
choice.

``` r
library(stratEst)
strategies.mixture = list("nash" = strategies.RPS$nash, "imitate" = strategies.RPS$imitate)
model.mixture <- stratEst.model(data.WXZ2014,strategies.mixture)
```

## References

  - Dal Bo P, Frechette GR (2011). “The Evolution of Cooperation in
    Infinitely Repeated Games: Experimental Evidence.” American Economic
    Review, 101(1), 411-429.
  - Dempster A, Laird N, Rubin DB (1977). “Maximum Likelihood from
    Incomplete Data via the EM Algorithm.” Journal of the Royal
    Statistical Society Series B, 39(1), 1-38.
  - Eddelbuettel D, Francois R (2011). “Rcpp: Seamless R and C++
    Integration.” Journal of Statistical Software, 40(8), 1-18.
  - Sanderson C, Curtin R (2016). “Armadillo: A Template-Based C++
    Library for Linear Algebra.”
  - Wang Z, Xu B, Zhou HJ (2014). “Social Cycling and Conditional
    Responses in the Rock-Paper-Scissors Game.” Scientific Reports,
    4(1), 2045-2322.
  - Wickham H, Hester J, Chang W (2020). “devtools: Tools to Make
    Developing R Packages Easier.” R package version 2.3.0.
