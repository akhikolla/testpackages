
<!-- README.md is generated from README.Rmd. Please edit that file -->
Valuer - VA value with R
========================

Valuer aims at pricing a type of life insurance contract called variable annuity. The package implements the valuation framework and algorithms described in [BMOP2011](#BMOP2011) where Monte Carlo methods are adapted to the life insurance case. It's written using [R6](https://CRAN.R-project.org/package=R6) and comes with classes which describe the variable annuity contracts and other classes, called pricing engines, which are used to price those contracts.

### Example

The following code prices a 10 years VA with GMAB guarantee by means of a pricing engine which models the underlying fund as a geometric Brownian motion. Please check the introductory [vignette](https://CRAN.R-project.org/package=valuer) for an explanation of this example and a description of the package structure.

``` r
library(valuer)
#> Loading required package: orthopolynom
#> Loading required package: polynom

rate <- constant_parameters$new(0.01)

premium <- 100
rollup <- payoff_rollup$new(premium, rate)

#Ten years time-line
begin <- timeDate::timeDate("2016-01-01")
end <- timeDate::timeDate("2025-12-31")

#Age of the policyholder.
age <- 60
# A constant fee of 4% per year (365 days)
fee <- constant_parameters$new(0.04)

#Barrier for a state-dependent fee. The fee will be applied only if
#the value of the account is below the barrier
barrier <- Inf
#Withdrawal penalty applied in case the insured surrenders the contract
#It is a constant penalty in this case
penalty <- penalty_class$new(type = 1, 0.01)
#Sets up the contract with GMAB guarantee
contract <- GMAB$new(rollup, t0 = begin, t = end, age = age, fee = fee, barrier = barrier, penalty = penalty)
#Interest rate
r <- constant_parameters$new(0.03)
#Initial value of the underlying fund
spot <- 100
#Volatility
vol <- constant_parameters$new(0.2)
#Dividend rate
div <- constant_parameters$new(0.0)
#Gatherer for the MC point estimates
the_gatherer <- mc_gatherer$new()
#Number of paths to simulate
no_of_paths <- 1e3

#Sets up the pricing engine specifying the va_contract, the interest rate
#the parameters of the Weibull intensity of mortality, the initial fund
#value, the volatility and dividends rate
engine <- va_bs_engine$new(contract, r, c1=90.43, c2=10.36, spot,
volatility=vol, dividends=div)

#Estimates the contract value by means of the static approach.

engine$do_static(the_gatherer, no_of_paths)
the_gatherer$get_results()
#>       mean       se
#> 1 91.84034 1.009428
```

### Release status

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/valuer)](https://cran.r-project.org/package=valuer)

### Installation

Get valuer from CRAN:

``` r
install.packages("valuer")
```

Get the development release from GitHub:

``` r
# install.packages("devtools")

devtools::install_github("IvanZoccolan/valuer")
```

### Build status

[![Travis-CI Build Status](https://travis-ci.org/IvanZoccolan/valuer.svg?branch=master)](https://travis-ci.org/IvanZoccolan/valuer)

### References

<a name="BMOP2011"></a> BMOP2011 - Bacinello A.R., Millossovich P., Olivieri A. e Pitacco E. "Variable annuities: unifying valuation approach." In: Insurance: Mathematics andEconomics 49 (2011), pp. 285-297.
