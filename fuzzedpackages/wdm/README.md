
wdm
===

[![Travis build status](https://travis-ci.org/tnagler/wdm-r.svg?branch=master)](https://travis-ci.org/tnagler/wdm-r) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/tnagler/wdm-r?branch=master&svg=true)](https://ci.appveyor.com/project/tnagler/wdm-r) [![Coverage status](https://codecov.io/gh/tnagler/wdm-r/branch/master/graph/badge.svg)](https://codecov.io/github/tnagler/wdm-r?branch=master) [![CRAN status](https://www.r-pkg.org/badges/version/wdm)](https://cran.r-project.org/package=wdm)

R interface to the [wdm](https://github.com/tnagler/wdm) C++ library, which provides efficient implementations of weighted dependence measures and related independence tests:

-   Pearsons's rho
-   Spearmans's rho
-   Kendall's tau
-   Blomqvist's beta
-   Hoeffding's D

All measures are computed in *O(n* log *n)* time, where *n* is the number of observations.

For a detailed description of the functionality, see the [API documentation](https://tnagler.github.io/wdm-r/).

### Installation

-   the stable release from CRAN:

``` r
install.packages("wdm")
```

-   the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
install_submodule_git <- function(x, ...) {
  install_dir <- tempfile()
  system(paste("git clone --recursive", shQuote(x), shQuote(install_dir)))
  devtools::install(install_dir, ...)
}
install_submodule_git("https://github.com/tnagler/wdm-r")
```

### Cloning

This repo contains [wdm](https://github.com/tnagler/wdm) as a submodule. For a full clone use

``` shell
git clone --recurse-submodules <repo-address>
```

### Examples

``` r
library(wdm)
```

##### Dependence between two vectors

``` r
x <- rnorm(100)
y <- rpois(100, 1)  # all but Hoeffding's D can handle ties
w <- runif(100)
wdm(x, y, method = "kendall")               # unweighted
#> [1] -0.002378163
wdm(x, y, method = "kendall", weights = w)  # weighted
#> [1] 0.08897884
```

##### Dependence in a matrix

``` r
x <- matrix(rnorm(100 * 3), 100, 3)
wdm(x, method = "spearman")               # unweighted
#>            [,1]        [,2]        [,3]
#> [1,]  1.0000000 -0.23383138 -0.05200120
#> [2,] -0.2338314  1.00000000 -0.03560756
#> [3,] -0.0520012 -0.03560756  1.00000000
wdm(x, method = "spearman", weights = w)  # weighted
#>            [,1]        [,2]        [,3]
#> [1,]  1.0000000 -0.13603502 -0.12951931
#> [2,] -0.1360350  1.00000000 -0.07473239
#> [3,] -0.1295193 -0.07473239  1.00000000
```

##### Independence test

``` r
x <- rnorm(100)
y <- rpois(100, 1)  # all but Hoeffding's D can handle ties
w <- runif(100)
indep_test(x, y, method = "kendall")               # unweighted
#>      estimate  statistic   p_value n_eff  method alternative
#> 1 -0.05555546 -0.6291182 0.5292717   100 kendall   two-sided
indep_test(x, y, method = "kendall", weights = w)  # weighted
#>      estimate  statistic   p_value   n_eff  method alternative
#> 1 -0.09298105 -0.9139628 0.3607364 75.5839 kendall   two-sided
```
