
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->
genodds
=======

genodds calculates Agresti's generalized odds ratios for a two-sample dataset. This measure calculates the odds that, were a pair of observation were to be randomly selected from two groups, the observation in one group would have a higher score than those in the other group. This measure tests the hypothesis that this odds ratio is different than 1 (no difference between groups). This measure can also be reported as Number Needed to Treat, a common outcome measure used in health economics.

Installation
------------

You can install the released version of genodds from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("genodds")
```

Example
-------

genodds can be run as follows:

``` r
library(genodds)
df <- alteplase
genodds(df$mRS,df$treat,df$time)
#>   Agresti's Generalized odds ratios
#> 
#>   0-90         Odds: 0.772 (0.599, 0.996)      p=0.0468
#>   91-180       Odds: 0.850 (0.709, 1.020)      p=0.0805
#>   181-270      Odds: 0.862 (0.771, 0.963)      p=0.0088
#>   271-360      Odds: 1.019 (0.891, 1.164)      p=0.7868
#> 
#> Test of H0: odds ratios are equal among strata:
#>   X-squared = 5.63, df= 3     p=0.1311
#> 
#> Test of H0: pooled odds = 1:
#>   Pooled odds: 0.897 (0.833,0.966)  p=0.0039
```
