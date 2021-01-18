
<!-- README.md is generated from README.Rmd. Please edit that file -->
adpss
=====

The goal of adpss is to provide the functions for planning and conducting a clinical trial with adaptive sample size determination. Maximal statistical efficiency will be exploited even when dramatic or multiple adaptations are made. Such a trial consists of adaptive determination of sample size at an interim analysis and implementation of frequentist statistical test at the interim and final analysis with a prefixed significance level. The required assumptions for the stage-wise test statistics are independent and stationary increments and normality. Predetermination of adaptation rule is not required.

Installation
------------

You can install adpss from github with:

``` r
# install.packages("devtools")
devtools::install_github("ca4wa/R-adpss")
```

Example
-------

This is a basic example which shows you how to solve a common problem: A confirmatory randomized clinical trial is to be planned, but sample size determination is not straightforward because of scarsity of available data. In such circumstances, adaptive sample size determination is a useful option; the maximum sample size (or more generally the maximum information level) can be determined at an interim analysis without violating the prespecified significance level. In the example below, suppose that four interim analysis and one final analysis are planned. However, how many patients is required at each analysis is left unspecified in advance. The timing of each analysis will be determined adaptively. The maximum sample size at which the final analysis will be conducted will be determined at the forth interim analysis, if the trial continues beyond it without interim stopping for efficacy. This package provides a way to implement such an adaptation via the conditional error function approach with maximal statistical efficiency.

``` r
## basic example code
library(adpss)
init_work_test <- work_test_norm_global(min_effect_size = -log(0.65))
sample_size_norm_global(
  initial_test = init_work_test,
  effect_size = 11.1110 / 20.02, # effect size for which the desired level of power is ensured
  time = 20.02, # time of the forth interim analysis
  target_power = 0.75,
  sample_size = TRUE
  )
#> [1] 25.88036
```
