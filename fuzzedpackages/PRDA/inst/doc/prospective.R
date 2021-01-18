## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE-----------------------------------------------------------
library(PRDA)

## ---- eval=FALSE, echo = T----------------------------------------------------
#  retrospective(effect_size, power, ratio_n = 1,
#                test_method = c("pearson", "two_sample", "welch",
#                                          "paired", "one_sample"),
#                alternative = c("two_sided","less","greater"),
#                sig_level = .05, ratio_sd = 1, B = 1e4,
#                tl = -Inf, tu = Inf, B_effect = 1e3,
#                sample_range = c(2, 1000), tol = .01,
#                display_message = TRUE)

## ---- example1----------------------------------------------------------------
set.seed(2020) # set seed to make results reproducible

prospective(effect_size = .25, power = .60, test_method = "pearson",
            display_message = TRUE)

## ---- example2----------------------------------------------------------------
prospective(effect_size = .35, power = .8, test_method = "paired",
            ratio_n = 1, display_message = FALSE)

## ---- example3----------------------------------------------------------------
prospective(effect_size = .35, power = .80, ratio_n = .5, 
            test_method = "two_sample", alternative = "great", sig_level = .10, 
            display_message = FALSE)

## ---- example4----------------------------------------------------------------
prospective(effect_size = .35, power = .80, ratio_n = .5, test_method = "welch",
            ratio_sd = 1.5, alternative = "great", sig_level = .10, 
            display_message = FALSE)

## ---- example5----------------------------------------------------------------
prospective(effect_size = function(n) rnorm(n, .3, .1), power = .60, 
            test_method = "pearson", tl = .15, tu = .45, B_effect = 500, 
            B = 500, display_message = FALSE)

## ---- data_plot---------------------------------------------------------------
da_fit <- prospective(effect_size = function(n) rnorm(n, .3, .1), power = .60,
                      test_method = "pearson", tl = .15, tu = .45, 
                      B_effect = 500, B = 500, display_message = FALSE)

str(da_fit, max.level = 1)

