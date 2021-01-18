## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")

## ---- message=FALSE-----------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(PRDA)

## ---- eval=FALSE, echo = T----------------------------------------------------
#  retrospective(effect_size, sample_n1, sample_n2 = NULL,
#                test_method = c("pearson", "two_sample", "welch",
#                                "paired", "one_sample")
#                alternative = c("two_sided","less","greater"),
#                sig_level = .05, ratio_sd = 1, B = 1e4,
#                tl = -Inf, tu = Inf, B_effect = 1e3,
#                display_message = TRUE)

## ---- example1----------------------------------------------------------------
set.seed(2020) # set seed to make results reproducible

retrospective(effect_size = .25, sample_n1 = 30, test_method = "pearson")

## ---- example2----------------------------------------------------------------
retrospective(effect_size = .35, sample_n1 = 25, sample_n2 = 25,
              test_method = "paired")

## ---- example3----------------------------------------------------------------
retrospective(effect_size = .35, sample_n1 = 25, sample_n2 = 35,
              test_method = "two_sample", alternative = "great", 
              sig_level = .10, B = 1e5)

## ---- example4----------------------------------------------------------------
retrospective(effect_size = .35, sample_n1 = 25, sample_n2 = 35,
              test_method = "welch", ratio_sd = 1.5, alternative = "great", 
              sig_level = .10, B = 1e5)

## ---- example5----------------------------------------------------------------
retrospective(effect_size = function(n) rnorm(n, .3, .1), sample_n1 = 30,
              test_method = "pearson", tl = .15, tu = .45, B_effect = 1e3, 
              B = 1e3, display_message = TRUE)

## ---- data_plot---------------------------------------------------------------
da_fit <- retrospective(effect_size = function(n) rnorm(n, .3, .1), 
                        sample_n1 = 30, test_method = "pearson",
                        tl = .15, tu = .45, B_effect = 1e3, B = 1e3, 
                        display_message = FALSE)

str(da_fit, max.level = 1)

## -----------------------------------------------------------------------------
data_plot <- da_fit$retrospective_res %>%
  mutate(effect = da_fit$effect_info$effect_samples)

## ---- fig.dim= c(4, 3), dev='png'---------------------------------------------
ggplot(data_plot)+
  geom_histogram(aes(effect, y = ..density..),
                 col = "black", fill = "#00BFC4", alpha = .8,
                 breaks=seq(.15,.45,.02))+
  scale_x_continuous(breaks = seq(.1,.5,.05), limits = c(.1,.5))+
  labs(x = "Sampled Effects",
       y = "Density")+
  theme_bw()

## ---- fig.dim = c(7.23, 2.5), dev='png'---------------------------------------
data_plot %>%
  pivot_longer(cols = c("power", "typeM", "typeS"), 
               names_to = "Criteria", values_to = "Value") %>%
  mutate(Criteria = recode(Criteria, power = "Power", typeM = "Type M",  typeS = "Type S")) %>%
  ggplot(aes(x = Value, y = ..density.., fill = Criteria)) +
  geom_histogram(col = "black", alpha = .7, bins = 15) + 
  facet_wrap(.~ Criteria, scales = "free") +
  labs(y = "Density") +
  theme_bw() +
  theme(legend.position = "none") 
 
  

