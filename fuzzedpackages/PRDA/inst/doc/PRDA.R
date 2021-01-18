## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(PRDA)

## ----setup, eval = F----------------------------------------------------------
#  # If devtools is not installed yet:
#  # install.packages( "devtools" )
#  devtools::install_github("CaludioZandonella/PRDA",
#                           build_vignettes = TRUE)
#  library(PRDA)

## ---- retro-------------------------------------------------------------------
set.seed(2020) # set seed to make results reproducible

retrospective(effect_size = .25, sample_n1 = 13, test_method = "pearson")

## ---- pro---------------------------------------------------------------------
prospective(effect_size = .25, power = .8, test_method = "pearson", 
            display_message = FALSE)

