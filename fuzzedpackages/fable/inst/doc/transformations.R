## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----sqrt---------------------------------------------------------------------
library(fable)
library(tsibble)
library(dplyr)
tourism %>%
  filter(Region == "Melbourne") %>% 
  model(ETS(sqrt(Trips)))

## ----combine------------------------------------------------------------------
library(tsibble)
tourism %>%
  filter(Region == "Melbourne") %>% 
  model(ETS(log(Trips + 1)))

## ----scaled-logit-------------------------------------------------------------
scaled_logit <- function(x, lower=0, upper=1){
  log((x-lower)/(upper-x))
}
inv_scaled_logit <- function(x, lower=0, upper=1){
  (upper-lower)*exp(x)/(1+exp(x)) + lower
}
my_scaled_logit <- new_transformation(scaled_logit, inv_scaled_logit)

## ----custom-transformation----------------------------------------------------
cbind(mdeaths, fdeaths) %>%
  as_tsibble(pivot_longer = FALSE) %>% 
  model(ETS(my_scaled_logit(mdeaths, 750, 3000) ~
              error("A") + trend("N") + season("A"))) %>%
  report()

