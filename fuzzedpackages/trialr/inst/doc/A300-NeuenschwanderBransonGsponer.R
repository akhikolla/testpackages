## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(trialr)

## -----------------------------------------------------------------------------
dose <- c(1, 2.5, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 250)
d_star <- 250

## -----------------------------------------------------------------------------
target <- 0.30

## -----------------------------------------------------------------------------
outcomes <- '1NNN 2NNNN 3NNNN 4NNNN 7TT'

## -----------------------------------------------------------------------------
fit <- stan_nbg(outcome_str = outcomes, real_doses = dose, d_star = d_star,
                target = target, alpha_mean = 2.15, alpha_sd = 0.84,
                beta_mean = 0.52, beta_sd = 0.8, seed = 2020, refresh = 0)

## -----------------------------------------------------------------------------
fit

## -----------------------------------------------------------------------------
skeleton = c(0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.10, 0.17, 0.30)
fit2 <- stan_crm(outcomes, skeleton = skeleton, target = target, 
                 model = 'empiric', beta_sd = 1.34, seed = 2020, refresh = 0)
fit2

## -----------------------------------------------------------------------------
skeleton = c(0.03, 0.06, 0.12, 0.20, 0.30, 0.40, 0.50, 0.59, 0.67, 0.74)
# Obtained using dfcrm::getprior(0.05, 0.3, 5, 10)
fit3 <- stan_crm(outcomes, skeleton = skeleton, target = target, 
                 model = 'empiric', beta_sd = 1.34, seed = 2020, refresh = 0)
fit3

