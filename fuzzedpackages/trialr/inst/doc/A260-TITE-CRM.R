## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- results='hide', message=FALSE-------------------------------------------
library(trialr)

fit <- stan_crm(skeleton = c(0.05, 0.12, 0.25, 0.40, 0.55), target = 0.25,
                doses_given = c(3, 3, 3, 3),
                tox = c(0, 0, 0, 0),
                weights = c(73, 66, 35, 28) / 126,
                model = 'empiric', beta_sd = sqrt(1.34), seed = 123)

## -----------------------------------------------------------------------------
fit

