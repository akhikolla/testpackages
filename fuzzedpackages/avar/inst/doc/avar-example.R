## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include=T----------------------------------------------------------
library(avar)

## ---- fig.height = 7, fig.width = 7, fig.align = 'center', fig.cap = "Allan Variance Representation."----
set.seed(2710)

# Simulate data
n = 1e5
data = rnorm(n, 0, 0.01) +  cumsum(rnorm(n, 0, 3.162278e-05))

# Compute the Maximum-Overlap Allan Variance
allan_variance = avar(data, type = "mo")

# Log-Log representation of the Allan Variance
plot(allan_variance)

## ---- warning = F--------------------------------------------------------
# Specify the scale at which we want to fit the WN and RW processes
wn = 1:7
rw = 13:15

# Compute the Allan Variance Linear Regression Estimator (AVLR)
fit = avlr(allan_variance, wn = wn, rw = rw)
fit

## ---- fig.height = 7, fig.width = 7, fig.align = 'center', fig.cap = "Empirical AV with AV implied by the latent model"----
plot(fit)
plot(fit, decomp = TRUE)

## ---- warning = F, eval = F----------------------------------------------
#  # AVLR estimator with 95% confidence intervals
#  fit_ci = avlr(allan_variance, wn = 1:7, rw = 13:15, ci = TRUE, B = 100)
#  fit_ci$ci

## ---- warning = F, echo = F, eval = T------------------------------------
load("fit_ci.rda")
fit_ci$ci

## ---- fig.height = 8, fig.width = 8, fig.align = 'center', fig.cap = "Allan Variance Representation."----
data("navchip_av")
plot(navchip_av)

## ------------------------------------------------------------------------
fit2 = avlr(navchip_av, qn_gyro = 1:4, wn_gyro = 6:8, rw_gyro = 10:13,
           wn_acc = 1:6, rw_acc = 14:16)
fit2

## ---- fig.height = 8, fig.width = 8, fig.align = 'center', fig.cap = "Empirical AV with AV implied by the latent model"----
plot(fit2)

