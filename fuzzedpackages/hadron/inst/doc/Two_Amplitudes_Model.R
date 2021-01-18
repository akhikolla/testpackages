## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(hadron)

## -----------------------------------------------------------------------------
scf <- bootstrap.cf(samplecf)
plot(scf, log = 'y')

## -----------------------------------------------------------------------------
fit_sample <- new_matrixfit(scf, 8, 22, model = 'single')
plot(fit_sample, log = 'y')
residual_plot(fit_sample, ylim = c(1/1.05, 1.05))

## -----------------------------------------------------------------------------
fit_sample_2 <- new_matrixfit(scf, 8, 22, model = 'two_amplitudes')
plot(fit_sample_2, log = 'y')
residual_plot(fit_sample_2, ylim = c(1/1.05, 1.05))

## -----------------------------------------------------------------------------
mapply(tex.catwitherror, fit_sample$t0, fit_sample$se, with.dollar = FALSE)
mapply(tex.catwitherror, fit_sample_2$t0, fit_sample_2$se, with.dollar = FALSE)

## -----------------------------------------------------------------------------
extent_time <- 48
time <- seq(0, extent_time - 1, by = 1)
model_E <- 0.015
model_A1 <- 0.35
model_A2 <- 0.4
val <- 0.5 * model_A1^2 * exp(-model_E * time) + 0.5 * model_A2^2 * exp(-model_E * (extent_time - time))

## -----------------------------------------------------------------------------
plot(time, val,
     main = 'Model data',
     xlab = 't',
     ylab = 'C(t)')

## -----------------------------------------------------------------------------
measurements <- do.call(cbind, lapply(val, function (v) rnorm(400, v, 0.01)))

cf <- cf_orig(cf_meta(Time = extent_time), cf = measurements)
cf <- symmetrise.cf(cf)
cf_boot <- bootstrap.cf(cf)

plot(cf_boot, log = 'y')

## -----------------------------------------------------------------------------
fit <- new_matrixfit(cf_boot, 2, 23, model = 'two_amplitudes')
plot(fit, log = 'y')
residual_plot(fit)

## -----------------------------------------------------------------------------
print(c(model_E, model_A1, model_A2))
mapply(tex.catwitherror, fit$t0, fit$se, with.dollar = FALSE)

