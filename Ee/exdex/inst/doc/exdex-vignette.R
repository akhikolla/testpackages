## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(comment = "#>", collapse = TRUE)

## ------------------------------------------------------------------------
library(exdex)
theta <- spm(newlyn, 20)
# Estimates: BB2018b is BB2018 - 1/b
theta
# Estimates, SEs and bias-adjustments
summary(theta)

## ---- fig.width=7, fig.height= 6-----------------------------------------
# Sliding maxima, symmetric intervals
conf <- confint(theta)
# Sliding maxima, likelihood-based intervals
conf <- confint(theta, interval_type = "lik")
plot(conf)

## ---- fig.width=7, fig.height= 6-----------------------------------------
theta <- spm(sp500, 225)
summary(theta)
conf <- confint(theta, interval_type = "lik")
plot(conf)

## ---- fig.width=7, fig.height= 6-----------------------------------------
# Plot like the top left of Northrop (2015)
# We remove the 14 values because 2880 has lots of factors
b_vals <- c(2,3,4,5,6,8,9,10,12,15,16,18,20,24,30,32,36,40,45,48,54,60)
res <- choose_b(newlyn[1:2880], b_vals)
# Some b are too small for the sampling variance of the sliding blocks
# estimator to be estimated
plot(res, ylim = c(0, 1))

## ---- fig.width=7, fig.height= 6-----------------------------------------
b_vals <- c(10, seq(from = 25, to = 350, by = 25), 357)
res500 <- choose_b(sp500, b_vals)
plot(res500, ylim = c(0, 1))

## ------------------------------------------------------------------------
u <- quantile(sp500, probs = 0.60)
theta <- kgaps(sp500, u, k = 1)
summary(theta)

## ---- fig.width=7, fig.height= 6-----------------------------------------
u <- quantile(sp500, probs = seq(0.1, 0.9, by = 0.1))
imt_theta <- choose_uk(sp500, u = u, k = 1:5)
plot(imt_theta, uprob = TRUE, alpha = 0.05)

## ------------------------------------------------------------------------
u <- quantile(newlyn, probs = 0.90)
theta <- iwls(newlyn, u)
theta

