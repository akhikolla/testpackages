## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  fig.width = 6, fig.height = 5
)

## ----simplereg----------------------------------------------------------------
library(gfilinreg)
set.seed(666L)
x <- rgamma(6L, shape = 10, rate = 1)
y <- rnorm(6L, mean = x, sd = 2)
dat <- data.frame(x = x, y = y)
fidsamples <- gfilinreg(y ~ x, data = dat, distr = "normal", L = 150L)

## ----simplereg_summary--------------------------------------------------------
gfiSummary(fidsamples)

## ----simplereg_lm-------------------------------------------------------------
lmfit <- lm(y ~ x, data = dat)
coefficients(lmfit)
sigma(lmfit)
confint(lmfit)

## ----simplereg_predictions----------------------------------------------------
new <- data.frame(x = c(9, 10, 11))
fidpred <- gfilinregPredictive(fidsamples, newdata = new)
gfiSummary(fidpred)
predict(lmfit, newdata = new, interval = "prediction")

## ----simulations, eval=FALSE--------------------------------------------------
#  library(gfilinreg)
#  library(heavy)
#  library(data.table)
#  
#  nsims <- 500L
#  MAXLHD <- matrix(NA_real_, nrow = nsims, ncol = 3L)
#  colnames(MAXLHD) <- c("group1", "group2", "sigma")
#  FIDlist <- vector("list", length = nsims)
#  
#  group <- gl(2L, 6L)
#  set.seed(666L)
#  
#  for(i in 1L:nsims){
#    # simulated dataset
#    dat <- data.frame(
#      group = group,
#      y = c(rcauchy(6L), 2 + rcauchy(6L))
#    )
#    # max-likelihood estimates
#    hfit <- heavyLm(y ~ 0 + group, data = dat, family = Cauchy())
#    MAXLHD[i, ] <- c(hfit[["coefficients"]], sqrt(hfit[["sigma2"]]))
#    # fiducial stuff
#    fidsamples <- gfilinreg(y ~ 0 + group, data = dat, L = 100L, distr = "cauchy")
#    FIDlist[[i]] <- cbind(
#      parameter = c("group1", "group2", "sigma"),
#      as.data.table(gfiSummary(fidsamples))
#    )
#  }
#  FID <- rbindlist(FIDlist)

## ----load_simulations---------------------------------------------------------
library(data.table)
data("FID")
data("MAXLHD")

## ----estimate_group1----------------------------------------------------------
library(kde1d)
group1_maxlhd     <- MAXLHD[, "group1"]
group1_fid_mean   <- FID[parameter == "group1"][["mean"]]
group1_fid_median <- FID[parameter == "group1"][["median"]]

kfit_maxlhd     <- kde1d(group1_maxlhd, mult = 4)
kfit_fid_mean   <- kde1d(group1_fid_mean, mult = 4)
kfit_fid_median <- kde1d(group1_fid_median, mult = 4)

curve(
  dkde1d(x, kfit_maxlhd), from = -4, to = 4, axes = FALSE,
  lwd = 2, col = "red", lty = "dashed", xlab = "beta1", ylab = NA
)
axis(1)
curve(
  dkde1d(x, kfit_fid_mean), add = TRUE, 
  lwd = 2, col = "green", lty = "dashed"
)
curve(
  dkde1d(x, kfit_fid_median), add = TRUE, 
  lwd = 2, col = "blue", lty = "dashed"
)

## ----estimate_group2----------------------------------------------------------
group2_maxlhd     <- MAXLHD[, "group2"]
group2_fid_mean   <- FID[parameter == "group2"][["mean"]][-333L]
group2_fid_median <- FID[parameter == "group2"][["median"]][-333L]

kfit_maxlhd     <- kde1d(group2_maxlhd, mult = 4)
kfit_fid_mean   <- kde1d(group2_fid_mean, mult = 4)
kfit_fid_median <- kde1d(group2_fid_median, mult = 4)

curve(
  dkde1d(x, kfit_maxlhd), from = -2, to = 8, axes = FALSE,
  lwd = 2, col = "red", lty = "dashed", xlab = "beta2", ylab = NA
)
axis(1)
curve(
  dkde1d(x, kfit_fid_mean), add = TRUE, 
  lwd = 2, col = "green", lty = "dashed"
)
curve(
  dkde1d(x, kfit_fid_median), add = TRUE, 
  lwd = 2, col = "blue", lty = "dashed"
)

## ----estimate_sigma-----------------------------------------------------------
sigma_maxlhd     <- MAXLHD[, "sigma"]
sigma_fid_mean   <- FID[parameter == "sigma"][["mean"]][-333L]
sigma_fid_median <- FID[parameter == "sigma"][["median"]][-333L]

kfit_maxlhd     <- kde1d(sigma_maxlhd, xmin = 0, mult = 4)
kfit_fid_mean   <- kde1d(sigma_fid_mean, xmin = 0, mult = 4)
kfit_fid_median <- kde1d(sigma_fid_median, xmin = 0, mult = 4)

curve(
  dkde1d(x, kfit_maxlhd), from = 0, to = 4, axes = FALSE,
  lwd = 2, col = "red", xlab = "sigma", ylab = NA
)
axis(1)
abline(v = median(sigma_maxlhd), col = "red", lwd = 2, lty = "dashed")
curve(
  dkde1d(x, kfit_fid_mean), add = TRUE, 
  lwd = 2, col = "green"
)
abline(v = median(sigma_fid_mean), col = "green", lwd = 2, lty = "dashed")
curve(
  dkde1d(x, kfit_fid_median), add = TRUE, 
  lwd = 2, col = "blue"
)
abline(v = median(sigma_fid_median), col = "blue", lwd = 2, lty = "dashed")
# true value:
abline(v = 1, col = "yellow", lwd = 3)

## ----coverages----------------------------------------------------------------
# group1
group1 <- FID[parameter == "group1"]
mean(group1[["lwr"]] < 0)
mean(0 < group1[["upr"]])
mean(group1[["lwr"]] < 0 & 0 < group1[["upr"]])
# group2
group2 <- FID[parameter == "group2"]
mean(group2[["lwr"]] < 2)
mean(2 < group2[["upr"]])
mean(group2[["lwr"]] < 2 & 2 < group2[["upr"]])
# sigma
sigma <- FID[parameter == "sigma"]
mean(sigma[["lwr"]] < 1)
mean(1 < sigma[["upr"]])
mean(sigma[["lwr"]] < 1 & 1 < sigma[["upr"]])

