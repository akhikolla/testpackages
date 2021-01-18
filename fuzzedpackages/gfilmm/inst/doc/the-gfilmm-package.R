## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  fig.width = 5, 
  fig.height = 4.5
)

## ----simulations_linreg-------------------------------------------------------
set.seed(666L)
n <- 30L
x <- 1L:n
y <- rnorm(n, mean = x, sd = 2)
y_rounded <- round(y, digits = 1L)
dat <- data.frame(
  ylwr = y_rounded - 0.05,
  yupr = y_rounded + 0.05,
  x = x
)

## ----gfi_linreg---------------------------------------------------------------
library(gfilmm)
fidSims <- gfilmm(
  y = ~ cbind(ylwr, yupr), # interval data
  fixed = ~ x,             # fixed effects
  random = NULL,           # random effects
  data = dat,              # data
  N = 30000L,              # number of simulations
  nthreads = 2L
)

## ----gfiSummary_linreg--------------------------------------------------------
gfiSummary(fidSims)

## ----lm_linreg----------------------------------------------------------------
lmfit <- lm(y ~ x)
confint(lmfit)

## ----gfiCDF_linreg------------------------------------------------------------
Fslope <- gfiCDF(~ x, fidSims)
plot(Fslope, main = "Slope", ylab = expression("Pr("<="x)"))

## ----gfiDensity_linreg--------------------------------------------------------
library(kde1d)
kfit <- kde1d(fidSims$VERTEX[["x"]], weights = fidSims$WEIGHT, mult = 4)
curve(dkde1d(x, kfit), from = 0.7, to = 1.1, col = "red", lwd = 2)
# observe the resemblance with the distribution of the 
# frequentist estimate of the slope:
curve(
  dnorm(x, coef(lmfit)["x"], sqrt(vcov(lmfit)["x","x"])), add = TRUE, 
  col = "blue", lwd = 2, lty = "dashed"
)

## ----scatterplot--------------------------------------------------------------
indcs <- sample.int(30000L, replace = TRUE, prob = fidSims$WEIGHT)
pseudoSims <- fidSims$VERTEX[indcs,]
library(GGally)
ggpairs(
  pseudoSims,
  upper = list(continuous = ggally_density),
  lower = list(continuous = wrap("points", alpha = 0.1))
)

## ----ellipse------------------------------------------------------------------
library(car)
dataEllipse(
  pseudoSims[["(Intercept)"]], pseudoSims[["x"]], 
  levels = c(0.5,0.95), col = c("black", "red"),
  xlab = "Intercept", ylab = "Slope"
)
confidenceEllipse(
  lmfit, 1:2, levels = c(0.5,0.95), add = TRUE, 
  col = "blue", lty = "dashed"
)

## ----gfilmmPredictive---------------------------------------------------------
fpd <- gfilmmPredictive(fidSims, newdata = data.frame(x = c(1, 30)))
gfiSummary(fpd)

## ----lmpredict----------------------------------------------------------------
predict(lmfit, newdata = data.frame(x = c(1, 30)), interval = "prediction")

## ----simulations_AOV1R--------------------------------------------------------
mu           <- 10000 # grand mean
sigmaBetween <- 2
sigmaWithin  <- 3
I            <- 6L # number of groups
J            <- 5L # sample size per group

set.seed(31415926L)
groupmeans <- rnorm(I, mu, sigmaBetween)
y          <- c(
  vapply(groupmeans, function(gmean) rnorm(J, gmean, sigmaWithin), numeric(J))
)
y_rounded  <- round(y, digits = 1L)
dat        <- data.frame(
                ylwr = y_rounded - 0.05,
                yupr = y_rounded + 0.05,
                group = gl(J, I)
              )

## ----gfilmm_AOV1R-------------------------------------------------------------
fidSims <- gfilmm(
  ~ cbind(ylwr, yupr), ~ 1, ~ group, data = dat, N = 30000L, nthreads = 2L
)

## ----gfiSummary_AOV1R---------------------------------------------------------
gfiSummary(fidSims)

## ----KenwardRoger, message=FALSE----------------------------------------------
library(lmerTest)
library(emmeans)
fit <- lmer(y ~ (1|group), data = dat)
emmeans(fit, ~ 1)

## ----AOV1R--------------------------------------------------------------------
library(AOV1R)
aovfit <- aov1r(y ~ dat$group)
confint(aovfit)

## ----gfiConfInt_CV_AOV1R------------------------------------------------------
gfiConfInt(~ sqrt(sigma_group^2 + sigma_error^2)/`(Intercept)`, fidSims)

## ----parallel-----------------------------------------------------------------
gfs <- lapply(c(40000L, 50000L), function(N){
  gfilmm(~ cbind(ylwr, yupr), ~ 1, ~ group, data = dat, N = N, nthreads = 2L)  
})
lapply(gfs, gfiSummary)

