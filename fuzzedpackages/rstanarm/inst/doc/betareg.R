## ---- SETTINGS-knitr, include=FALSE-------------------------------------------
stopifnot(require(knitr))
opts_chunk$set(
  comment=NA, 
  message = FALSE, 
  warning = FALSE, 
  eval = identical(Sys.getenv("NOT_CRAN"), "true"),
  dev = "png",
  dpi = 150,
  fig.asp = 0.618,
  fig.width = 5,
  out.width = "60%",
  fig.align = "center"
)

## ---- SETTINGS-gg, include=TRUE-----------------------------------------------
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default())

## ----simulated-data, fig.height=5---------------------------------------------
SEED <- 1234
set.seed(SEED)
eta <- c(1, -0.2)
gamma <- c(1.8, 0.4)
N <- 200
x <- rnorm(N, 2, 2)
z <- rnorm(N, 0, 2)
mu <- binomial(link = logit)$linkinv(eta[1] + eta[2]*x)
phi <- binomial(link = log)$linkinv(gamma[1] + gamma[2]*z)
y <- rbeta(N, mu * phi, (1 - mu) * phi)
dat <- data.frame(cbind(y, x, z))
hist(dat$y, col = "darkgrey", border = F, main = "Distribution of Outcome Variable", xlab = "y", breaks = 20, freq = F)

## ----simulated-fit, results = "hide"------------------------------------------
library(rstanarm)
fit1 <- stan_betareg(y ~ x | z, data = dat, link = "logit", link.phi = "log",
                     cores = 2, seed = 12345)
fit2 <- stan_betareg(y ~ -1 + x , data = dat, link = "logit", link.phi = "log",
                     cores = 2, seed = 12345)
round(coef(fit1), 2)
round(coef(fit2), 2)

## ----simulated-fit-print, echo=FALSE------------------------------------------
round(coef(fit1), 2)
round(coef(fit2), 2)

## ----print-priors-------------------------------------------------------------
prior_summary(fit1)

## ----simulated-analysis, fig.height=5-----------------------------------------
library(ggplot2)
library(bayesplot)
bayesplot_grid(
  pp_check(fit1), pp_check(fit2), 
  xlim = c(0,1),  
  ylim = c(0,4), 
  titles = c("True Model: y ~ x | z", "False Model: y ~ x - 1"),
  grid_args = list(ncol = 2)
)

## ----simulated-loo------------------------------------------------------------
loo1 <- loo(fit1)
loo2 <- loo(fit2)
loo_compare(loo1, loo2)

## ---- gas-fit, results="hide"-------------------------------------------------
library(rstanarm)
data("GasolineYield", package = "betareg")
gas_fit1 <- stan_betareg(yield ~ temp + batch, data = GasolineYield, link = "logit",
                         seed = 12345)
gas_fit2 <- stan_betareg(yield ~ temp + batch | pressure,
                         data = GasolineYield, link = "logit",
                         seed = 12345)
round(coef(gas_fit1), 2)
round(coef(gas_fit2), 2)

## ---- gas-print, echo=FALSE---------------------------------------------------
round(coef(gas_fit1), 2)
round(coef(gas_fit2), 2)

## ----gas-analysis, fig.height=5-----------------------------------------------
library(ggplot2)
bayesplot_grid(
  pp_check(gas_fit1), pp_check(gas_fit2), 
  xlim = c(0,1),  
  ylim = c(0,5), 
  titles = c("gas_fit1", "gas_fit2"),
  grid_args = list(ncol = 2)
)

## ---- gas-loo-----------------------------------------------------------------
gas_loo1 <- loo(gas_fit1)
gas_loo2 <- loo(gas_fit2)
loo_compare(gas_loo1, gas_loo2)

