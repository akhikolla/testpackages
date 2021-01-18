## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(SLOPE)

x <- heart$x
y <- heart$y

fit <- SLOPE(x, y, family = "binomial", lambda = "bh")

## ---- fig.cap = "Regularization path for a binomial regression model fit to the heart data set.", fig.width = 6, fig.height = 5----
plot(fit)

## -----------------------------------------------------------------------------
set.seed(924)

x <- bodyfat$x
y <- bodyfat$y

tune <- trainSLOPE(x,
                   y,
                   q = c(0.1, 0.2),
                   number = 5,
                   solver = "admm",
                   repeats = 2)

## ---- fig.cap = "Model tuning results from Gaussian SLOPE on the bodyfat dataset.", fig.width = 5.5, fig.height = 3----
plot(tune, measure = "mae") # plot mean absolute error

## -----------------------------------------------------------------------------
tune

## ---- fig.cap = "Control of false discovery rate using SLOPE.", fig.width = 4----
# proportion of real signals
q <- seq(0.05, 0.5, length.out = 20)
fdr <- double(length(q))
set.seed(1)

for (i in seq_along(q)) {
  n <- 1000
  p <- n/2
  alpha <- 1
  problem <- SLOPE:::randomProblem(n, p, q = q[i], alpha = alpha)
  
  x <- problem$x
  y <- problem$y
  signals <- problem$nonzero
  
  fit <- SLOPE(x,
               y,
               lambda = "gaussian",
               solver = "admm",
               q = 0.1,
               alpha = alpha/sqrt(n))
  
  selected_slope <- which(fit$nonzeros)
  V <- length(setdiff(selected_slope, signals))
  R <- length(selected_slope)
  fdr[i] <- V/R
}

library(lattice)
xyplot(fdr ~ q, type = "b", ylab = "FDR",
       panel = function(...) {
         panel.refline(h = 0.1)
         panel.xyplot(...)
       })

