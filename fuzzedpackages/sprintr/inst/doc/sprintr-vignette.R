## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(sprintr)
set.seed(123)
n <- 100
p <- 200
x <- matrix(data = rnorm(n * p), nrow = n, ncol = p)
y <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5] + rnorm(100)

## ------------------------------------------------------------------------
mod <- sprinter(x = x, y = y, square = FALSE, nlam = 100, lam_min_ratio = 0.01)

## ------------------------------------------------------------------------
mod$idx[(p + 1) : nrow(mod$idx), ]

## ------------------------------------------------------------------------
estimate <- mod$coef[, 30]
cb <- cbind(mod$idx, estimate)
cb[cb[, 3] != 0, ]

## ------------------------------------------------------------------------
mod_cv <- cv.sprinter(x = x, y = y, square = FALSE, nlam = 100, lam_min_ratio = 0.01)

## ------------------------------------------------------------------------
mod_cv$compact

## ------------------------------------------------------------------------
newdata <- matrix(rnorm(20 * p), nrow = 20, ncol = p)
pred <- predict(mod_cv, newdata = newdata)

