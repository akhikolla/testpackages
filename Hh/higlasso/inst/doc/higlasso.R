## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = F----------------------------------------------------------------
#  install.packages("higlasso")

## ---- eval = F----------------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("umich-cphds/higlasso")

## ---- eval = T----------------------------------------------------------------
library(higlasso)

set.seed(48109)

X <- as.matrix(higlasso.df[, paste0("V", 1:7)])
Y <- higlasso.df$Y
Z <- matrix(1, nrow(X))

# This can take a bit of time, so we run the example with nondefault
# parameters to speed up the process
fit <- cv.higlasso(Y, X, Z)

print(fit)

