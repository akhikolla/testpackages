## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 5
)

## ----setup, message=TRUE------------------------------------------------------
library(gfiExtremes)

## ----gfigpd1------------------------------------------------------------------
set.seed(666L)
X <- rgpareto(200L, mu = 10, gamma = 0.3, sigma = 2)
gf1 <- gfigpd1(
  X, beta = c(0.99, 0.995, 0.999), threshold = 10, iter = 10000L,
  nchains = 4L, nthreads = 2L
)

## ----summary_gfigpd1----------------------------------------------------------
summary(gf1)
# compare with the true quantiles:
qgpareto(c(0.99, 0.995, 0.999), mu = 10, gamma = 0.3, sigma = 2)

## ----hpd_gfigpd1--------------------------------------------------------------
HPDinterval(joinMCMCchains(gf1))

## ----gfigpd2------------------------------------------------------------------
set.seed(666L)
X <- rlnorm(400L, meanlog = 1)
gf2 <- gfigpd2(
  X, beta = c(0.99, 0.995, 0.999), iter = 10000L, burnin = 10000L,
  nchains = 4L, nthreads = 2L
)
summary(gf2)
# compare with the true quantiles:
qlnorm(c(0.99, 0.995, 0.999), meanlog = 1)

