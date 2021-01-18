## ----echo=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(
    message = FALSE,
    warning = FALSE,
    error = FALSE,
    tidy = FALSE,
    cache = FALSE
)

## -----------------------------------------------------------------------------
library(cubature)
m <- 3
sigma <- diag(3)
sigma[2,1] <- sigma[1, 2] <- 3/5 ; sigma[3,1] <- sigma[1, 3] <- 1/3
sigma[3,2] <- sigma[2, 3] <- 11/15
logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
my_dmvnorm <- function (x, mean, sigma, logdet) {
    x <- matrix(x, ncol = length(x))
    distval <- stats::mahalanobis(x, center = mean, cov = sigma)
    exp(-(3 * log(2 * pi) + logdet + distval)/2)
}

## -----------------------------------------------------------------------------
cubintegrate(f = my_dmvnorm, lower = rep(-0.5, 3), upper = c(1, 4, 2), method = "pcubature",
             mean = rep(0, m), sigma = sigma, logdet = logdet)

## -----------------------------------------------------------------------------
cubintegrate(f = my_dmvnorm, lower = rep(-0.5, 3), upper = c(1, 4, 2), method = "cuhre",
             mean = rep(0, m), sigma = sigma, logdet = logdet)

## -----------------------------------------------------------------------------
cubintegrate(f = my_dmvnorm, lower = rep(-0.5, 3), upper = c(1, 4, 2), method = "cuhre",
             mean = rep(0, m), sigma = sigma, logdet = logdet,
             flags = list(verbose = 2))

## -----------------------------------------------------------------------------
default_args()

## -----------------------------------------------------------------------------
my_dmvnorm_v <- function (x, mean, sigma, logdet) {
    distval <- stats::mahalanobis(t(x), center = mean, cov = sigma)
    exp(matrix(-(3 * log(2 * pi) + logdet + distval)/2, ncol = ncol(x)))
}

## -----------------------------------------------------------------------------
cubintegrate(f = my_dmvnorm_v, lower = rep(-0.5, 3), upper = c(1, 4, 2), method = "pcubature",
             mean = rep(0, m), sigma = sigma, logdet = logdet,
             nVec = 128)

## -----------------------------------------------------------------------------
cubintegrate(f = my_dmvnorm_v, lower = rep(-0.5, 3), upper = c(1, 4, 2), method = "cuhre",
             mean = rep(0, m), sigma = sigma, logdet = logdet,
             nVec = 128)

## -----------------------------------------------------------------------------
sessionInfo()

