## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("bama")
#  

## ---- eval = FALSE------------------------------------------------------------
#  # install.packages(devtools)
#  devtools::install_github("umich-cphds/bama", built_opts = c())

## -----------------------------------------------------------------------------
library(bama)
# print just the first 10 columns
head(bama.data[,1:10])

## -----------------------------------------------------------------------------

Y <- bama.data$y
A <- bama.data$a

# grab the mediators from the example data.frame
M <- as.matrix(bama.data[, paste0("m", 1:100)])

# We just include the intercept term in this example.
C1 <- matrix(1, nrow(M), 1)
C2 <- matrix(1, nrow(M), 1)

# Initial guesses for coefficients
beta.m  <- rep(0, ncol(M))
alpha.a <- rep(0, ncol(M))

set.seed(12345)
# It is recommended to pick a larger number for burnin.
bama.out <- bama(Y, A, M, C1, C2, beta.m, alpha.a, burnin = 3000, ndraws = 100)

# rank the mediators by PIP and return the highest 10
summary(bama.out, rank = T)[1:10,]

