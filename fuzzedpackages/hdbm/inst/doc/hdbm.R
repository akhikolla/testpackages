## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE-------------------------------------------------------
#  install.packages("hdbm")
#  

## ---- eval = FALSE-------------------------------------------------------
#  # install.packages(devtools)
#  devtools::install_github("umich-cphds/hdbm", built_opts = c())

## ------------------------------------------------------------------------
library(hdbm)
# print just the first 10 columns
head(hdbm.data[,1:10])

## ------------------------------------------------------------------------

Y <- hdbm.data$y
A <- hdbm.data$a

# grab the mediators from the example data.frame
M <- as.matrix(hdbm.data[, paste0("m", 1:100)], nrow(hdbm.data))

# We just include the intercept term in this example.
C <- matrix(1, nrow(M), 1)

# Initial guesses for coefficients
beta.m  <- rep(0, ncol(M))
alpha.a <- rep(0, ncol(M))

set.seed(12345)
# It is recommended to pick a larger number for burnin.
hdbm.out <- hdbm(Y, A, M, C, C, beta.m, alpha.a,
                   burnin = 1000, ndraws = 100)

# Which mediators are active?
active <- which(colSums(hdbm.out$r1 * hdbm.out$r3) > 100 / 2)
colnames(M)[active]

