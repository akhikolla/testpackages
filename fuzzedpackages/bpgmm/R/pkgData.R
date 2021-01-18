# set.seed(111)
#
# source("./R/sourceR.R")
#
#
# n <- 500 # sample size
# p <- 10 # number of variables
# q <- 4 # number of factors
# K <- 10 # true number of clusters
# # sINV_diag =  # diagonal of inverse variance of errors
#
# syntheticDataset <- simData(
#   sameLambda = TRUE, sameSigma = TRUE, K.true = K, n = n, q = q, p = p,
#   sINV_values = 1 / ((1:p))
# )
#
# ##
# nsim <- 10
# burn <- 20
#
# X <- t(syntheticDataset$data)
#
# # delta = 2; ggamma = 2
# # dVec = c(1,1,1)
# # sVec = c(1,1,1)
#
# qnew <- 4
#
# Mstep <- 1
# Vstep <- 1
# constraint <- c(0, 0, 0)
#
# mInit <- 20
# mVec <- c(1, 20)
#
# res <- pgmmRJMCMC(X, mInit, mVec, qnew,
#   niter = nsim, burn = burn, constraint = constraint,
#   Mstep = Mstep, Vstep = Vstep
# )
#
# summerizePgmmRJMCMC(res, syntheticDataset$class)
