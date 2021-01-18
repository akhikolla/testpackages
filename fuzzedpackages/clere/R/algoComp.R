#' Performances SEM algorithm versus MCEM
#' 
#' This data contains four matrices corresponding to four performance indictors
#' used to compare SEM algorithm and three versions of the MCEM algorithm
#' (MCEMA: with 5 MC interations; MCEMB: with 25 MC iterations and MCEMC: 125
#' MC iterations) as described in the package vignette. The first matrix
#' \code{Pred} contains prediction errors; matrix \code{Bias} contains the bias
#' over all model parameters, matrix \code{Time} contains execution times for
#' the four methods and matrix \code{Liks} the log-likelihood reached by each
#' method. These data were used to generate the Table 1. in the package
#' vignette. For more details, please refer to the package vignette.  The R
#' script used to create this dataset is
#' clere/inst/doc/SEM_vs_MCEM_simulations.R.
#' 
#' @format A list containing four 200 x 4/5 matrices.
#' 
#' @seealso Overview : \code{\link{clere-package}} \cr 
#' Classes : \code{\linkS4class{Clere}}, \code{\linkS4class{Pacs}} \cr 
#' Methods : \code{\link{plot}}, \code{\link{clusters}}, \code{\link{predict}}, \code{\link{summary}} \cr 
#' Functions : \code{\link{fitClere}}, \code{\link{fitPacs}} 
#' Datasets : \code{\link{numExpRealData}}, \code{\link{numExpSimData}}, \code{\link{algoComp}}
#' 
"algoComp"

# library("MASS")
# library("parallel")
# library("clere")
# 
# Seed <- 1234
# set.seed(Seed)
# 
# nsim <- 200 # set nsim > 1
# plotit <- FALSE
# sparse <- FALSE
# nItEM <- 2000
# nsamp <- 5000
# analysis <- "fit"
# maxit <- 10000
# tol <- 1e-8
# n <- 25
# p <- 50
# intercept <- 0
# sigma <- 1
# gamma <- 0
# 
# g <- 3
# probs <- c(0.36 + 0.28, 0.20, 0.12 + 0.04)
# Eff <- p * probs
# a <- 4
# B <- a**(0:(g - 1)) - 1
# Beta <- rep(B, Eff)
# u <- rep(NA, nsim)
# 
# Pred <- cbind(MCEMA = u, MCEMB = u, MCEMC = u, SEM = u, ORACLE = u)
# Time <- cbind(MCEMA = u, MCEMB = u, MCEMC = u, SEM = u)
# Bias <- cbind(MCEMA = u, MCEMB = u, MCEMC = u, SEM = u)
# Liks <- cbind(MCEMA = u, MCEMB = u, MCEMC = u, SEM = u, ORACLE = u)
# 
# 
# mA <- 5
# mB <- 25
# mC <- 125
# 
# perms <- cbind(
#   c(1, 2, 3),
#   c(1, 3, 2),
#   c(2, 1, 3),
#   c(2, 3, 1),
#   c(3, 1, 2),
#   c(3, 2, 1)
# )
# 
# thetaTrue <- c(intercept, B, probs, sigma^2, gamma^2)
# ThetaA <- matrix(NA, nrow = nsim, ncol = length(thetaTrue))
# colnames(ThetaA) <- c("beta0", paste("b", 1:g, sep = ""), paste("pi", 1:g, sep = ""), "sigma2", "gamma2")
# ThetaB <- ThetaA
# ThetaC <- ThetaA
# ThetaS <- ThetaA
# 
# 
# gtheta <- function(theta) {
#   inds <- c(1, 8, 9)
#   subtheta <- theta[-inds]
#   j <- which.min(apply(perms, 2, function(j) sum((thetaTrue[-inds] - c(subtheta[1:3][j], subtheta[4:6][j]))^2)))
#   return(c(theta[1], theta[2:4][perms[, j]], theta[5:7][perms[, j]], theta[8:9]))
# }
# 
# ## Prediction
# N <- ifelse(n * nsim < 1000, 1000, n * nsim)
# Eps <- rnorm(N, mean = 0, sd = sigma)
# xpop <- matrix(rnorm(N * p), nrow = N, ncol = p)
# ypop <- as.numeric(intercept + xpop %*% Beta + Eps)
# 
# ## Parallel
# nCPU <- 10 ## was run with nCPU = 10 in the submitted article
# parallel <- TRUE
# 
# for (sim in 1:nsim) {
#   cat(paste("\tSimulation n° ", sim, ".\n", sep = ""))
#   lsim <- (1 + (sim - 1) * n):(sim * n)
#   X <- xpop[+lsim, ]
#   Y <- ypop[+lsim]
#   Xv <- xpop[-lsim, ]
#   Yv <- ypop[-lsim]
# 
#   ## MCEM A
#   mcttA <- system.time(mcmodA <- fitClere(
#     y = Y, x = X, g = g,
#     analysis = analysis, parallel = parallel, nstart = nCPU,
#     algorithm = "MCEM", seed = Seed,
#     plotit = plotit, sparse = sparse,
#     nItEM = nItEM, nBurn = max(1, 0.2 * mA),
#     nItMC = 1, nsamp = mA
#   ))
#   mcthetaA <- gtheta(c(mcmodA@intercept, mcmodA@b, mcmodA@pi, mcmodA@sigma2, mcmodA@gamma2))
#   zstartA <- clusters(mcmodA)
#   misgrpA <- which(!(1:g) %in% zstartA)
#   if (length(misgrpA)) {
#     tabA <- table(zstartA)
#     zstartA[which(zstartA == as.numeric(names(tabA[which.max(tabA)])))[1:length(misgrpA)]] <- misgrpA
#   }
#   tmpModA <- fitClere(
#     y = Y, x = X, g = g, analysis = analysis, parallel = parallel, nstart = nCPU,
#     algorithm = "SEM", plotit = plotit, sparse = sparse,
#     nItEM = 1, nBurn = 0, nItMC = 0, nsamp = nsamp, seed = Seed,
#     theta0 = mcthetaA, Z0 = zstartA - 1
#   )
#   ## MCEM B
#   mcttB <- system.time(mcmodB <- fitClere(
#     y = Y, x = X, g = g,
#     analysis = analysis, parallel = parallel, nstart = nCPU,
#     algorithm = "MCEM", seed = Seed,
#     plotit = plotit, sparse = sparse,
#     nItEM = nItEM, nBurn = max(1, 0.2 * mB),
#     nItMC = 1,
#     nsamp = mB
#   ))
#   mcthetaB <- gtheta(c(mcmodB@intercept, mcmodB@b, mcmodB@pi, mcmodB@sigma2, mcmodB@gamma2))
#   zstartB <- clusters(mcmodB)
#   misgrpB <- which(!(1:g) %in% zstartB)
#   if (length(misgrpB)) {
#     tabB <- table(zstartB)
#     zstartB[which(zstartB == as.numeric(names(tabB[which.max(tabB)])))[1:length(misgrpB)]] <- misgrpB
#   }
#   tmpModB <- fitClere(
#     y = Y, x = X, g = g, analysis = analysis, parallel = parallel, nstart = nCPU,
#     algorithm = "SEM", plotit = plotit, sparse = sparse,
#     nItEM = 1, nBurn = 0, nItMC = 0, nsamp = nsamp, seed = Seed,
#     theta0 = mcthetaB, Z0 = zstartB - 1
#   )
# 
#   ## MCEM C
#   mcttC <- system.time(mcmodC <- fitClere(
#     y = Y, x = X, g = g,
#     analysis = analysis, parallel = parallel, nstart = nCPU,
#     algorithm = "MCEM", seed = Seed,
#     plotit = plotit, sparse = sparse,
#     nItEM = nItEM, nBurn = max(1, 0.2 * mC),
#     nItMC = 1,
#     nsamp = mC
#   ))
#   mcthetaC <- gtheta(c(mcmodC@intercept, mcmodC@b, mcmodC@pi, mcmodC@sigma2, mcmodC@gamma2))
#   zstartC <- clusters(mcmodC)
#   misgrpC <- which(!(1:g) %in% zstartC)
#   if (length(misgrpC)) {
#     tabC <- table(zstartC)
#     zstartC[which(zstartC == as.numeric(names(tabC[which.max(tabC)])))[1:length(misgrpC)]] <- misgrpC
#   }
#   tmpModC <- fitClere(
#     y = Y, x = X, g = g, analysis = analysis, parallel = parallel, nstart = nCPU,
#     algorithm = "SEM", plotit = plotit, sparse = sparse,
#     nItEM = 1, nBurn = 0, nItMC = 0, nsamp = nsamp, seed = Seed,
#     theta0 = mcthetaC, Z0 = zstartC - 1
#   )
# 
#   ## SEM
#   tt2 <- system.time(mod2 <- fitClere(
#     y = Y, x = X, g = g,
#     analysis = analysis, parallel = parallel, nstart = nCPU,
#     algorithm = "SEM", seed = Seed,
#     plotit = plotit, sparse = sparse,
#     nItEM = nItEM, nBurn = nItEM / 2,
#     nItMC = 10, nsamp = nsamp
#   ))
#   theta2 <- gtheta(c(mod2@intercept, mod2@b, mod2@pi, mod2@sigma2, mod2@gamma2))
# 
#   ## oracle
#   thetaOracle <- thetaTrue
#   thetaOracle[9] <- 1e-16
#   zOracle <- rep(0:2, Eff)
#   oracle <- fitClere(
#     y = Y, x = X, g = g, analysis = analysis, parallel = parallel, nstart = nCPU,
#     algorithm = "SEM", plotit = plotit, sparse = sparse,
#     nItEM = 1, nBurn = 0, nItMC = 0, nsamp = nsamp, seed = Seed,
#     theta0 = thetaOracle, Z0 = zOracle
#   )
# 
#   Pred[sim, "MCEMA"] <- mean((Yv - predict(mcmodA, Xv))^2, na.rm = TRUE)
#   Pred[sim, "MCEMB"] <- mean((Yv - predict(mcmodB, Xv))^2, na.rm = TRUE)
#   Pred[sim, "MCEMC"] <- mean((Yv - predict(mcmodC, Xv))^2, na.rm = TRUE)
#   Pred[sim, "SEM"] <- mean((Yv - predict(mod2, Xv))^2, na.rm = TRUE)
#   Pred[sim, "ORACLE"] <- mean((Yv - predict(oracle, Xv))^2, na.rm = TRUE)
# 
#   Time[sim, "MCEMA"] <- mcttA["elapsed"]
#   Time[sim, "MCEMB"] <- mcttB["elapsed"]
#   Time[sim, "MCEMC"] <- mcttC["elapsed"]
#   Time[sim, "SEM"] <- tt2["elapsed"]
# 
#   Bias[sim, "MCEMA"] <- sum((thetaTrue - mcthetaA)^2)
#   Bias[sim, "MCEMB"] <- sum((thetaTrue - mcthetaB)^2)
#   Bias[sim, "MCEMC"] <- sum((thetaTrue - mcthetaC)^2)
#   Bias[sim, "SEM"] <- sum((thetaTrue - theta2)^2)
# 
#   Liks[sim, "MCEMA"] <- tmpModA@likelihood
#   Liks[sim, "MCEMB"] <- tmpModB@likelihood
#   Liks[sim, "MCEMC"] <- tmpModC@likelihood
#   Liks[sim, "SEM"] <- mod2@likelihood
#   Liks[sim, "ORACLE"] <- oracle@likelihood
# 
#   ThetaA[sim, ] <- mcthetaA
#   ThetaB[sim, ] <- mcthetaB
#   ThetaC[sim, ] <- mcthetaC
#   ThetaS[sim, ] <- theta2
# }
# 
# ## Output
# algoComp <- list(Bias = Bias, Pred = Pred, Time = Time, Liks = Liks)
# usethis::use_data(algoComp, overwrite = TRUE)
