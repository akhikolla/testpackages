# code for revision
# library("SASxport")
# library("matrixcalc")
# library("GSE")


## The cellHandler method
#########################


cellHandler <- function(X, mu, Sigma, quant = 0.99) {
  # This is the cellFlagger algorithm starting from a given mu/Sigma. 
  # It returns a binary matrix W where 1 indicates a flagged & imputed cell.
  # Ximp is the imputed matrix, where the imputed values were RECALCULATED
  # using EM after W was constructed.
  # 
  # Arguments:
  #   X:     data matrix
  #   mu:    center
  #   Sigma: covariance matrix
  #   quant: quantile used for flagging cells
  # Returns: 
  #   Ximp:       imputed data
  #   W:          matrix indicating which cells were flagged/imputed
  #               (1 = flagged & imputed)
  #   Zres:       matrix with cellwise residuals for the flagged cells
  #   cellPaths:  matrix with paths of individual regressions
  #   cellOrder:  matrix with order of each cell in its path
  #   Zres_num:   numerator of the standardized residuals
  #   Zres_denom: denominator of the standardized residuals
  #
  n <- dim(X)[1]
  d <- dim(X)[2]
  inv.out    <- mpinv(Sigma)
  Sigmai     <- inv.out$Inv
  Sigmaisqrt <- inv.out$InvSqrt
  scales     <- sqrt(diag(Sigma))
  predictors <- Sigmaisqrt
  Ximp <- X
  naMask <- is.na(X) + 0
  W <- matrix(0, n, d)
  Zres <- matrix(0, n, d)
  cellPaths <- Zres_num <- Zres_denom <- matrix(0, n, d)
  #
  for (i in 1:n) {
    x        <- X[i, ] 
    indNA    <- which(naMask[i, ] == 1)
    x[indNA] <- mu[indNA]
    response <- Sigmaisqrt %*% (x - mu)
    weights  <- huberweights(x = (x - mu) / scales, b = 1.5)
    larOut   <- findCellPath_cpp(predictors = predictors,
                                            response = response,
                                            weights = weights,
                                            Sigmai = Sigmai,
                                            naMask = naMask[i, ])
    cellPaths[i, ] <- larOut$ordering
    deltas   <- abs(diff(larOut$RSS))
    badCells <- which(deltas > qchisq(quant, 1))
    if (length(indNA) > 0) {
      badCells <- unique(indNA, badCells)
    }
    #
    if (length(badCells) > 0) {
      badinds  <- larOut$ordering[1:max(badCells)] # maxDelta rule
      # now calculate residuals:
      if (length(badinds) == d) {
        stdresid <- (x - mu) / sqrt(diag(Sigma))
      } else {
        stdresid <- rep(0, d)
        stdresid[badinds] <- abs(larOut$beta[length(badinds) + 1, badinds ]) / 
          sqrt(diag(Sigma[badinds, badinds] - Sigma[badinds, -badinds] %*%
                      solve(Sigma[-badinds, -badinds]) %*% Sigma[-badinds, badinds]))
      }
      badinds <- which(abs(stdresid) > sqrt(qchisq(quant, 1)))
      
      if (length(indNA) > 0) {
        badinds <- unique(indNA, badinds)
      }
      
      W[i, badinds] <- 1
      #
      if (length(badinds) == d) {
        Ximp[i, ] <- mu
        Zres_num[i, ] <- (X[i, ] - mu)
        Zres_denom[i, ] <- sqrt(diag(Sigma))
        Zres[i, ] <-  Zres_num[i, ] / Zres_denom[i, ]
      } else {
        replacement <- X[i, ]
        replacement[badinds] <- mu[badinds] +  Sigma[badinds, -badinds] %*%
          solve(Sigma[-badinds, -badinds]) %*% (replacement[-badinds] - mu[-badinds])
        Ximp[i, ] <- replacement
        residual  <- X[i, ] - replacement
        Zres_num[i, badinds]   <- residual[badinds]
        Zres_denom[i, badinds] <- sqrt(diag(Sigma[badinds, badinds] - Sigma[badinds, -badinds] %*%
                                              solve(Sigma[-badinds, -badinds]) %*% Sigma[-badinds, badinds]))
        Zres[i, badinds] <- Zres_num[i, badinds] / Zres_denom[i, badinds] 
      }
    }
  }
  
  indcells <- which(abs(Zres) > sqrt(qchisq(quant, 1)))
  # W[badinds] <- 1
  # 
  # indcells <- which(W == 1)
  indNAs   <- which(naMask == 1)
  indcells <- setdiff(indcells, indNAs) # take NA inds out of indcells
  
  return(list(Ximp = Ximp, indcells = indcells,
              indNAs = indNAs, Zres = Zres,
              cellPaths = cellPaths,
              Zres_denom = Zres_denom))
}



mpinv <- function(X, tol = sqrt(.Machine$double.eps)) {
  ## Moore-Penrose generalized inverse of a matrix and its square root.
  ## This is the function ginv from the package MASS.
  #
  # Arguments:
  #   X:   a square matrix
  #   tol: tolerance parameter
  # Returns:
  #   Inv: Moore-Penrose inverse of X
  #   InvSqrt: square root of Inv
  #
  dnx = dimnames(X)
  if (is.null(dnx)) dnx <- vector("list", 2)
  s = svd(X)
  nz = s$d > tol * s$d[1]
  if (any(nz)) {
    outInv <- s$v[, nz] %*% (t(s$u[, nz])/s$d[nz])
    outInvsqrt <- s$v[, nz] %*% (t(s$u[, nz])/sqrt(s$d[nz]))
  } else { 
    outInv <- outInvsqrt <- X
  }
  return(list(Inv = outInv, InvSqrt = outInvsqrt))
}


huberweights <- function(x, b) {
  # vectorized Huber weight function
  # 
  # Arguments:
  #   x: a vector
  #   b: tuning constant
  # Returns:
  #   weights
  #
  highind <- which(abs(x) > b)
  result  <- rep(1, length(x))
  result[highind] <- abs(b / x[highind])
  return(result)
}


## Initial estimators:

DDCWcov <- function(X, maxCol = 0.25) {
  # computes a cellwise robust initial estimator for the covariance matrix
  # based on DetectDeviatingDatacells in combination with Wrapping
  # 
  # Arguments:
  #   X: data matrix
  #   maxcol: maximum number of flagged cells per column
  # Returns:
  #   locScale: robust estimates of the location and scale of each variable
  #             in X
  #   mu:       initial estimate of location for the standardized data
  #   cov:      initial estimate of the covariance matrix for the
  #             standardized data
  #   Z:        the robustly standardized data
  #   indrows:  indices of the outlying rows
  #
  # Functions only used here:
  
  DDC_controlled <- function(X, tolProbCell, maxCol = 0.25) {
    # Executes DDC with a given tolProbCell, which ensures that
    # no more than maxCol*n cells are flagged in any variable
    n <- dim(X)[1]
    d <- dim(X)[2]
    DDCout <- DDC(X, list(fastDDC = FALSE, silent = TRUE,
                                    tolProbCell = tolProbCell, standType = "wrap"))
    Wna <- matrix(0, n, d); Wna[DDCout$indcells] <- 1
    overflag <- which(colSums(Wna) > maxCol * n)
    if (length(overflag) > 0) {
      for (i in 1:length(overflag)) {
        ind <- overflag[i]
        replacement <- rep(0, n)
        replacement[order(abs(DDCout$stdResid[, ind]),
                          decreasing = TRUE)[1:(floor(maxCol*n))]] <- 1
        Wna[, ind] <- replacement
      }
      DDCout$indcells <- which(Wna == 1)
      DDCout$Ximp <- X; DDCout$Ximp[DDCout$indcells] <- DDCout$Xest[DDCout$indcells]
    }
    return(DDCout)
  }
  
  iDDC9I.O.Wrap <- function(X, maxCol = 0.25) {
    n <- dim(X)[1]
    d <- dim(X)[2]
    DDCout   <- DDC_controlled(X, tolProbCell = 0.9, maxCol = maxCol)
    locScale <- list(loc = DDCout$locX, scale = DDCout$scaleX)
    Z        <- scale(X, locScale$loc, locScale$scale)
    Zimp     <- scale(DDCout$Ximp, locScale$loc, locScale$scale)
    Zorig    <- Z
    Zimporig <- Zimp
    Znaorig  <- Z; Znaorig[DDCout$indcells] <- NA
    if (length(DDCout$indrows) > 0) {
      Z                   <- Z[-DDCout$indrows, ]
      Zimp                <- Zimp[-DDCout$indrows, ]
    }
    # Calculate first eigenvector estimate on imputed data, project and estimate scale
    eigenvectors  <- eigen(cov(Zimp), symmetric = TRUE)$vectors
    Zimp_proj     <- Zimp %*% eigenvectors
    locscale_proj <- estLocScale(Zimp_proj)
    #
    # Calculate final estimate
    Zimp_proj_w <- wrap(Zimp_proj, locscale_proj$loc, locscale_proj$scale)$Xw
    cov         <- eigenvectors %*% cov(Zimp_proj_w) %*%  t(eigenvectors)
    cov         <- t(cov2cor(cov) * locScale$scale) * locScale$scale
    return(list(locScale = locScale,
                mu = rep(0, dim(X)[2]),
                cov = cov,
                Z = Z,
                Zorig = Zorig,
                Zimporig = Zimporig,
                Znaorig = Znaorig,
                indrows = DDCout$indrows))
  }
  
  RR <- function(Z, cov, b = 2, quant = 0.99) {
    d <- dim(Z)[2]
    MDs <- mahalanobis(pmin(pmax(Z, -b), b), rep(0, d), cov)
    rowinds <- which(MDs / median(MDs) * qchisq(0.5, d) > qchisq(quant, d))
    return(rowinds)
  }
  
  iDDC9I.O.Wrap.RR <- function(X, maxCol = 0.25) {
    d <- dim(X)[2]
    result   <- iDDC9I.O.Wrap(X, maxCol = maxCol)
    cov      <- result$cov
    locScale <- result$locScale
    Z        <- result$Zorig
    Zimp     <- result$Zimporig
    Zna      <- result$Znaorig
    rowinds <- RR(Z, cov, 2, 0.99)
    if (length(rowinds) > 0) {
      Z <- Z[-rowinds, ]
      Zimp <- Zimp[-rowinds,]
    }
    return(list(locScale = locScale,
                mu = rep(0, dim(X)[2]),
                cov = cov,
                Z = Z,
                Zimp = Zimp,
                Zna = Zna,
                indrows = rowinds))
  }
  
  d <- dim(X)[2]
  result   <- iDDC9I.O.Wrap.RR(X, maxCol = maxCol)
  cov      <- result$cov
  locScale <- result$locScale
  Z        <- result$Z
  Zimp     <- result$Zimp
  indrows  <- result$indrows
  #
  # orthogonalize
  eigenvectors  <- eigen(cov(Zimp), symmetric = TRUE)$vectors
  Zimp_proj     <- Zimp %*% eigenvectors
  locscale_proj <- estLocScale(Zimp_proj)
  #
  # Wrap and Calculate final estimate
  Zimp_proj_w <- wrap(Zimp_proj, locscale_proj$loc, locscale_proj$scale)$Xw
  Zcov.raw <- eigenvectors %*% cov(Zimp_proj_w) %*%  t(eigenvectors)
  Zcov <- cov2cor(Zcov.raw)
  cov <- t(t(Zcov) * locScale$scale) * locScale$scale
  
  return(list(center = locScale$loc,
              cov = cov, 
              locScale = locScale,
              Z = Z,
              Zcov = Zcov,
              Zcov.raw = Zcov.raw))
}


TwoSGS <- function(X) {
  # Wrapper around GSE::TSGS
  # 
  tsgs.out <- GSE::TSGS(X)
  locScale <- list(loc = tsgs.out@mu, scale = sqrt(diag(tsgs.out@S)))
  Zcov      <- cov2cor(tsgs.out@S)
  Z        <- scale(X, locScale$loc, locScale$scale)
  return(list(locScale = locScale,
              Zcov = Zcov,
              Zcov.raw = Zcov,
              cov = tsgs.out@S,
              Z = Z))
}


## The Detection Imputation (DI) method:


DI = function(X,
              initEst = "DDCWcov",
              crit = 0.01,
              maxits = 10,
              quant = 0.99,
              maxCol = 0.25,
              checkPars = list()){
  # Computes a covariance matrix on data with possibly
  # both cellwise and casewise outliers.
  # 
  # Arguments:
  #   X:          data matrix
  #   initEst:    the initial estimator used
  #   crit:       criterion for convergence
  #   maxits:     maximum number of iterations
  #   quant:      quantile used for flagging cells
  #   maxCol:     maximum number of flagged cells per column
  # Returns: 
  #   center:      final estimate of the center
  #   cov:         final estimate of the covariance matrix
  #   center_init: initial estimate of the center
  #   cov_init:    initial estimate of the covariance matrix
  #   allSigmas:   3-dim array with the covariance estimates of every iteration
  #   Ximp:        imputed data
  #   nbimps:      final number of imputes for each column
  #   W:           matrix indicating which cells were imputed
  #
  # DI uses its own version of CellFlagger since it needs to take 
  # into account maxCol, biascorrections, etc.
  # 
  
  # Check inputs
  
  if (!is.data.frame(X) & !is.matrix(X)) {
    stop("The input data must be a matrix or a data frame")
  }
  X <- as.matrix(X)
  
    # parameters for checkDataSet
  if (!"coreOnly" %in% names(checkPars)) {
    checkPars$coreOnly <- FALSE
  }
  if (!"silent" %in% names(checkPars)) {
    checkPars$silent <- FALSE
  }
  if (!"numDiscrete" %in% names(checkPars)) {
    checkPars$numDiscrete <- 5
  }
  if (!"precScale" %in% names(checkPars)) {
    checkPars$precScale <- 1e-12
  }
  if (!"fracNA" %in% names(checkPars)) {
    checkPars$fracNA <- 0.15
  }
  
  CD_out <- list()
  if (!checkPars$coreOnly) {
    # Check the data set and set aside columns and rows that do
    # not satisfy the conditions:
    CD_out <- checkDataSet(X,
                        fracNA = checkPars$fracNA,
                        numDiscrete = checkPars$numDiscrete,
                        precScale = checkPars$precScale,
                        silent = checkPars$silent)
    X <- CD_out$remX
  } 
  
  mS_cov <- function(X, distances, nbimps) {
    return(list(mu = colMeans(X), cov = cov(X)))
  }
  
  # Step 1A: initial estimate & standardization
  if (is.list(initEst)) { # initEst is list with mu and Sigma
    locScale_init <- list(loc = initEst$mu,
                          scale = sqrt(diag(initEst$Sigma)))
    mu_init       <- rep(0, dim(X)[2])
    cov_init      <- cov2cor(initEst$Sigma)
    Z             <- scale(X, initEst$mu, sqrt(diag(initEst$Sigma)))
    out_init      <- list()
  } else {
    if (initEst == "TSGS") {
      initEst <- TwoSGS
    } else {
      initEst <- DDCWcov
    }
    out_init      <- initEst(X)
    locScale_init <- out_init$locScale # locScale of X
    mu_init       <- rep(0, ncol(X)) # initial location of Z (should be zeroes)
    cov_init      <- out_init$Zcov.raw # initial cov of Z
    Z             <- out_init$Z
  }
  
  # Step 1B: initialization
  nbits      <- 0
  convcrit   <- 1
  n          <- dim(Z)[1]
  d          <- dim(Z)[2]
  M          <- floor(maxCol * n) # max number of imputed cells per variable
  mu         <- mu_init
  Sigma      <- cov_init
  invOut     <- mpinv(cov_init)
  Sigmai     <- invOut$Inv
  Sigmaisqrt <- invOut$InvSqrt
  Zimp       <- Z
  naMask     <- is.na(Z) + 0
  
  # Step 1C: containers for simulation
  mus           <- array(0, dim = c(maxits + 1, d))
  Sigmas        <- array(0, dim = c(maxits + 1, d, d))
  mus[1,]       <- mu_init
  Sigmas[1, , ] <- cov_init
  
  # Step 2: iteration step
  while ((nbits < maxits) && (convcrit > crit)) {
    # Step 2A: flag 
    
    # 2A Stage 1: Univariate regressions
    predictors <- Sigmaisqrt
    betamat    <- array(0, c(n, d + 1, d))
    Bmat       <- array(0, c(n, d + 1, d, d)) # matrix containing bias terms
    orderings  <- matrix(0, n, d)
    distances  <- matrix(0, n, d + 1)
    deltas     <- matrix(0, n, d)
    for (i in 1:n) {
      z        <- Z[i, ] 
      indNA    <- which(naMask[i, ] == 1)
      z[indNA] <- mu[indNA]
      Z[i, ]   <- z # needed to calculate the imputed values later on
      response <- Sigmaisqrt %*% (z - mu)
      weights  <- huberweights(x = z - mu, b = 1.5)
      
      larOut <- findCellPath_cpp(predictors = predictors,
                                 response = response,
                                 weights = weights,
                                 Sigmai = Sigmai,
                                 naMask = naMask[i, ])
      deltas[i, larOut$ordering] <- abs(diff(larOut$RSS))
      deltas[i, larOut$ordering] <- rev(cummax(rev(deltas[i, larOut$ordering])))
      
      betamat[i, , ] <- larOut$beta
      Bmat[i, , , ]  <- larOut$biasMat
      distances[i, ] <- larOut$RSS
      orderings[i, ] <- larOut$ordering
    }
    
    
    # 2A stage 2: sort and iterate through p-values
    # to determine the actual flagged cells
    # Goal of this stage is to take maxCol into account
    
    tiebraker    <- t(apply(orderings, 1, function(y) order(y))) + 1:n * d
    deltas_order <- order(deltas, tiebraker, # order of deltas in decreasin order
                          decreasing = TRUE) # second argument of order() solves ties
    
    # Fore imputation of NAs first (unique() keeps order of other cells):
    deltas_order <- unique(c(which(naMask == 1), deltas_order))
    
    NBimps_col   <- rep(0, d)
    cutpoints    <- rep(1, n) # where to stop the paths (1 = no imputes)
    droppedPaths <- rep(0, n) # which paths are dropped (="locked")
    W            <- matrix(0, n , d)
    for (i in 1:length(deltas_order)) {
      idx   <- deltas_order[i]
      delta <- deltas[idx]
      rownb <- (idx - 1) %% n + 1
      if (delta > qchisq(quant, 1)) {
        if (!droppedPaths[rownb]) {# check if case still has an open path
          colnb <- ((idx - 1) %/% n) + 1
          if (NBimps_col[colnb] < M) {
            cutpoints[rownb]  <- cutpoints[rownb] + 1
            NBimps_col[colnb] <- NBimps_col[colnb] + 1
            W[rownb, colnb]   <- 1
          } else {
            droppedPaths[rownb] <- 1
          }
        }
      } else {
        droppedPaths[rownb] <- 1
      }
    }
    
    #  Step 2B: impute cells
    finalBetas     <- matrix(0, n, d)
    finalBias      <- matrix(0, d, d)
    finalDistances <- rep(0, n)
    finalNbimps    <- cutpoints - 1
    finalW         <- matrix(0, n, d)
    for (i in 1:n) {
      finalBetas[i, ]   <- betamat[i, cutpoints[i], ]
      finalBias         <- finalBias + Bmat[i, cutpoints[i], , ]
      finalDistances[i] <- distances[i, cutpoints[i]]
      finalW[i, ]       <- (abs(betamat[i, cutpoints[i], ]) > 1e-10) + 0
    }
    
    Zimp <- Z - finalBetas
    
    # Step 2C: re-estimate the covariance matrix
    muSigmaNew <- mS_cov(Zimp, finalDistances, finalNbimps) 
    mu         <- muSigmaNew$mu
    Sigma      <- muSigmaNew$cov + finalBias / n # add bias matrix
    
    # Step 2D: bookkeeping and setting up for next iteration
    
    mus[nbits + 2, ]      <- mu
    Sigmas[nbits + 2, , ] <- Sigma
    invOut     <- mpinv(Sigma)
    Sigmai     <- invOut$Inv
    Sigmaisqrt <- invOut$InvSqrt
    
    convcrit <- sum((Sigmas[nbits + 2, , ] - Sigmas[nbits + 1, , ])^2) + 
      sum((mus[nbits + 2, ] - mus[nbits + 1, ])^2)
    nbits <- nbits + 1
  } # end of Step2: iteration
  
  # unstandardize and clean:
  Sigmas <- Sigmas[1:(nbits + 1), , ]
  
  Ximp   <- scale(Zimp, FALSE, 1 / locScale_init$scale)
  Ximp   <- scale(Ximp, -locScale_init$loc, FALSE)
  
  center_out  <- locScale_init$loc + mu * locScale_init$scale 
  cov_out     <- t(t(Sigma) * locScale_init$scale) * locScale_init$scale
  # icenter_out <- locScale_init$loc
  # icov_out    <- t(t(cov_init) * locScale_init$scale) * locScale_init$scale
  
  Sigmas      <- apply(Sigmas, 1,
                       function(y) t(t(y) * locScale_init$scale) * locScale_init$scale)
  dim(Sigmas) <- c(d, d, nbits + 1)
  Sigmas      <- aperm(Sigmas, perm = c(3, 1, 2))
  
  
  
  # Do one final run of the cellHandler algorithm
  CH.out <- cellHandler(X, center_out, cov_out, quant = quant)
  
  return(list(center = center_out,
              cov = cov_out,
              nits = nbits,
              Ximp = CH.out$Ximp,
              indcells = CH.out$indcells,
              indNAs = CH.out$indNAs,
              Zres = CH.out$Zres,
              cellPaths = CH.out$cellPaths,
              Zres_denom = CH.out$Zres_denom,
              checkDataSet_out = CD_out))
}
