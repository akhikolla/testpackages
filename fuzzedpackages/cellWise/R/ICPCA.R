ICPCA <- function(X, k, scale = FALSE, maxiter = 20, tol = 0.005,
                 tolProb = 0.99, distprob = 0.99) {
  #  
  #   This function is based on a Matlab function from
  #   Missing Data Imputation Toolbox v1.0
  #   A. Folch-Fortuny, F. Arteaga and A. Ferrer
  #   
  # Its inputs are:
  #
  # X : the input data, which must be a matrix or a data frame.
  #     It may contain NA's. It must always be provided.
  # k : the desired number of principal components.
  # scale     : a value indicating whether and how the original
  #             variables should be scaled. If scale=FALSE (default)
  #             or scale=NULL no scaling is performed (and a vector
  #             of 1s is returned in the $scaleX slot).
  #             If scale=TRUE the data are scaled to have a
  #             standard deviation of 1.
  #             Alternatively scale can be a function like mad,
  #             or a vector of length equal to the number of columns
  #             of x.
  #             The resulting scale estimates are returned in the
  #             $scaleX slot of the ICPCA output.
  # maxiter   : maximum number of iterations. Default is 20.
  # tol       : tolerance for iterations. Default is 0.005.
  # tolProb   : tolerance probability for residuals. Defaults to 0.99.
  # distprob  : probability determining the cutoff values for
  #             orthogonal and score distances. Default is 0.99.
  #
  # The outputs are:
  #
  # scaleX      : the scales of the columns of X.
  # k           : the number of principal components.
  # loadings    : the columns are the k loading vectors.
  # eigenvalues : the k eigenvalues.
  # center      : vector with the fitted center.
  # covmatrix   : estimated covariance matrix.
  # n.obs       : number of cases.
  # It          : number of iteration steps.
  # diff        : convergence criterion.
  # X.NAimp     : data with all NA's imputed by MacroPCA.
  # scores      : scores of X.NAimp
  # OD          : orthogonal distances of the rows of X.NAimp
  # cutoffOD    : cutoff value for the OD.
  # SD          : score distances of the rows of X.NAimp
  # cutoffSD    : cutoff value for the SD.
  # indrows     : row numbers of rowwise outliers.
  # residScale  : scale of the residuals.
  # stdResid    : standardized residuals. Note that these are NA
  #               for all missing values of the data X.
  # indcells    : indices of cellwise outliers.
  
  X <- as.matrix(X)
  n <- dim(X)[1]
  p <- dim(X)[2]
  mis <- is.na(X)
  
  # Determine scale
  if (is.logical(scale)) {
    if (!scale) scaleX <- vector("numeric", p) + 1
    if (scale) scaleX <- apply(X, 2, FUN = sd, na.rm = TRUE)
    scale <- sd # scale is set to the function sd
  } else if (is.function(scale)) {
    scaleX <- apply(X, 2, FUN = scale, na.rm = TRUE)
  } else {
    scaleX <- scale
    scale <- sd # scale is set to the function sd
  }
  
  X <- sweep(X, 2, scaleX, "/") # standardization
  XO <- X # still has the NA's
  rm(X)  # to save space
  Xnai <- XO  # initialize the NA-imputed Xnai
  
  diff <- 100
  It <- 0
  if (any(mis)) { # if there are missings
    misrc <- which(is.na(XO), arr.ind = TRUE)
    r <- misrc[, 1]
    c <- misrc[, 2]
    meanc <- colMeans(XO, na.rm = TRUE) # Mean vector
    for (i in 1:length(r)) {
      Xnai[r[i], c[i]] <- meanc[c[i]];  
      # Impute missing data with mean values 
    }
    while (It < maxiter & diff > tol) { # Iterate
      It <- It + 1;
      # Xmis <- Xnai[mis]        # current imputations
      mXnai <- colMeans(Xnai)       # mean vector
      Xnaic <- sweep(Xnai,2,mXnai)  # centered Xnai
      if (n < p) {
        XnaicSVD <- svd(t(Xnaic))
        Pr <- as.matrix(XnaicSVD$u[, 1:k]) 
        # reduced loadings matrix
      } else {
        XnaicSVD <- svd(Xnaic)
        Pr <- as.matrix(XnaicSVD$v[, 1:k])
        # reduced loadings matrix
      }
      Tr <- Xnaic %*% Pr          # reduced scores matrix     
      Xnaihat <- Tr %*% t(Pr)     # fit to Xnaic
      Xnaihat <- sweep(Xnaihat, 2, mXnai, "+") # fit to Xnai
      Xnai[mis] <- Xnaihat[mis]         # impute missings
      # d <- (Xnai[mis]-Xmis)^2       
      if (It > 1) diff <- maxAngle(Pr,PrPrev)
      PrPrev <- Pr
    } 
  } else {# if there are no missings
    mXnai <- colMeans(Xnai)      # mean vector
    Xnaic <- sweep(Xnai, 2, mXnai) # centered Xnai
    if (n < p) {
      XnaicSVD <- svd(t(Xnaic))
      Pr <- as.matrix(XnaicSVD$u[, 1:k]) 
      # reduced loadings matrix
    } else {
      XnaicSVD <- svd(Xnaic)
      Pr <- as.matrix(XnaicSVD$v[, 1:k]) 
      # reduced loadings matrix
    }
    Tr <- Xnaic %*% Pr        # reduced scores matrix     
    Xnaihat <- Tr %*% t(Pr)   # fit to Xnaic       
    Xnaihat <- sweep(Xnaihat, 2, mXnai, "+") # fit to Xnai
  }
  
  rank   <- rankMM(Xnaic, sv = XnaicSVD$d)
  S      <- cov(Xnai) # is not the EM covariance
  center <- colMeans(Xnai)
  
  res <- list(loadings = Pr, eigenvalues = (XnaicSVD$d ^ 2) / (n - 1),
              center = center, k = k, h = n, alpha = 1,  scores = Tr)
  NAimp <- pca.distances.classic(res, Xnai, rank, distprob)
  NAimp$indrowsnai <- which(NAimp$OD > NAimp$cutoffOD)
  
  # Compute standardized residuals with NA's
  stdResid   <- XO - Xnaihat
  residScale <- apply(stdResid, 2, FUN = scale, na.rm = TRUE)
  stdResid   <- sweep(stdResid, 2, residScale, "/")
  indcells   <- which(abs(stdResid) > sqrt(qchisq(tolProb, 1)))
  
  Xnai = sweep(Xnai, 2, scaleX, "*") # unstandardize
  center = center * scaleX
  
  return(list(scaleX = scaleX,
              k = k,
              loadings = Pr,
              eigenvalues = res$eigenvalues[1:k],
              center = center,
              covmatrix = S,
              It = It,
              diff = diff,
              X.NAimp = Xnai,
              scores = Tr,
              OD = NAimp$OD,
              cutoffOD = NAimp$cutoffOD,
              SD = NAimp$SD,
              cutoffSD = NAimp$cutoffSD,              
              indrows = NAimp$indrowsnai,
              residScale = residScale,
              stdResid = stdResid,
              indcells = indcells))
} # ends ICPCA



pca.distances.classic <- function(obj, data, r, crit = 0.99) {   
  # based on rrcov:::.distances
  n <- nrow(data)
  q <- ncol(obj$scores)
  smat <- diag(obj$eigenvalues, ncol = q)
  nk <- min(q, rankMM(smat))
  if (nk < q) warning(paste("The smallest eigenvalue is ", 
                           obj$eigenvalues[q], 
                           " so the diagonal matrix of the eigenvalues",
                           "cannot be inverted!", sep = ""))  
  mySd <- sqrt(mahalanobis(as.matrix(obj$scores[, 1:nk]),
                           rep(0, nk), diag(obj$eigenvalues[1:nk], ncol = nk)))
  cutoffSD <- sqrt(qchisq(crit, obj$k))
  out <- list(SD = mySd, cutoffSD = cutoffSD)
  out$OD <- apply(data - matrix(rep(obj$center, times = n), nrow = n, 
                                byrow = TRUE) - obj$scores %*% t(obj$loadings),
                  1, vecnorm)
  if (is.list(dimnames(obj$scores))) {
    names(out$OD) <- dimnames(obj$scores)[[1]]
  }
  out$cutoffOD <- 0
  if (obj$k != r) {
    out$cutoffOD <- critOD(out$OD, crit = crit, classic = TRUE) }
  return(out)
}

