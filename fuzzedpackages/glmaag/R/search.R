setClass(
  Class = 'search',
  slots = c(
    input = 'matrix',
    lambda1 = 'vector',
    lambda2 = 'vector',
    ns = 'matrix',
    coefs = 'array',
    intercept = 'matrix',
    loglik = 'matrix',
    fam = 'character'
  )
)


##' @title Fit glmaag model
##' @description Fit the glmaag model with given tuning parameters without cross validation or stability selection
##' @param y outcome
##' @param x predictors matrix
##' @param L Laplacian matrix for the network
##' @param tune whether to tune with an estimated network, default to be FALSE
##' @param est whether to estimate a network from the data
##' @param adaptl1 whether to adapt the l1 penalty, default to be TRUE
##' @param adaptl2 whether to adapt the sign for quadratic penalty, default to be TRUE
##' @param lam1 The tuning parameters for L1 penalty. If not defined, searched by default
##' @param lam2 The tuning parameters for quadratic penalty. If not defined, searched by default
##' @param dfmax maximum number of parameters allowed in the model, default to be p/2
##' @param w0 Weights for l1 penalty. If not defined, estimated via quadratic penalyzed regression
##' @param nfolds number of folds used in cross validation to obatin network sign estimate and l1 weight estimate, default to be five
##' @param pind indicator vector whether to put L1 penalty on the feature, 1 means penalyzed while 0 means not penalyzed, default to be all ones (all penalyzed)
##' @param intercept whether to include intercept. Ignore for Cox regression
##' @param standardize whether to standardize predictors
##' @param gam the parameter for l1 adaptive weight, default to be ones
##' @param fam family for the outcome, can be "Gaussian", "Logistic", and "Cox"
##' @param type1se whether to use one standard error or maximum rule when estimate network sign and l1 weight, default to be one standard error rule
##' @param measdev Whether to use deviance to tune when estimate l1 weight and network sign, default to be deviance. If not, use mean absolue error, area under ROC curve, or concordance index for Gaussian, Logistic, and Cox
##' @param maxiter maximum number of iterations, default to be 500
##' @param cri stoppint criterion, default to be 0.001
##' @param parallel whether to do parallel computing at each lambda2, need to set up parallel first, default to be FALSE
##' @return \item{input}{input predictors}
##' @return \item{lambda1}{l1 penalty parameter search sequence}
##' @return \item{lambda2}{quadratic penalty parameter search sequence}
##' @return \item{ns}{number of parameters selected given provided tuning parameter}
##' @return \item{coefs}{coefficients estimated}
##' @return \item{intercept}{intercepts estimated}
##' @return \item{loglik}{log likelihood estimated}
##' @return \item{fam}{family of the outcome}
##' @useDynLib glmaag, .registration = TRUE
##' @importFrom Rcpp evalCpp
##' @importFrom foreach foreach %dopar%
##' @importFrom stats coef lm glm residuals resid binomial
##' @importFrom methods new
##' @examples data(sampledata)
##' data(L0)
##' y <- sampledata$Y_Gau
##' x <- sampledata[, -(1:3)]
##' mod <- glmaag(y, x, L0)
##' @export
glmaag <- function(y, x, L, tune = F, est = T, gam = 1, lam1, lam2, nfolds = 5, dfmax, w0 , adaptl1 = T, adaptl2 = T, pind, intercept = T, standardize = T, maxiter = 10000, cri = .001, fam = 'Gaussian', measdev = T, type1se = T, parallel = F) {
  p <- ncol(x)
  n <- nrow(x)
  x <- as.matrix(x)
  input <- as.matrix(x)
  if (fam == 'Gaussian') {
    if (intercept) {
      ycen <- scale(y, scale = F)
      my <- mean(y)
    } else {
      ycen <- y
      my <- 0
    }
    cvwhich <- sample(rep(0:(nfolds - 1), length.out = n))
  } else if (fam == 'Logistic') {
    if (sum(y == 1) < nfolds | sum(y == 0) < nfolds) {
      stop('Too unbalance!')
    }
    pos <- which(y == 1)
    neg <- which(y == 0)
    poswhich <- sample(rep_len(0:(nfolds - 1), length(pos)))
    negwhich <- sample(rep_len(0:(nfolds - 1), length(neg)))
    cvwhich <- c()
    cvwhich[pos] <- poswhich
    cvwhich[neg] <- negwhich
    my <- mean(y)
    if(intercept) {
      b0 <- log(my / (1 - my))
    } else {
      b0 <- 0
    }
  } else {
    if(nnzero(y[, 2]) < nfolds) {
      stop('Too many censoring!')
    }
    if (intercept) {
      message('No intercept for Cox')
      intercept <- F
    }
    yord <- order(y)
    y <- y[yord]
    x <- x[yord, ]
    if (which(y[, 2] == 1)[1] > 1) {
      y <- y[-(1:(which(y[, 2] == 1)[1] - 1))]
      x <- x[-(1:(which(y[, 2] == 1)[1] - 1)), ]
    }
    tie <- any(duplicated(y[, 1]))
    if(tie) {
      tie1 <- which(y[, 1] %in% y[, 1][duplicated(y[, 1])])
      tie2 <- rank(y[, 1], ties.method = 'min')[tie1] - 1
      tie3 <- rank(y[, 1], ties.method = 'max')[tie1] - 1
      tie1 <- tie1 - 1
      ntie <- length(tie1)
    } else {
      tie1 <- tie2 <- tie3 <- 0
      ntie <- 0
    }
    obs <- which(y[, 2] == 1)
    cen <- which(y[, 2] == 0)
    if (length(cen) < nfolds) {
      cvwhich <- sample(rep(0:(nfolds - 1), length.out = n))
    } else {
      obswhich <- sample(rep(0:(nfolds - 1), length.out = length(obs)))
      cenwhich <- sample(rep(0:(nfolds - 1), length.out = length(cen)))
      cvwhich <- c()
      cvwhich[obs] <- obswhich
      cvwhich[cen] <- cenwhich
      for (i in 0:(nfolds - 1)) {
        if (which(cumsum(y[cvwhich == i, 2]) == 1)[1] > 1) {
          cvwhich[which(cvwhich == i)[1:(which(cumsum(y[cvwhich == i, 2]) == 1)[1] - 1)]] <- nfolds + 1
        }
      }
    }
  }
  
  
  if (missing(dfmax)) {
    dfmax <- round(.5*p)
    minlam10 <- F
  } else if (dfmax >= p) {
    minlam10 <- T
  } else {
    minlam10 <- F
  }
  
  
  if(standardize & intercept) {
    meanx <- colMeans(x)
    sdx <- apply(x,2,sd)
    xstd <- scale(x)
  } else if(standardize & !intercept) {
    meanx <- rep(0, p)
    sdx <- apply(x, 2, sd)
    xstd <- sweep(x, 2, sdx, '/')
  } else {
    meanx <- rep(0, p)
    sdx <- rep(1, p)
    xstd <- sweep(x, 2, sdx, '/')
  }
  
  
  
  names(xstd) <- NULL
  
  
  
  if(missing(pind)) {
    pind <- rep(1, p)
  }
  
  
  if (missing(lam1)) {
    if (fam == 'Gaussian'){
      if (all(pind == 1)) {
        lammax <- max(abs(crossprod(xstd, ycen))) / n
      } else {
        lammax <- max(abs(residuals(lm(ycen ~ xstd[, -which(pind == 0)])) %*% xstd[, -which(pind == 0)])) / n
      }
      
    } else if (fam == 'Logistic'){
      if (all(pind == 1)) {
        lammax <- max(abs(crossprod(xstd, y - my))) / n
      } else {
        lammax <- max(abs(residuals(glm(y ~ xstd[, -which(pind == 0)], family = binomial()), 'response') %*% xstd[, -which(pind == 0)])) / n
      }
    } else {
      if (all(pind == 1)) {
        lammax <- max(abs(resid(coxph(y ~ 1)) %*% xstd)) / n
      } else {
        lammax <- max(abs(resid(coxph(y ~ xstd[, -which(pind == 0)])) %*% xstd)) / n
      }
    }
    if (minlam10) {
      if (n < p) {
        lam1 <- c(exp(seq(log(lammax), log(.01*lammax), len = 100)), 0)
      } else {
        lam1 <- c(exp(seq(log(lammax), log(.0001*lammax), len = 100)), 0)
      }
    } else {
      if (n < p) {
        lam1 <- exp(seq(log(lammax), log(.01*lammax), len = 100))
      } else {
        lam1 <- exp(seq(log(lammax), log(.0001*lammax), len = 100))
      }
    }
    
  }
  if (missing(lam2)) {
    lam2 <- c(.01*2^(7:0), 0)
  }
  
  if(n < p) {
    lam20 <- lam2[lam2 != 0]
  } else {
    lam20 <- lam2
  }
  
  if(missing(L)) {
    if (tune) {
      message('Cannot tune without an input network!')
    }
    gototune <- F
    if (est) {
      message('Use estimated network')
      L <- getS(xstd)
      dl <- rep(1, p)
      Lad <- L
      diag(Lad) <- 0
    } else {
      message('Use elastic net')
      dl <- rep(1, p)
      Lad <- matrix(0, p, p)
      L <- diag(p)
    }
  } else if(any(eigen(L)$valu < .001)){
    dl <- rep(1, p)
    Lad <- matrix(0, p, p)
    L <- diag(p)
    warning('Networks is not positive definite, change to identity matrix!')
    if (tune) {
      if (est) {
        L2 <- getS(xstd)
        gototune <- T
      } else {
        gototune <- F
      }
    } else {
      gototune <- F
    }
  } else {
    dl <- diag(L)
    Lad <- L
    diag(Lad) <- 0
    gototune <- tune
    if (tune) {
      if (est) {
        L2 <- getS(xstd)
      } else {
        L2 <- diag(p)
      }
    }
  }
  
  tuneweight <- c(1, 0)
  
  if(minlam10 | adaptl1) {
    if (parallel) {
      cvmod0 <- foreach(i = 0:(nfolds - 1), .combine = 'cbind') %dopar% {
        if (fam == 'Gaussian') {
          cvlrnet2_pal(xstd[cvwhich != i, ], xstd[cvwhich == i, ], ycen[cvwhich != i], ycen[cvwhich == i], L, lam20, measdev)
        } else if (fam == 'Logistic') {
          cvloginet2_pal(xstd[cvwhich != i, ], xstd[cvwhich == i, ], y[cvwhich != i], y[cvwhich == i], rep(0, p), L, lam20, intercept, measdev, maxiter, cri)
        } else {
          ytr <- y[cvwhich != i]
          tietr <- any(duplicated(ytr[, 1]))
          if(tietr) {
            tie1tr <- which(ytr[, 1] %in% ytr[, 1][duplicated(ytr[, 1])])
            tie2tr <- rank(ytr[, 1], ties.method = 'min')[tie1tr] - 1
            tie3tr <- rank(ytr[, 1], ties.method = 'max')[tie1tr] - 1
            tie1tr <- tie1tr - 1
            ntietr <- length(tie1tr)
          } else {
            tie1tr <- tie2tr <- tie3tr <- 0
            ntietr <- 0
          }
          cvCoxnet2_pal(rep(0, p), xstd, xstd[cvwhich != i, ], xstd[cvwhich == i, ], y[cvwhich == i, 1], tie, ntie, tie1, tie2, tie3, tietr, ntietr, tie1tr, tie2tr, tie3tr, y[, 2], y[cvwhich != i, 2], y[cvwhich == i, 2], L, lam20, measdev, maxiter, cri)
        }
      }
      
    } else {
      cvmod0 <- switch (fam,
        Gaussian = cvlrnet2(nfolds, xstd, ycen, L, lam20, cvwhich, measdev),
        Logistic = cvloginet2(nfolds, xstd, y, rep(0, p), L, lam20, intercept, cvwhich, measdev, maxiter, cri),
        Cox = cvCoxnet2(nfolds, rep(0, p), xstd, y[, 1], tie, ntie, tie1, tie2, tie3, y[, 2], L, lam20, cvwhich, measdev, maxiter, cri)
      )
    }
    
    cvmean0 <- rowMeans(cvmod0, na.rm = T)
    cvse0 <- apply(cvmod0, 1, sd, na.rm = T)/sqrt(nfolds)
    
    
    if (type1se) {
      b <- as.vector(switch (fam,
        Gaussian = lrnet(xstd, ycen, lam20[which(cvmean0 >= max(cvmean0, na.rm = T) - cvse0[which.max(cvmean0)])[1]]*L),
        Logistic = loginet(b0, rep(0, p), xstd, y, lam20[which(cvmean0 >= max(cvmean0, na.rm = T) - cvse0[which.max(cvmean0)])[1]]*L, intercept, maxiter, cri)[-1],
        Cox = Coxnet(rep(0, p), xstd, tie, tie1, tie2, tie3, y[, 2], lam20[which(cvmean0 >= max(cvmean0, na.rm = T) - cvse0[which.max(cvmean0)])[1]]*L, maxiter, cri)
      ))
    } else {
      b <- as.vector(switch (fam,
        Gaussian = lrnet(xstd, ycen, lam20[which.max(cvmean0)]*L),
        Logistic = loginet(b0, rep(0, p), xstd, y, lam20[which.max(cvmean0)]*L, intercept, maxiter, cri)[-1],
        Cox = Coxnet(rep(0, p), xstd, tie, tie1, tie2, tie3, y[, 2], lam20[which.max(cvmean0)]*L, maxiter, cri)
      ))
    }
  }
  
  
  
  if (adaptl1) {
    w <- (abs(b)^(-gam))*pind
  } else {
    w <- rep(1, p)*pind
  }
  
  if (!missing(w0)) {
    w <- w0
  }
  lam1 <- lam1/min(w)
  lam10 <- lam1[lam1 != 0]
  if (parallel) {
    outlist <- foreach(i = lam2) %dopar% {
      switch (fam,
        Gaussian = lraagg_search_pal(dfmax, rep(0, p), w, xstd, ycen, Lad, dl, lam1, i, maxiter, cri),
        Logistic = logiaagg_search_pal(dfmax, rep(0, p), xstd, y, Lad, dl, w, lam1, i, intercept, maxiter, cri),
        Cox = Coxaagg_search_pal(ntie, dfmax, rep(0, p), w, xstd, tie, tie1, tie2, tie3, y[, 2], Lad, dl, lam1, i, maxiter, cri)
      )
    }
    coefs <- simplify2array(lapply(outlist, '[[', 'coef'))
    ns <- sapply(outlist, '[[', 'n')
    loglik <- sapply(outlist, '[[', 'loglik')
    intermat <- sapply(outlist, '[[', 'interc')
  } else {
    mod <- switch (fam,
      Gaussian = lraagg_search(dfmax, rep(0, p), w, xstd, ycen, Lad, dl, lam1, lam2, maxiter, cri),
      Logistic = logiaagg_search(dfmax, rep(0, p), xstd, y, Lad, dl, w, lam1, lam2, intercept, maxiter, cri),
      Cox = Coxaagg_search(ntie, dfmax, rep(0, p), w, xstd, tie, tie1, tie2, tie3, y[, 2], Lad, dl, lam1, lam2, maxiter, cri)
    )
    ns <- mod$n
    coefs <- mod$coef
    loglik <- mod$loglik
    intermat <- mod$interc
  }
  
  
  nacoef <- which(ns > dfmax, arr.ind = T)
  if (nrow(nacoef) > 0) {
    for (i in 1:nrow(nacoef)) {
      coefs[, nacoef[i, 1], nacoef[i, 2]] <- NA
    }
    ns[nacoef] <- NA
    loglik[nacoef] <- NA
  }
  
  coefs <- sweep(coefs, c(1, 3), sdx, '/')
  abb <- apply(coefs, 2:3, function(x) sum(meanx*x))
  intermat <- switch (fam,
    Gaussian = matrix(my, length(lam1), length(lam2)) - abb,
    Logistic = intermat - abb,
    Cox = matrix(0, length(lam1), length(lam2))
  )
  
  intermat[nacoef] <- NA
  if (min(lam1) == 0 & min(lam2) == 0 & n <= p) {
    intermat[length(lam1), length(lam2)] <- NA
  }

  
  rownames(ns) <- rownames(intermat) <- rownames(loglik) <- lam1
  colnames(ns) <- colnames(intermat) <- colnames(loglik) <- lam2
  if (dim(coefs)[2] == 1) {
    dimnames(coefs)[[2]] <- list(lam1)
  } else {
    dimnames(coefs)[[2]] <- lam1
  }
  
  if (dim(coefs)[3] == 1) {
    dimnames(coefs)[[3]] <- list(lam2)
  } else {
    dimnames(coefs)[[3]] <- lam2
  }
  
  
  res <- new('search', input = input, lambda1 = lam1, lambda2 = lam2, ns = ns, coefs = coefs, intercept = intermat, loglik = loglik, fam = fam)
  structure(res, class = 'glmaag')
}


##' @title Coefficients for glmaag
##' @description Get coefficients for glmaag objects
##' @param object fitted glmaag object
##' @param lam1 lambda1 sequence need coefficients, must be within the fitted model
##' @param lam2 lambda2 sequence need coefficients, must be within the fitted model
##' @param ... \dots
##' @return coefficients
##' @method coef glmaag
##' @examples data(sampledata)
##' data(L0)
##' y <- sampledata$Y_Gau
##' x <- sampledata[, -(1:3)]
##' mod <- glmaag(y, x, L0)
##' cc <- coef(mod)
##' @export
coef.glmaag <- function(object, lam1, lam2, ...) {
  if (missing(lam1) & missing(lam2)) {
    return(object@coefs)
  } else if (missing(lam1)) {
    if (all(lam2 %in% object@lambda2)) {
      return(object@coefs[, as.character(lam2)])
    } else {
      stop('Tuning parameter not in search values!')
    }
  } else if (missing(lam2)) {
    if (all(lam1 %in% object@lambda1)) {
      return(object@coefs[as.character(lam1), ])
    } else {
      stop('Tuning parameter not in search values!')
    }
  } else {
    if(all(lam2 %in% object@lambda2) & all(lam1 %in% object@lambda1)) {
      return(object@coefs[as.character(lam1), as.character(lam2)])
    } else {
      stop('Tuning parameter not in search values!')
    }
  }
}

##' @title Prediction for glmaag
##' @description Prediction using glmaag model
##' @param object fitted glmaag object
##' @param x The new dataset to be predicted, do training prediction if x is missing
##' @param lam1 lambda1 sequence for prediction, must be within the fitted model
##' @param lam2 lambda2 sequence for prediction, must be within the fitted model
##' @param type type of prediction (can be "link", "reponse"), ignored for Gaussian model. "link" is the linear predicted score, "response" is the predicted probability for logistic model and relative risk for Cox model
##' @param cutp the cut off value for binary outcome, default to be 0.5
##' @param ... \dots
##' @return predicted values
##' @method predict glmaag
##' @examples data(sampledata)
##' data(L0)
##' y <- sampledata$Y_Gau
##' x <- sampledata[, -(1:3)]
##' mod <- glmaag(y, x, L0)
##' pp <- predict(mod)
##' @export
predict.glmaag <- function(object, x, lam1, lam2, type = 'link', cutp = .5, ...) {
  if (missing(x)) {
    x <- object@input
  }
  x <- as.matrix(x)
  if (missing(lam1) & missing(lam2)) {
    ypre <- array(dim = c(nrow(x), dim(object@coefs)[-1]))
    dimnames(ypre) <- dimnames(object@coefs)
    for (i in as.character(object@lambda1)) {
      for (j in as.character(object@lambda2)) {
        ypre[, i, j] <- object@intercept[i, j] + x %*% object@coefs[, i, j]
      }
    }
  } else if (missing(lam1)) {
    if (all(lam2 %in% object@lambda2)) {
      ypre <- array(dim = c(nrow(x), dim(object@coefs)[2], length(lam2)))
      dimnames(ypre)[2] <- dimnames(object@coefs)[2]
      dimnames(ypre)[3] <- as.character(lam2)
      for (i in as.character(object@lambda1)) {
        for (j in as.character(lam2)) {
          ypre[, i, j] <- object@intercept[i, j] + x %*% object@coefs[, i, j]
        }
      }
    } else {
      stop('Tuning parameter not in search values!')
    }
  } else if (missing(lam2)) {
    if (all(lam1 %in% object@lambda1)) {
      ypre <- array(dim = c(nrow(x), length(lam1), dim(object@coefs)[3]))
      dimnames(ypre)[2] <- as.character(lam1)
      dimnames(ypre)[3] <- dimnames(object@coefs)[3]
      for (i in as.character(lam1)) {
        for (j in as.character(object@lambda2)) {
          ypre[, i, j] <- object@intercept[i, j] + x %*% object@coefs[, i, j]
        }
      }
    } else {
      stop('Tuning parameter not in search values!')
    }
  } else {
    if(all(lam2 %in% object@lambda2) & all(lam1 %in% object@lambda1)) {
      ypre <- array(dim = c(nrow(x), length(lam1), length(lam2)))
      dimnames(ypre)[2] <- as.character(lam1)
      dimnames(ypre)[3] <- as.character(lam2)
      for (i in as.character(lam1)) {
        for (j in as.character(lam2)) {
          ypre[, i, j] <- object@intercept[i, j] + x %*% object@coefs[, i, j]
        }
      }
    } else {
      stop('Tuning parameter not in search values!')
    }
  }
  
  if(type == 'link' | object@fam == 'Gaussian') {
    return(ypre)
  } else {
    if (object@fam == 'Logistic') {
      if (type == 'response') {
        return(1 / (1 + exp(-ypre)))
      } else {
        return(as.numeric(1 / (1 + exp(-ypre)) > cutp))
      }
    } else {
      return(exp(ypre))
    }
  }
}

##' @importFrom  ggplot2 ggplot aes_ geom_line theme_bw theme xlab ggtitle
plotcoef <- function(x, datplot) {
  lambda2 <- NULL
  ggplot(datplot[lambda2 == x], aes_(x =~ lambda1, y =~ Coefficient, color =~ variable)) + geom_line() + theme_bw() + xlab(expression(lambda[1])) + ggtitle(substitute(lambda[2] == p, list(p = x))) + theme(plot.title = element_text(hjust = .5), legend.position = 'none', axis.title.y = element_blank())
}



##' @title Paths for glmaag object
##' @description Generates coefficients, lok likelihood, or number of parameters paths for glmaag models
##' @param x glmaag object
##' @param col_count number of columns shown in the plot (when type = 'coef)
##' @param type can be "coef" (coefficients paths), "loglik" (log likelihood paths), or "n" (number of parameters paths)
##' @param ... \dots
##' @importFrom data.table as.data.table melt setnames := rbindlist
##' @importFrom ggplot2 ggplot geom_line labs theme_bw theme aes_
##' @importFrom gridExtra grid.arrange
##' @return plots
##' @method plot glmaag
##' @examples data(sampledata)
##' data(L0)
##' y <- sampledata$Y_Gau
##' x <- sampledata[, -(1:3)]
##' mod <- glmaag(y, x, L0)
##' gg <- plot(mod, type = 'loglik')
##' @export
plot.glmaag <- function(x, col_count = 3, type = 'coef', ...) {
  lambda1 <- lambda2 <- NULL
  if (type == 'coef') {
    coefcube <- x@coefs
    coeftolist <- lapply(seq(dim(coefcube)[3]), function(x) coefcube[ , , x])
    names(coeftolist) <- dimnames(coefcube)[3][[1]]
    plottolist <- list(); length(plottolist) <- length(coeftolist)
    for(i in 1:length(coeftolist)) {
      plottolist[[i]] <- as.data.table(t(coeftolist[[i]]))
      plottolist[[i]][, lambda1 := as.numeric(colnames(coeftolist[[i]]))]
      plottolist[[i]][, lambda2 := names(coeftolist)[i]]
    }
    datplot <- melt(rbindlist(plottolist), id = c('lambda1', 'lambda2'))
    setnames(datplot, 'value', 'Coefficient')
    plotlist <- lapply(x@lambda2, plotcoef, datplot = datplot)
    grid.arrange(grobs = plotlist, ncol = col_count, left = 'Coefficients')
  } else if (type == 'loglik') {
    loglikmat <- as.data.table(x@loglik)
    setnames(loglikmat, as.character(x@lambda2))
    loglikmat[, lambda1 := x@lambda1]
    loglikdat <- melt(loglikmat, 'lambda1')
    setnames(loglikdat, c('lambda1', 'lambda2', 'loglik'))
    ggplot(loglikdat, aes_(x =~ lambda1, y =~ loglik, color =~ lambda2)) + geom_line() + labs(title = 'Log Likelihood', color = expression(lambda[2]~''), x = expression(lambda[1])) + theme_bw() + theme(axis.title.y = element_blank(), plot.title = element_text(hjust = .5), legend.position = 'bottom')
  } else {
    nmat <- as.data.table(x@ns)
    setnames(nmat, as.character(x@lambda2))
    nmat[, lambda1 := x@lambda1]
    ndat <- melt(nmat, 'lambda1')
    setnames(ndat, c('lambda1', 'lambda2', 'n'))
    ggplot(ndat, aes_(x =~ lambda1, y =~ n, color =~ lambda2)) + geom_line() + labs(title = 'Number of Parameters', color = expression(lambda[2]~''), x = expression(lambda[1])) + theme_bw() + theme(axis.title.y = element_blank(), plot.title = element_text(hjust = .5), legend.position = 'bottom')
  }
}
