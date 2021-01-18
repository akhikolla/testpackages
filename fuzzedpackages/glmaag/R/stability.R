setClass(
  Class = 'ss',
  slots = c(
    input = 'matrix',
    lambda1 = 'vector',
    lambda2 = 'vector',
    lambda1_ss = 'numeric',
    lambda2_ss = 'numeric',
    n_ss = 'numeric',
    ssm = 'matrix',
    ssf = 'array',
    intercept_ss = 'numeric',
    coef_ss = 'vector',
    fam = 'character'
  )
)


findlam <- function(x, beta) {
  sum(x <= beta)
}



##' @title Stability selection for glmaag
##' @description Do stability selection for glmaag
##' @param y outcome
##' @param x predictors matrix
##' @param L Laplacian matrix for the first network
##' @param nfolds number of folds used in cross validation to estimate the l1 weights or network tuning, default to be five
##' @param subn number of samples in each subset, default to be n/2 if n<400 and 10sqrt(10) if n>400
##' @param nsam number of subsets, default to be 100
##' @param gam The power of weights of l1 penalty, default to be ones
##' @param beta the cut off for instability score
##' @param tune whether to tune the input network with estimated network or identity matrix, ignored if no input network
##' @param est when there is no input network whether to use estimated network or identiy matrix (elastic net) or mixed the network with estimated network or identity matrix, default to be estimated network
##' @param lam1 The tuning parameters for l1 penalty. If not defined, searched by default
##' @param lam2 The tuning parameters for quadratic penalty. If not defined, searched by default
##' @param w0 Weights for l1 penalty. If not defined, estimated via quadratic penalyzed regression
##' @param adaptl1 whether to adapt the l1 penalty, default to be TRUE
##' @param adaptl2 whether to adapt the sign for quadratic penalty, default to be TRUE
##' @param pind indicator vector whether to put l1 penalty on the feature, 1 means penalyzed while 0 means not penalyzed, default to be all ones (all penalyzed)
##' @param intercept whether to include intercept. Ignore for Cox regression
##' @param standardize whether to standardize predictors
##' @param fam family for the outcome, can be "Gaussian", "Logistic", and "Cox"
##' @param type1se whether to use one standard error or maximum rule for l1 weight estimation and network sign, default to be one standard error rule
##' @param measdev Whether to use deviance to tune when estimate l1 weight and network sign, default to be deviance. If not, use mean absolue error, area under ROC curve, or concordance index for Gaussian, Logistic, and Cox
##' @param maxiter maximum number of iterations, default to be 500
##' @param cri stoppint criterion, default to be 0.001
##' @param parallel whether to do parallel computing at each subset, need to set up parallel first, default to be FALSE
##' @return \item{input}{input matrix for predictors}
##' @return \item{lambda1}{searching sequence for l1 penalty parameters}
##' @return \item{lamdba2}{searching sequence for quadratic penalty parameters}
##' @return \item{lambda1_ss}{optimal l1 parameter}
##' @return \item{lambda2_ss}{optimal quadratic parameter}
##' @return \item{n_ss}{number of parameters obtained by the optimal model}
##' @return \item{ssm}{instability score paths}
##' @return \item{ssf}{selection probability paths}
##' @return \item{intercept_ss}{intercept estimated by the optimal model}
##' @return \item{coef_ss}{coefficients estimated by the optimal model}
##' @return \item{fam}{the family of the outcome}
##' @useDynLib glmaag, .registration = TRUE
##' @importFrom Rcpp evalCpp
##' @importFrom foreach foreach %dopar%
##' @importFrom stats coef residuals resid glm lm binomial
##' @references Meinshausen, N., & B{\"u}hlmann, P. (2010). Stability selection. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 72(4), 417-473.
##' @references Liu, H., Roeder, K., & Wasserman, L. (2010). Stability approach to regularization selection (stars) for high dimensional graphical models. In Advances in neural information processing systems (pp. 1432-1440).
##' @examples data(sampledata)
##' data(L0)
##' y <- sampledata$Y_Gau
##' x <- sampledata[, 4:6]
##' mod <- ss_glmaag(y, x, L0[seq_len(3), seq_len(3)], nsam = 3)
##' @export
ss_glmaag <- function(y, x, L, nfolds = 5, subn, nsam = 100, beta = .15, gam = 1, tune = F, est = T, lam1, lam2, w0 , adaptl1 = T, adaptl2 = T, pind, intercept = T, standardize = T, maxiter = 10000, cri = .001, fam = 'Gaussian', measdev = T, type1se = T, parallel = F) {
  p <- ncol(x)
  n <- nrow(x)
  x <- as.matrix(x)
  input <- x
  if (missing(subn)) {
    subn <- ifelse(n < 400, floor(.5*n), floor(10*sqrt(n)))
  }
  if (fam == 'Gaussian') {
    if (intercept) {
      ycen <- scale(y, scale = F)
      my <- mean(y)
    } else {
      ycen <- y
      my <- 0
    }
    cvwhich <- sample(rep(0:(nfolds - 1), length.out = n))
    sswhich <- replicate(nsam, sort(sample.int(n, subn)))
  } else if (fam == 'Logistic') {
    if (sum(y == 1) < nfolds | sum(y == 0) < nfolds) {
      stop('Too unbalance!')
    }
    pos <- which(y == 1)
    neg <- which(y == 0)
    poswhich <- sample(rep_len(0:(nfolds - 1), length(pos)))
    negwhich <- sample(rep_len(0:(nfolds - 1), length(neg)))
    posn <- round(mean(y) * subn)
    negn <- subn - posn
    cvwhich <- c()
    cvwhich[pos] <- poswhich
    cvwhich[neg] <- negwhich
    sswhich <- replicate(nsam, sort(c(sample(pos, posn), sample(neg, negn))))
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
      obswhich <- sample(rep_len(0:(nfolds - 1), length(obs)))
      cenwhich <- sample(rep_len(0:(nfolds - 1), length(cen)))
      cvwhich <- c()
      cvwhich[obs] <- obswhich
      cvwhich[cen] <- cenwhich
    }
    for (i in 0:(nfolds - 1)) {
      if (which(cumsum(y[cvwhich == i, 2]) == 1)[1] > 1) {
        cvwhich[which(cvwhich == i)[1:(which(cumsum(y[cvwhich == i, 2]) == 1)[1] - 1)]] <- nfolds + 1
      }
    }
    if (floor((1 - mean(y[, 2]))*subn) < 1) {
      sswhich <- replicate(nsam, sort(sample.int(n, subn)))
    } else {
      obsn <- round(mean(y[, 2]) * subn)
      cenn <- subn - obsn
      sswhich <- replicate(nsam, sort(c(sample(obs, obsn), sample(cen, cenn))))
    }
  }

  sswhich <- sswhich - 1
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
    if (n < p) {
      lam1 <- exp(seq(log(lammax), log(.01*lammax), len = 100))
    } else {
      lam1 <- exp(seq(log(lammax), log(.0001*lammax), len = 100))
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

  if (gototune) {
    modL <- switch (fam,
                    Gaussian = tune_network(ycen, xstd, L, L2, adaptl2 = adaptl2, nfolds = nfolds, intercept = F, standardize = F, fam = fam, type1se = type1se, measdev = measdev, maxiter = maxiter, cri = cri, parallel = parallel),
                    Logistic = tune_network(y, xstd, L, L2, adaptl2 = adaptl2, nfolds = nfolds, intercept = intercept, standardize = F, fam = fam, type1se = type1se, measdev = measdev, maxiter = maxiter, cri = cri, parallel = parallel),
                    Cox = tune_network(y, xstd, L, L2, adaptl2 = adaptl2, nfolds = nfolds, intercept = F, standardize = F, fam = fam, type1se = type1se, measdev = measdev, maxiter = maxiter, cri = cri, parallel = parallel)
    )
    L <- modL@est
    tuneweight <- modL@weight
    dl <- diag(L)
    Lad <- L
    diag(Lad) <- 0
  } else if(adaptl2) {
    II <- diag(p)
    if (parallel) {
      cvmod0 <- foreach(i = 0:(nfolds - 1), .combine = 'cbind') %dopar% {
        if (fam == 'Gaussian') {
          cvlrnet2_pal(xstd[cvwhich != i, ], xstd[cvwhich == i, ], ycen[cvwhich != i], ycen[cvwhich == i], II, lam20, measdev)
        } else if (fam == 'Logistic') {
          cvloginet2_pal(xstd[cvwhich != i, ], xstd[cvwhich == i, ], y[cvwhich != i], y[cvwhich == i], rep(0, p), II, lam20, intercept, measdev, maxiter, cri)
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
          cvCoxnet2_pal(rep(0, p), xstd, xstd[cvwhich != i, ], xstd[cvwhich == i, ], y[cvwhich == i, 1], tie, ntie, tie1, tie2, tie3, tietr, ntietr, tie1tr, tie2tr, tie3tr, y[, 2], y[cvwhich != i, 2], y[cvwhich == i, 2], II, lam20, measdev, maxiter, cri)
        }
      }

    } else {
      cvmod0 <- switch (fam,
                        Gaussian = cvlrnet2(nfolds, xstd, ycen, II, lam20, cvwhich, measdev),
                        Logistic = cvloginet2(nfolds, xstd, y, rep(0, p), II, lam20, intercept, cvwhich, measdev, maxiter, cri),
                        Cox = cvCoxnet2(nfolds, rep(0, p), xstd, y[, 1], tie, ntie, tie1, tie2, tie3, y[, 2], II, lam20, cvwhich, measdev, maxiter, cri)
      )
    }

    cvmean0 <- rowMeans(cvmod0, na.rm = T)

    if (type1se) {
      cvse0 <- apply(cvmod0, 1, sd, na.rm = T)/sqrt(nfolds)
      b <- as.vector(switch (fam,
                             Gaussian = lrnet(xstd, ycen, lam20[which(cvmean0 >= max(cvmean0, na.rm = T) - cvse0[which.max(cvmean0)])[1]]*II),
                             Logistic = loginet(0, rep(0, p), xstd, y, lam20[which(cvmean0 >= max(cvmean0, na.rm = T) - cvse0[which.max(cvmean0)])[1]]*II, intercept, maxiter, cri)[-1],
                             Cox = Coxnet(rep(0, p), xstd, tie, tie1, tie2, tie3, y[, 2], lam20[which(cvmean0 >= max(cvmean0, na.rm = T) - cvse0[which.max(cvmean0)])[1]]*II, maxiter, cri)
      ))
    } else {
      b <- as.vector(switch (fam,
                             Gaussian = lrnet(xstd, ycen, lam20[which.max(cvmean0)]*II),
                             Logistic = loginet(0, rep(0, p), xstd, y, lam20[which.max(cvmean0)]*II, intercept, maxiter, cri)[-1],
                             Cox = Coxnet(rep(0, p), xstd, tie, tie1, tie2, tie3, y[, 2], lam20[which.max(cvmean0)]*II, maxiter, cri)
      ))
    }

    ss <- diag(sign(b))
    Lad <- ss %*% -abs(Lad) %*% ss
    L <- ss %*% L %*% ss
  }

  if(adaptl1) {
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
        Cox = Coxnet(rep(0, p), xstd, y[, 1], y[, 2], lam20[which.max(cvmean0)]*L, maxiter, cri)
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

  if (parallel) {
    outlist <- simplify2array(
      foreach(i = 1:nsam) %dopar% {
        if (fam == 'Gaussian') {
          sslraagg_pal(rep(0, p), w, xstd[sswhich[, i], ], ycen[sswhich[, i]], Lad, dl, lam1, lam2, maxiter, cri)
        } else if (fam == 'Logistic') {
          sslogiaagg_pal(rep(0, p), xstd[sswhich[, i], ], y[sswhich[, i]], Lad, dl, w, lam1, lam2, intercept, maxiter, cri)
        } else {
          deltr <- y[sswhich[, i] + 1, 2]
          if (deltr[1] == 0){
            evest <- which(deltr == 1)[1]
            xtr <- xstd[sswhich[-(1:(evest - 1)), i] + 1, ]
            ytr <- y[sswhich[-(1:(evest - 1)), i] + 1]
          } else {
            xtr <- xstd[sswhich[, i] + 1, ]
            ytr <- y[sswhich[, i] + 1]
          }
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
          ssCoxaagg_pal(ntietr, rep(0, p), w, xtr, tietr, tie1tr, tie2tr, tie3tr, ytr[, 2], Lad, dl, lam1, lam2, maxiter, cri)
        }
      }
    )
  } else {
    outlist <- switch (fam,
                       Gaussian = sslraagg(rep(0, p), w, xstd, ycen, Lad, dl, lam1, lam2, sswhich, maxiter, cri),
                       Logistic = sslogiaagg(rep(0, p), xstd, y, Lad, dl, w, lam1, lam2, sswhich, intercept, maxiter, cri),
                       Cox = ssCoxaagg(ntie, rep(0, p), w, xstd, y[, 1], tie, tie1, tie2, tie3, y[, 2], Lad, dl, lam1, lam2, sswhich, maxiter, cri)
    )

  }

  out2 <- apply(outlist, 1:3, mean)
  ssmod <- 2 / p * apply(out2 * (1 - out2), 2:3, sum)


  ssmodbar <- apply(ssmod, 2, cummax)

  ncut <- ssmodbar
  ncut[ncut > beta] <- NA
  lamloc <- which(ncut == max(ncut, na.rm = T), arr.ind = T)[1, ]
  lambda1final <- lam1[lamloc[1]]
  lambda2final <- lam2[lamloc[2]]


  if (fam == 'Gaussian') {
    coeffinal <- lraagg(rep(0, p), lambda1final, lambda2final, w, xstd, ycen, Lad, dl, maxiter, cri)
    b0_final <- my-sum(meanx/sdx*coeffinal)
    b_final <- coeffinal/sdx
  } else if (fam == 'Logistic') {
    modfinal <- logiaagg(b0, rep(0, p), lambda1final, lambda2final, w, xstd, y, Lad, dl, intercept, maxiter, cri)
    b0 <- modfinal[1]
    b <- modfinal[-1]
    b0_final <- b0-sum(meanx/sdx*b)
    b_final <- b/sdx
  } else {
    modfinal <- Coxaagg(rep(0, p), lambda1final, lambda2final, ntie, w, xstd, tie, tie1, tie2, tie3, y[, 2], Lad, dl, maxiter, cri)
    b_final <- modfinal/sdx
    b0_final <- 0
  }
  n_final <- sum(b_final != 0)
  rownames(ssmodbar) <- lam1
  colnames(ssmodbar) <- lam2
  res <- new('ss', input = input, lambda1 = lam1, lambda2 = lam2, lambda1_ss = lambda1final,
             lambda2_ss = lambda2final, n_ss = n_final, ssm = ssmodbar, ssf = out2,
             intercept_ss = b0_final, coef_ss = b_final, fam = fam)
  structure(res, class = 'ss_glmaag')
}

##' @title Coefficients for ss_glmaag
##' @description Get the coefficients tuned by stability selection
##' @param object the model estimated via stability selection
##' @param ... \dots
##' @return the optimal coefficients get from stability selection including intercept (except for Cox)
##' @method coef ss_glmaag
##' @examples data(sampledata)
##' data(L0)
##' y <- sampledata$Y_Gau
##' x <- sampledata[, 4:6]
##' mod <- ss_glmaag(y, x, L0[seq_len(3), seq_len(3)], nsam = 3)
##' cc <- coef(mod)
##' @export
coef.ss_glmaag <- function(object, ...) {
  if(object@fam != 'Cox') {
    return(c(object@intercept_ss, object@coef_ss))
  } else {
    return(object@coef_ss)
  }
}

##' @title Prediction via stability selection
##' @description Predict using the model tuned by stability selection
##' @param object the ss_glmaag object
##' @param x the new dataset to be predicted, do training prediction if x is missing
##' @param type type of prediction (can be "link", or "reponse" ), ignored for Gaussian model. "link" is the linear predicted score, "response" is the predicted probability for logistic model and relative risk for Cox model
##' @param ... \dots
##' @return the predicted values
##' @method predict ss_glmaag
##' @examples data(sampledata)
##' data(L0)
##' y <- sampledata$Y_Gau
##' x <- sampledata[, 4:6]
##' mod <- ss_glmaag(y, x, L0[seq_len(3), seq_len(3)], nsam = 3)
##' pp <- predict(mod)
##' @export
predict.ss_glmaag <- function(object, x, type = 'link', ...) {
  if (missing(x)) {
    x <- object@input
  }
  x <- as.matrix(x)
  ypre <- as.vector(object@intercept_ss + x %*% object@coef_ss)
  if(type == 'link' | object@fam == 'Gaussian') {
    return(ypre)
  } else {
    if (object@fam == 'Logistic') {
      return(1 / (1 + exp(-ypre)))
    } else {
      return(exp(ypre))
    }
  }
}


##' @title Instability plot
##' @description Instability path plot
##' @param x the input ss_glmagrph object
##' @param ... \dots
##' @return the instability path
##' @importFrom data.table as.data.table melt setnames
##' @importFrom ggplot2 ggplot geom_line labs theme_bw theme aes_
##' @method plot ss_glmaag
##' @examples data(sampledata)
##' data(L0)
##' y <- sampledata$Y_Gau
##' x <- sampledata[, 4:6]
##' mod <- ss_glmaag(y, x, L0[seq_len(3), seq_len(3)], nsam = 3)
##' gg <- plot(mod)
##' @export
plot.ss_glmaag <- function(x, ...) {
  lambda1 <- NULL
  meanmat <- as.data.table(x@ssm)
  setnames(meanmat, as.character(x@lambda2))
  meanmat[, lambda1 := x@lambda1]
  meandat <- melt(meanmat, 'lambda1')
  setnames(meandat, c('lambda1', 'lambda2', 'Instability'))
  ggplot(meandat, aes_(x =~ lambda1, y =~ Instability, col =~ lambda2)) + geom_line() + labs(title = 'Instability', color = expression(lambda[2]~''), x = expression(lambda[1])) + theme_bw() + theme(axis.title.y = element_blank(), plot.title = element_text(hjust = .5), legend.position = 'bottom')
}

##' @title the results of the stability selection model
##' @description print fitted information
##' @param x the fitted ss_glmaag object
##' @param ... \dots
##' @method print ss_glmaag
##' @examples data(sampledata)
##' data(L0)
##' y <- sampledata$Y_Gau
##' x <- sampledata[, 4:6]
##' mod <- ss_glmaag(y, x, L0[seq_len(3), seq_len(3)], nsam = 3)
##' print(mod)
##' @export
print.ss_glmaag <- function(x, ...) {
  cat('The solution is: lambda1 = ', x@lambda1_ss, ', lambda2 = ', x@lambda2_ss, '.', sep = '')
}

