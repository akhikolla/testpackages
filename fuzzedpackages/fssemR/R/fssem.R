### FSSEM solver for multiple conditions and multiple eQTLs, specified as two.
##' @title multiFSSEMiPALM2
##' @description Implementing FSSELM algorithm for network inference. If Xs is identify for different conditions, multiFSSEMiPALM will be use, otherwise, please
##' use \code{multiFSSEMiPALM2} for general cases
##' @param Xs  eQTL matrices
##' @param Ys  Gene expression matrices
##' @param Bs  initialized GRN-matrices
##' @param Fs  initialized eQTL effect matrices
##' @param Sk  eQTL index of genes
##' @param sigma2 initialized noise variance from ridge regression
##' @param lambda Hyperparameter of lasso term in FSSEM
##' @param rho Hyperparameter of fused-lasso term in FSSEM
##' @param Wl  weight matrices for adaptive lasso terms
##' @param Wf  weight matrix for adaptive fused lasso term
##' @param p   number of genes
##' @param maxit maximum iteration number. Default 100
##' @param inert inertial function for iPALM. Default as k-1/k+2
##' @param threshold convergence threshold. Default 1e-6
##' @param verbose Default TRUE
##' @param sparse Sparse Matrix or not
##' @param strict Converge strictly or not. Default False
##' @param trans  Fs matrix is transposed to k x p or not. If Fs from ridge regression, trans = TRUE, else, trans = FALSE
##' @param B2norm B2norm matrices generated from ridge regression. Default NULL.
##' @return fit List of FSSEM model
##' \describe{
##' \item{Bs}{ coefficient matrices of gene regulatory networks}
##' \item{Fs}{ coefficient matrices of eQTL-gene effect}
##' \item{mu}{ Bias vector}
##' \item{sigma2}{ estimate of covariance in SEM}
##' }
##' @examples
##' seed = 1234
##' N = 100                                           # sample size
##' Ng = 5                                            # gene number
##' Nk = 5 * 3                                        # eQTL number
##' Ns = 1                                            # sparse ratio
##' sigma2 = 0.01                                     # sigma2
##' set.seed(seed)
##' library(fssemR)
##' data = randomFSSEMdata(n = N, p = Ng, k = Nk, sparse = Ns, df = 0.3, sigma2 = sigma2,
##'                        u = 5, type = "DG", nhub = 1, dag = TRUE)
##' ## If we assume that different condition has different genetics perturbations (eQTLs)
##' data$Data$X = list(data$Data$X, data$Data$X)
##' ## gamma = cv.multiRegression(data$Data$X, data$Data$Y, data$Data$Sk, ngamma = 20, nfold = 5,
##' ##                            N, Ng, Nk)
##' gamma = 0.6784248     ## optimal gamma computed by cv.multiRegression
##' fit   = multiRegression(data$Data$X, data$Data$Y, data$Data$Sk, gamma, N, Ng, Nk,
##'                       trans = FALSE)
##' Xs    = data$Data$X
##' Ys    = data$Data$Y
##' Sk    = data$Data$Sk
##'
##'
##' cvfitc <- cv.multiFSSEMiPALM2(Xs = Xs, Ys = Ys, Bs = fit$Bs, Fs = fit$Fs, Sk = Sk,
##'                              sigma2 = fit$sigma2, nlambda = 5, nrho = 5,
##'                              nfold = 5, p = Ng, q = Nk, wt = TRUE)
##'
##' fitc0 <- multiFSSEMiPALM2(Xs = Xs, Ys = Ys, Bs = fit$Bs, Fs = fit$Fs, Sk = Sk,
##'                           sigma2 = fit$sigma2, lambda = cvfitc$lambda, rho = cvfitc$rho,
##'                           Wl = inverseB(fit$Bs), Wf = flinvB(fit$Bs),
##'                           p = Ng, maxit = 100, threshold = 1e-5, sparse = TRUE, 
##'                           verbose = TRUE, trans = TRUE, strict = TRUE)
##'
##'
##' (TPR(fitc0$Bs[[1]], data$Vars$B[[1]]) + TPR(fitc0$Bs[[2]], data$Vars$B[[2]])) / 2
##' (FDR(fitc0$Bs[[1]], data$Vars$B[[1]]) + FDR(fitc0$Bs[[2]], data$Vars$B[[2]])) / 2
##' TPR(fitc0$Bs[[1]] - fitc0$Bs[[2]], data$Vars$B[[1]] - data$Vars$B[[2]])
##' FDR(fitc0$Bs[[1]] - fitc0$Bs[[2]], data$Vars$B[[1]] - data$Vars$B[[2]])
##' @export
multiFSSEMiPALM2 = function(Xs, Ys, Bs, Fs, Sk, sigma2, lambda, rho,
                           Wl, Wf, p, maxit = 100, inert = inert_opt("linear"), threshold = 1e-6,
                           verbose = TRUE, sparse = TRUE, trans = FALSE, B2norm = NULL,
                           strict = FALSE) {
  ## condition numbers
  m   = length(Ys)
  centered = proc.centerFSSEM2(Xs, Ys)
  Xs  = centered[["Xs"]]
  Ys  = centered[["Ys"]]
  mX  = centered[["mX"]]
  mY  = centered[["mY"]]
  n   = sapply(Ys, ncol)
  ## pre-processing fssem
  if (trans) {
    Fs  = lapply(Fs, t)
  }                                                                             # k x p
  Fx  = vector("list", m)
  f0  = vector("list", m)
  f1  = vector("list", m)
  Y2norm = vector("list", m)

  if(is.null(B2norm)) {
    B2norm = vector("list", m)
    for (k in 1:m) {
      B2norm[[k]] = colSums(Bs[[k]]**2)
      Bs[[k]] = matrix(0, nrow = p, ncol = p)
    }
  }

  for (k in 1:m) {
    Fx[[k]] = Fs[[k]] %*% Xs[[k]]
    Y2norm[[k]] = vector("numeric", p)
    f0[[k]] = vector("list", p)
    f1[[k]] = vector("list", p)
  }

  for (i in 1:p) {
    for (k in 1:m) {
      Xi = Xs[[k]][Sk[[i]], , drop = F]
      P  = solve(tcrossprod(Xi)) %*% Xi                                         # sk x n[k]
      yi = Ys[[k]][i, , drop = F]                                               # 1 x n[k] of gene i
      Yi = Ys[[k]][-i,]                                                         # (p-1) x n[k] of candidate regulators
      f0[[k]][[i]] = tcrossprod(P, yi)                                          # sk x 1
      f1[[k]][[i]] = tcrossprod(P, Yi)                                          # sk x (p-1)
      Y2norm[[k]][[i]] = tcrossprod(yi)
    }
  }

  ## block coordinate descent of FSSEM(column-major)
  niter = 1
  ImBs  = lapply(Bs, function(B){ diag(p) - B })
  Dets  = sapply(ImBs, det)
  IBinv = lapply(ImBs, solve)
  sigma2 = sigma2FSSEM2(Xs, Ys, Bs, Fs, n, p, m)
  L      = logLikFSSEM(Bs, Wl, Wf, lambda, rho, sigma2, Dets, n, p)
  minit = 1
  nonincreasing = TRUE
  Bs_last = list(Bs, Bs)
  while (niter <= maxit) {
    t  = inert(niter)
    Bt = lapply(1:m, function(k) {
      Bs_last[[2]][[k]] + t * (Bs_last[[2]][[k]] - Bs_last[[1]][[k]])
    })
    Fs_last = Fs
    L_last  = L
    if (!nonincreasing) {
      if (minit <= 1e6) {
        minit = min(minit * 1.02, 1e6)
      } else {
        break
      }
    }
    for (i in 1:p) {
      ci = lapply(IBinv, function(IBi) { IBi[i, , drop = F] })
      bi = lapply(Bt, function(B) { B[, i, drop = F] })
      Ri = lapply(1:m, function(k) {
        Ys[[k]] - Bs[[k]][,-i] %*% Ys[[k]][-i,] - Fx[[k]]
      })
      gi = lapply(1:m, function(k) {
        grad = cwiseGradient4FSSEM(n[k], ci[[k]], Ys[[k]][i, ,drop = F], Ri[[k]], Y2norm[[k]][i],sigma2[1])
        grad(bi[[k]])[-i, , drop = F]
      })
      # `cwiseLipschitz4FSSEM` is the improved version of `lips_cwise_SML` and `lips_rwise_SML` in
      # ./inst/00_SparsemaximumLiklihood.R; `multiFSSEMiPALM2` is equivalent to `genSML_iPALM` in
      # deprecated version @ https://github.com/Ivis4ml/FSSEM
      ## Li = sapply(1:m, function(k) {
      ##   U = crossprod(ImBs[[k]][, -i])
      ##   cwiseLipschitzFSSEMv0(n[k], chol2inv(U),
      ##                    ImBs[[k]][-i, -i],
      ##                    ImBs[[k]][i, -i, drop = F],
      ##                    sum((ci[[k]] * Dets[[k]]) ^ 2),
      ##                    det(U),
      ##                    Y2norm[[k]][i],
      ##                    sigma2, p)[1]
      ##})
      Li = sapply(1:m, function(k) {
        z  = ci[[k]][,i] * Dets[[k]]
        c2 = sum((ci[[k]][, -i] * Dets[[k]])**2)
        b2 = B2norm[[k]][i]
        ## cwiseLipschitz4FSSEM(n[k], z, c2, b2, Y2norm[[k]][i], sigma2)
        cwiseLipschitz4FSSEM2B(n[k], z, c2, b2, Y2norm[[k]][i], sigma2, ImBs[[k]], i)
      })
      Li = norm(Li, "2")
      ## Armoji-scheme line-search
      tau  = minit
      while (TRUE) {
        ui = lapply(1:m, function(k) { bi[[k]][-i, ] - 1 / (tau * Li) * gi[[k]] })
        w  = list(Wl[[1]][-i, i], Wl[[2]][-i, i])
        r  = Wf[-i, i]
        xi = proxFLSA(ui, w, r, lambda, rho, tau * Li)
        det_update = sapply(1:m, function(k) { IBinv[[k]][i, i] - tcrossprod(xi[[k]], IBinv[[k]][i, -i])[1] })
        if (all(det_update != 0))
          break
        tau = tau * 1.02
      }
      for(k in 1:m) {
        Bs[[k]][-i, i] = xi[[k]]
        ImBs[[k]] = diag(p) - Bs[[k]]
        Dets[k]   = tcrossprod(ImBs[[k]][,i], IBinv[[k]][i, , drop=F])[1] * Dets[k]
        deltaB  = Bs_last[[2]][[k]][, i, drop=F] - Bs[[k]][, i, drop=F]
        IBinv[[k]] = IBinv[[k]] - IBinv[[k]] %*% deltaB %*%  IBinv[[k]][i, , drop=F] / (1 + IBinv[[k]][i, , drop=F] %*% deltaB)[1]
      }
    }
    ## update Fs
    for (i in 1:p) {
      for (k in 1:m) {
        Fs[[k]][i, Sk[[i]]] = f0[[k]][[i]] - f1[[k]][[i]] %*% Bs[[k]][i, -i]
      }
    }
    for (k in 1:m) {
      Fx[[k]] = Fs[[k]] %*% Xs[[k]]
    }
    ## update sigma2
    sigma2 = sigma2FSSEM2(Xs, Ys, Bs, Fs, n, p, m)
    ## converged
    Berr = sum(sapply(1:m, function(k) { norm(Bs[[k]] - Bs_last[[2]][[k]], "f") / (1 + norm(Bs_last[[2]][[k]], "f")) }))
    Ferr = sum(sapply(1:m, function(k) { norm(Fs[[k]] - Fs_last[[k]], "f") / (1 + norm(Fs_last[[k]], "f"))}))
    L = logLikFSSEM(Bs, Wl, Wf, lambda, rho, sigma2, Dets, n, p)
    Lerr = abs(L_last - L) / (1 + abs(L_last))
    nonincreasing = (L_last > L)
    converged = if (strict) {
      (Berr + Ferr) <= threshold && Lerr <= threshold
    } else {
      Lerr <= threshold
    }
    if (verbose) {
      cat(sprintf("FSSEM: niter = %d, err_param = %f, loglik = %f, sigma2 = %f\n", niter, Berr + Ferr, L, sigma2))
    }
    niter = niter + 1
    Bs_last = list(Bs_last[[2]], Bs)
    if (converged || niter > maxit || sigma2 < 1e-6) {
      mu = lapply(1:m, function(k) {
        (diag(p) - Bs[[k]]) %*% mY[[k]] - sapply(1:p, function(i) { mX[[k]][Sk[[i]]] %*% Fs[[k]][i, Sk[[i]]] })
      })
      if (sparse) {
        Bs = lapply(Bs, Matrix, sparse = T)
      }
      break
    }
  }
  list(Bs = Bs, Fs = Fs, mu = mu, sigma2 = sigma2, Dets = Dets, trans = TRUE)
}


## initialization of lambda parameter
##' @title initLambdaiPALM2
##' @param Xs      eQTL matrices
##' @param Ys      Gene expression matrices
##' @param Bs      initialized GRN-matrices
##' @param Fs      initialized eQTL effect matrices
##' @param Sk      eQTL index of genes
##' @param sigma2  initialized noise variance
##' @param Wl  weight matrices for adaptive lasso terms
##' @param Wf  weight matrix for adaptive fused lasso term
##' @param p   number of genes
##' @param k   number of eQTL
##' @return lambda_max
initLambdaiPALM2 = function(Xs, Ys, Bs, Fs, Sk, sigma2, Wl, Wf, p, k) {
  m   = length(Ys)
  centered = proc.centerFSSEM2(Xs, Ys)
  Xs  = centered[["Xs"]]
  Ys  = centered[["Ys"]]
  mX  = centered[["mX"]]
  mY  = centered[["mY"]]
  n   = sapply(Ys, ncol)
  Fi  = lapply(1:m, function(i){ matrix(0, p, k) })
  Bs  = lapply(1:m, function(i){ matrix(0, p, p) })
  Fx  = vector("list", m)
  f0  = vector("list", m)
  f1  = vector("list", m)
  Y2norm = vector("list", m)
  for (i in 1:p) {
    for (j in 1:m) {
      Xi = Xs[[j]][Sk[[i]], , drop = F]
      P  = solve(tcrossprod(Xi)) %*% Xi                                         # sk x n[k]
      yi = Ys[[j]][i, , drop = F]                                               # 1 x n[k] of gene i
      Fi[[j]][i, Sk[[i]]] = tcrossprod(P, yi)
    }
  }
  sigma2 = sigma2FSSEM2(Xs, Ys, Bs, Fi, n, p, m)
  res = vector("list", m)
  for (j in 1:m) {
    res[[j]] = abs(1 / sigma2 * tcrossprod(Ys[[j]] - Fi[[j]] %*% Xs[[j]], Ys[[j]])) / Wl[[j]]
  }
  lambda_max = max(sapply(res, max))
  while (TRUE) {
    fit = multiFSSEMiPALM2(Xs, Ys, Bs, Fs, Sk, sigma2, lambda = lambda_max, rho = 0, Wl, Wf, p = p,
                          maxit = 20, threshold = 1e-3, verbose = FALSE, sparse = FALSE, trans = TRUE)
    if (all(fit$Bs[[1]] == 0) && all(fit$Bs[[2]] == 0)) {
      lambda_max = lambda_max / 1.02
    } else {
      break
    }
  }
  lambda_max
}

## initialization of rho parameter
##' @title initRhoiPALM2
##' @param Xs      eQTL matrices
##' @param Ys      Gene expression matrices
##' @param Bs      initialized GRN-matrices
##' @param Fs      initialized eQTL effect matrices
##' @param Sk      eQTL index of genes
##' @param sigma2  initialized noise variance
##' @param Wl  weight matrices for adaptive lasso terms
##' @param Wf  weight matrix for adaptive fused lasso term
##' @param lambda  lambda w.r.t. rho_max
##' @param n       number of observations
##' @param p       number of genes
##' @return rho_max
initRhoiPALM2 = function(Xs, Ys, Bs, Fs, Sk, sigma2, Wl, Wf, lambda, n, p) {
  fit = multiFSSEMiPALM2(Xs, Ys, Bs, Fs, Sk, sigma2, lambda = lambda, rho = Inf, Wl, Wf, p = p,
                        maxit = 20, threshold = 1e-3, verbose = FALSE, sparse = FALSE, trans = TRUE)
  m = length(Ys)
  res = vector("list", m)
  for(k in 1:m) {
    ImB = diag(p) - fit$Bs[[k]]
    d = n[k] * solve(ImB) - tcrossprod(ImB %*% Ys[[k]] - fit$Fs[[k]] %*% Xs[[k]] - tcrossprod(fit$mu[[k]], rep(1, n[k])), Ys[[k]]) / fit$sigma2
    d = d + lambda * Wl[[k]] * sign(fit$Bs[[k]])
    res[[k]] = abs(d * sign(fit$Bs[[k]]))
    res[[k]] = res[[k]] / Wf
    diag(res[[k]]) = 0
  }
  max(sapply(res, max))
}


## cross-validation function for SEM model
##' @title cv.multiFSSEMiPALM2
##' @param Xs      eQTL matrices
##' @param Ys      Gene expression matrices
##' @param Bs      initialized GRN-matrices
##' @param Fs      initialized eQTL effect matrices
##' @param Sk      eQTL index of genes
##' @param sigma2  initialized noise variance
##' @param nlambda number of hyper-parameter of lasso term in CV
##' @param nrho    number of hyper-parameter of fused-lasso term in CV
##' @param nfold   CVfold number. Default 5/10
##' @param p       number of genes
##' @param q       number of eQTLs
##' @param wt  use adaptive lasso or not. Default TRUE.
##' @param plot    plot contour of cvmean or not. Default FALSE.
##' @return list of cross-validation result
##' @export
cv.multiFSSEMiPALM2 = function(Xs, Ys, Bs, Fs, Sk, sigma2, nlambda = 20, nrho = 20,
                              nfold = 5, p, q, wt = TRUE, plot = FALSE) {
  m   = length(Ys)
  n   = sapply(Ys, ncol)
  B2norm = vector("list", m)
  for (k in 1:m) {
    B2norm[[k]] = colSums(Bs[[k]] ** 2)
  }
  if (wt) {
    Wl = inverseB(Bs)
    Wf = flinvB(Bs)
  } else {
    Wl = inverseB(Bs)
    Wf = floneB(Bs)
  }

  lambda_max  = initLambdaiPALM2(Xs, Ys, Bs, Fs, Sk, sigma2, Wl, Wf, p, q)
  lambda_factors = 10 ** seq(0, -3, length.out = nlambda) * lambda_max
  rho_factors = vector("list", nlambda)
  for (i in 1:nlambda) {
    lambda  = lambda_factors[i]
    rho_max = initRhoiPALM2(Xs, Ys, Bs, Fs, Sk, sigma2, Wl, Wf, lambda, n, p)
    rho_factors[[i]] = 10 ** seq(0, -3, length.out = nrho) * rho_max
  }
  ## cv on nfold
  params = vector("list", nlambda * nrho)
  cverr  = vector("list", nlambda * nrho)
  cvfold = lapply(1:m, function(i) {
    sample(seq(1, nfold), size = n[i], replace = T)
  })
  Xtrain = vector("list", nfold)
  Xtest  = vector("list", nfold)
  Ytrain = vector("list", nfold)
  Ytest  = vector("list", nfold)
  for (i in 1:nfold) {
    Xtrain[[i]] = lapply(1:m, function(mi) {
      Xs[[mi]][, cvfold[[mi]] != i, drop = F]
    })
    Xtest[[i]]  = lapply(1:m, function(mi) {
      Xs[[mi]][, cvfold[[mi]] == i, drop = F]
    })
    Ytrain[[i]] = lapply(1:m, function(mi) {
      Ys[[mi]][, cvfold[[mi]] != i, drop = F]
    })
    Ytest[[i]]  = lapply(1:m, function(mi) {
      Ys[[mi]][, cvfold[[mi]] == i, drop = F]
    })
  }
  for (i in 1:nfold) {
    fit = list()
    j  = 1
    for (ilambda in 1:nlambda) {
      for (irho in 1:nrho) {
        lambda = lambda_factors[ilambda]
        rho    = rho_factors[[ilambda]][irho]
        cat(sprintf("lambda = %4f, rho = %4f, ncv = %d\n", lambda, rho, i))
        if (j %% nrho == 1) {
          fit[[j]] = multiFSSEMiPALM2(Xtrain[[i]], Ytrain[[i]], Bs, Fs, Sk, sigma2, lambda, rho, Wl, Wf, p,
                                     maxit = 50, threshold = 1e-5, sparse = FALSE, trans = TRUE, verbose = FALSE)
        } else {
          if (rho_factors[[ilambda]][irho] == rho_factors[[ilambda]][irho - 1]) {
            fit[[j]] = fit[[j-1]]
          } else {
            fit[[j]] = multiFSSEMiPALM2(Xtrain[[i]], Ytrain[[i]], Bs, Fs, Sk, fit[[j-1]]$sigma2, lambda, rho, Wl, Wf, p,
                                     maxit = 50, threshold = 1e-5, sparse = FALSE, trans = TRUE, verbose = FALSE)
          }
        }
        err = logLiklihood2(Xtest[[i]], Ytest[[i]], fit[[j]]$Bs, fit[[j]]$Fs, fit[[j]]$mu, fit[[j]]$Dets, fit[[j]]$sigma2, p)
        if (i == 1) {
          params[[j]] = c(lambda, rho)
        }
        cverr[[j]] = c(cverr[[j]], err)
        j = j + 1
      }
    }
  }
  cvmean = data.frame(lambda = sapply(params, `[`, 1), rho = sapply(params, `[`, 2),
                      cvmean = sapply(cverr, mean), cvsd = sapply(cverr, sd))
  cvmin  = which.min(cvmean$cvmean)
  cvms   = matrix(0, nrow = nlambda, ncol = nrho)
  for (i in 1:nlambda) {
    for (j in 1:nrho) {
      cvms[i, j] = cvmean$cvmean[(i - 1) * nrho + j]
    }
  }
  list(lambda = cvmean[cvmin, 1], rho = cvmean[cvmin, 2],
       cvm = cvms, res = cverr, params = params)
}

## optimized model selection of FSSEM selected from BIC/eBIC for multiple eQTL source
##' @title opt.multiFSSEMiPALM2
##' @description optimize multiFSSEMiPALM's parameters by minimize BIC,
##' when feature size is large (> 300), BIC methods will be much faster than Cross-validation
##' @param Xs      eQTL matrices
##' @param Ys      Gene expression matrices
##' @param Bs      initialized GRN-matrices
##' @param Fs      initialized eQTL effect matrices
##' @param Sk      eQTL index of genes
##' @param sigma2  initialized noise variance
##' @param nlambda number of hyper-parameter of lasso term in CV
##' @param nrho    number of hyper-parameter of fused-lasso term in CV
##' @param p       number of genes
##' @param q       number of eQTLs
##' @param wt  use adaptive lasso or not. Default TRUE.
##' @return list of model selection result
##' @export
##' @examples
##' seed = 1234
##' N = 100                                           # sample size
##' Ng = 5                                            # gene number
##' Nk = 5 * 3                                        # eQTL number
##' Ns = 1                                            # sparse ratio
##' sigma2 = 0.01                                     # sigma2
##' set.seed(seed)
##' library(fssemR)
##' data = randomFSSEMdata(n = N, p = Ng, k = Nk, sparse = Ns, df = 0.3, sigma2 = sigma2,
##'                        u = 5, type = "DG", nhub = 1, dag = TRUE)
##' ## If we assume that different condition has different genetics perturbations (eQTLs)
##' data$Data$X = list(data$Data$X, data$Data$X)
##' ## gamma = cv.multiRegression(data$Data$X, data$Data$Y, data$Data$Sk, ngamma = 20, nfold = 5,
##' ##                            N, Ng, Nk)
##' gamma = 0.6784248     ## optimal gamma computed by cv.multiRegression
##' fit   = multiRegression(data$Data$X, data$Data$Y, data$Data$Sk, gamma, N, Ng, Nk,
##'                       trans = FALSE)
##' Xs    = data$Data$X
##' Ys    = data$Data$Y
##' Sk    = data$Data$Sk
##'
##' fitm <- opt.multiFSSEMiPALM2(Xs = Xs, Ys = Ys, Bs = fit$Bs, Fs = fit$Fs, Sk = Sk,
##'                            sigma2 = fit$sigma2, nlambda = 10, nrho = 10,
##'                            p = Ng, q = Nk, wt = TRUE)
##'
##' fitc0 <- fitm$fit
##'
##' (TPR(fitc0$Bs[[1]], data$Vars$B[[1]]) + TPR(fitc0$Bs[[2]], data$Vars$B[[2]])) / 2
##' (FDR(fitc0$Bs[[1]], data$Vars$B[[1]]) + FDR(fitc0$Bs[[2]], data$Vars$B[[2]])) / 2
##' TPR(fitc0$Bs[[1]] - fitc0$Bs[[2]], data$Vars$B[[1]] - data$Vars$B[[2]])
##' FDR(fitc0$Bs[[1]] - fitc0$Bs[[2]], data$Vars$B[[1]] - data$Vars$B[[2]])
opt.multiFSSEMiPALM2 = function(Xs, Ys, Bs, Fs, Sk, sigma2, nlambda = 20, nrho = 20,
                               p, q, wt = TRUE) {
  m   = length(Ys)
  n   = sapply(Ys, ncol)
  B2norm = vector("list", m)
  for (k in 1:m) {
    B2norm[[k]] = colSums(Bs[[k]] ** 2)
  }
  if (wt) {
    Wl = inverseB(Bs)
    Wf = flinvB(Bs)
  } else {
    Wl = inverseB(Bs)
    Wf = floneB(Bs)
  }

  lambda_max  = initLambdaiPALM2(Xs, Ys, Bs, Fs, Sk, sigma2, Wl, Wf, p, q)
  lambda_factors = 10 ** seq(0,-3, length.out = nlambda) * lambda_max
  rho_factors = vector("list", nlambda)
  for (i in 1:nlambda) {
    lambda  = lambda_factors[i]
    rho_max = initRhoiPALM2(Xs, Ys, Bs, Fs, Sk, sigma2, Wl, Wf, lambda, n, p)
    rho_factors[[i]] = 10 ** seq(0,-3, length.out = nrho) * rho_max
  }
  params = vector("list", nlambda * nrho)
  fit = list()
  j   = 1
  for (ilambda in 1:nlambda) {
      for (irho in 1:nrho) {
        lambda = lambda_factors[ilambda]
        rho    = rho_factors[[ilambda]][irho]
        cat(sprintf("FSSEM@lambda = %4f, rho = %4f\n", lambda, rho))
        if (j %% nrho == 1) {
          fit[[j]] = multiFSSEMiPALM2(Xs, Ys, Bs, Fs, Sk, sigma2, lambda, rho, Wl, Wf, p,
                                     maxit = 50, threshold = 1e-5, sparse = FALSE, trans = TRUE, verbose = FALSE)
        } else {
          if (rho_factors[[ilambda]][irho] == rho_factors[[ilambda]][irho - 1]) {
            fit[[j]] = fit[[j-1]]
          } else {
            fit[[j]] = multiFSSEMiPALM2(Xs, Ys, Bs, Fs, Sk, sigma2, lambda, rho, Wl, Wf, p,
                                     maxit = 50, threshold = 1e-5, sparse = FALSE, trans = TRUE, verbose = FALSE)
          }
        }
        err = bayesianInfocriterion2(Xs, Ys, fit[[j]]$Bs, fit[[j]]$Fs, fit[[j]]$mu, fit[[j]]$Dets, fit[[j]]$sigma2, p)
        params[[j]] = c(lambda, rho, err)
        j = j + 1
      }
  }
  BICmat = data.frame(do.call(rbind, params))
  colnames(BICmat) = c("lambda", "rho", "BIC")
  BICmin = which.min(BICmat$BIC)
  fitmin = fit[[BICmin]]
  for(i in 1:m) {
    fitmin$Bs[[i]] = Matrix(fitmin$Bs[[i]], sparse = TRUE)
    fitmin$Fs[[i]] = Matrix(fitmin$Fs[[i]], sparse = TRUE)
  }
  list(lambda = BICmat[BICmin, 1], rho = BICmat[BICmin, 2], fit = fitmin,
       minBIC = min(BICmat$BIC), fit = fit, BICs = BICmat)
}

