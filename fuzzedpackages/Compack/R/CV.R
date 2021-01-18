#' @title
#' Cross-validation for FuncompCGL.
#'
#' @description
#' k-fold cross-validation for FuncompCGL; produce a plot and return
#' optimal values of \code{lam} and \code{k}.
#'
#' @usage
#' cv.FuncompCGL(y, X, Zc = NULL, lam = NULL, nlam = 100, k = 4:10, ref = NULL,
#'               foldid, nfolds = 10, W = rep(1,times = p - length(ref)),
#'               trim = 0, outer_maxiter = 1e+06, keep = FALSE, ...)
#'
#' @param X a data frame or matrix.
#'       \itemize{
#'       \item If \code{nrow(X)} > \eqn{n},
#'             \code{X} should be a data frame or matrix of the functional compositional
#'             predictors with \eqn{p} columns for the values of the compositional components,
#'             one column indicating the subject ID and one column of observed time points.
#'             The order of the Subject ID should be the SAME as that of \code{y}.
#'       \item If \code{nrow(X)[1]}=\eqn{n},
#'             \code{X} is considered as the integrated design matrix, a
#'             \code{n*(k*p - length(ref))} matrix.
#'       }
#'
#' @param k a vector of integer values of the degrees of freedom; default is 4:10.
#'
#' @param nfolds number of folds, default is 10. The smallest allowable value is \code{nfolds=3}.
#'
#' @param foldid an optional vector of values between 1 and the sample size \code{n}, providing the fold
#'             assignments. If supplied, \code{nfold} can be missing.
#'
#' @param trim percentage to be trimmed off the prediction errors from either side; default is 0.
#'
#' @param outer_maxiter maximum number of loops allowed for the augmented Lagrange method.
#'
#' @param keep If \code{keep=TRUE}, fitted models in cross validation are reported.
#'             Default is \code{keep=FALSE}.
#'
#' @param ... other arguments that can be passed to \code{\link{FuncompCGL}}.
#'
#' @inheritParams FuncompCGL
#'
#' @return An object of S3 class \code{"cv.FuncompCGL"} is return, which is a list
#' containing:
#' \item{FuncompCGL.fit}{a list of length \code{length(k)},
#'                       with elements being the fitted \code{\link{FuncompCGL}} objects of different
#'                       degrees of freedom.}
#'
#' \item{lam}{the sequence of \code{lam}.}
#'
#' \item{Ftrim}{a list for cross validation results with trim = 0.
#'                \itemize{
#'                \item \code{cvm} the mean cross-validated error
#'                                 - a matrix of dimension \code{length(k)*length(lam).}
#'                \item \code{cvsd} estimated standard error of \code{cvm}.
#'                \item \code{cvup} upper curve = \code{cvm + cvsd}.
#'                \item \code{cvlo} lower curve = \code{cvm - cvsd}.
#'                \item \code{lam.min} the optimal values of \code{k} and \code{lam}
#'                      that give minimum cross validation error \code{cvm}.
#'                \item \code{lam.1se} the optimal values of \code{k} and \code{lam}
#'                      that give cross validation error withnin 1 standard error of
#'                      the miminum \code{cvm}.
#'                }
#'              }
#'
#' \item{Ttrim}{a list of cross validation result with \code{trim*100\%}. The structure is the
#'              same as that for \code{Ftrim}.}
#'
#' \item{fit.preval, foldid}{\code{fit.preval} is the array of fitted models.
#'                           Only kept when \code{keep=TRUE}.}
#'
#'
#' @examples
#' \donttest{
#' ## generate training and testing data
#' df_beta = 5
#' p = 30
#' beta_C_true = matrix(0, nrow = p, ncol = df_beta)
#' beta_C_true[1, ] <- c(-0.5, -0.5, -0.5 , -1, -1)
#' beta_C_true[2, ] <- c(0.8, 0.8,  0.7,  0.6,  0.6)
#' beta_C_true[3, ] <- c(-0.8, -0.8 , 0.4 , 1 , 1)
#' beta_C_true[4, ] <- c(0.5, 0.5, -0.6  ,-0.6, -0.6)
#'
#' n_train = 50
#' n_test = 30
#' nfolds = 5
#' foldid <- sample(rep(seq(nfolds), length = n_train))
#' k_list <- c(4,5)
#'
#' Data <- Fcomp_Model(n = n_train, p = p, m = 0, intercept = TRUE,
#'                     SNR = 4, sigma = 3, rho_X = 0.2, rho_T = 0.5,
#'                     df_beta = df_beta, n_T = 20, obs_spar = 1, theta.add = FALSE,
#'                     beta_C = as.vector(t(beta_C_true)))
#' arg_list <- as.list(Data$call)[-1]
#' arg_list$n <- n_test
#' Test <- do.call(Fcomp_Model, arg_list)
#'
#' ## cv_cgl: Constrained group lasso
#' cv_cgl <-  cv.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
#'                          Zc = Data$data$Zc, intercept = Data$data$intercept,
#'                          k = k_list, foldid = foldid,
#'                          keep = TRUE)
#' plot(cv_cgl,k = k_list)
#' cv_cgl$Ftrim[c("lam.min", "lam.1se")]
#' beta <-  coef(cv_cgl, trim = FALSE, s = "lam.min")
#' k_opt <- cv_cgl$Ftrim$lam.min['df']
#' ## plot path against L2-norm of group coefficients
#' plot(cv_cgl$FuncompCGL.fit[[as.character(k_opt)]])
#' ## or plot path against L1-norm of group coefficients
#' plot(cv_cgl$FuncompCGL.fit[[as.character(k_opt)]], ylab = "L1")
#'
#' m1 <- ifelse(is.null(ncol(Data$data$Zc)), 0, ncol(Data$data$Zc))
#' m1 <- m1 + Data$data$intercept
#' if(k_opt == df_beta) {
#'   plot(Data$beta, col = "red", pch = 19,
#'        ylim = range(c(range(Data$beta), range(beta))))
#'   abline(v= seq(from = 0, to = (p*df_beta), by = df_beta ))
#'   abline(h = 0)
#'   points(beta)
#'   if(m1 > 0) points(p*df_beta + 1:m1, tail(Data$beta, m1),
#'                     col = "blue", pch = 19)
#' } else {
#'   plot(beta, ylim = range(c(range(Data$beta), range(beta))) )
#'   abline(v= seq(from = 0, to = (p*k_opt), by = k_opt ))
#'   abline(h = 0, col = "red")
#'   if(m1 > 0) points(p*k_opt + 1:m1, tail(Data$beta, m1),
#'                     col = "blue", pch = 19)
#' }
#'
#' beta_C <- matrix(beta[1:(p*k_opt)], byrow = TRUE, nrow = p)
#' ## satisfies zero-sum constraints
#' cat("colSums:", colSums(beta_C))
#' Nonzero <- (1:p)[apply(beta_C, 1, function(x) max(abs(x)) >0)]
#' cat("selected groups:", Nonzero)
#'
#' oldpar <- par(mfrow=c(2,1))
#' sseq <- Data$basis.info[, 1]
#' beta_curve_true <- Data$basis.info[, -1] %*%  t(beta_C_true)
#' Nonzero_true <- (1:p)[apply(beta_C_true, 1, function(x) max(abs(x)) >0)]
#' matplot(sseq, beta_curve_true, type = "l", ylim = range(beta_curve_true),
#'         ylab = "True coeffcients curves", xlab = "TIME")
#' abline(a = 0, b = 0, col = "grey", lwd = 2)
#' text(0, beta_curve_true[1, Nonzero_true], labels = Nonzero_true)
#'
#' beta_curve <- splines::bs(sseq, df = k_opt, intercept = TRUE) %*% t(beta_C)
#' matplot(sseq, beta_curve, type = "l", ylim = range(beta_curve_true),
#'         ylab = "Estimated coefficient curves", xlab = "TIME")
#' abline(a = 0, b = 0, col = "grey", lwd = 2)
#' text(0, beta_curve[1, Nonzero], labels = Nonzero)
#' par(oldpar)
#'
#' ## plot L1-norm of the estimated coefficients for each component of the composition
#' plot(apply(abs(beta_C),1,sum), ylab = "L1-norm", xlab = "Component index")
#' ## or plot L2-norm
#' plot(apply(abs(beta_C),1, function(x) sqrt(sum(x^2))),
#'      ylab = "L2-norm", xlab = "Component index")
#'
#' ## set a thresholding for variable selection via cross-validation model
#' ## example 1: cut by average L2-norm for estimated coefficient curves
#' Curve_L2 <- colSums(beta_curve^2)
#' Curve_L2 <- Curve_L2 - colSums(beta_curve[c(1, nrow(beta_curve)), ]^2) / 2
#' Curve_L2 <- Curve_L2 * (Data$basis.info[2,1] - Data$basis.info[1,1])
#' Curve_L2 <- sqrt(Curve_L2)
#' plot(Curve_L2, xlab = "Component index", ylab = "L2-norm for coefficient curves")
#' cutoff <- sum(Curve_L2) / p
#' Nonzero_cut <- (1:p)[which(Curve_L2 >= cutoff)]
#' Nonzero_cut
#' ## example 2: cut by average L2-norm for estimated coefficient vectors
#' cutoff <- sum(apply(beta_C, 1, function(x) norm(x, "2")))/p
#' Nonzero_cut2 <- (1:p)[apply(beta_C, 1, function(x, a) norm(x, "2") >= a, a = cutoff)]
#' ## example 3: cut by average L1-norm for estimated coefficient vectors
#' cutoff <- sum(abs(beta_C))/p
#' Nonzero_cut3 <- (1:p)[apply(beta_C, 1, function(x, a) sum(abs(x)) >= a, a = cutoff)]
#'
#' y_hat <- predict(cv_cgl, Data$data$Comp, Data$data$Zc, s = "lam.min")
#' MSE <- sum((drop(Data$data$y) - y_hat)^2) / n_train
#' y_hat <- predict(cv_cgl, Test$data$Comp, Test$data$Zc, s = "lam.min")
#' PRE <- sum((drop(Test$data$y) - y_hat)^2) / n_test
#' cgl_result <- list(cv.result = cv_cgl, beta = beta,
#'                    Nonzero = c("Original" = Nonzero, "Cut" = Nonzero_cut),
#'                    MSE = MSE, PRE = PRE)
#'
#' ## cv_naive: ignoring the zero-sum constraints
#' ## set mu_raio = 0 to identifying without linear constraints,
#' ## no outer_loop for Lagrange augmented multiplier
#' cv_naive <-  cv.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
#'                            Zc = Data$data$Zc, intercept = Data$data$intercept,
#'                            k = k_list, foldid = foldid, keep = TRUE,
#'                            mu_ratio = 0)
#' plot(cv_naive, k = k_list)
#' beta <-  coef(cv_naive, trim = FALSE, s = "lam.min")
#' k_opt <- cv_naive$Ftrim$lam.min['df']
#' beta_C <- matrix(beta[1:(p*k_opt)], byrow = TRUE, nrow = p)
#' ## does NOT satisfy zero-sum constraints
#' cat("colSums:", colSums(beta_C))
#' Nonzero <- (1:p)[apply(beta_C, 1, function(x) max(abs(x)) >0)]
#' beta_curve <- splines::bs(sseq, df = k_opt, intercept = TRUE) %*% t(beta_C)
#' Curve_L2 <- colSums(beta_curve^2) - colSums(beta_curve[c(1, nrow(beta_curve)), ]^2) / 2
#' Curve_L2 <- sqrt(Curve_L2 * (Data$basis.info[2,1] - Data$basis.info[1,1]))
#' cutoff <- sum(Curve_L2) / p
#' Nonzero_cut <- (1:p)[which(Curve_L2 >= cutoff)]
#' y_hat <- predict(cv_naive, Data$data$Comp, Data$data$Zc, s = "lam.min")
#' MSE <- sum((drop(Data$data$y) - y_hat)^2) / n_train
#' y_hat <- predict(cv_naive, Test$data$Comp, Test$data$Zc, s = "lam.min")
#' PRE <- sum((drop(Test$data$y) - y_hat)^2) / n_test
#' naive_result <- list(cv.result = cv_naive, beta = beta,
#'                      Nonzero = c("Original" = Nonzero, "Cut" = Nonzero_cut),
#'                      MSE = MSE, PRE = PRE)
#'
#' ## cv_base: random select a component as reference
#' ## mu_ratio is set to 0 automatically once ref is set to a integer
#' ref = sample(1:p, 1)
#' cv_base <- cv.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
#'                          Zc = Data$data$Zc, intercept = Data$data$intercept,
#'                          k = k_list, foldid = foldid, keep = TRUE,
#'                          ref = ref)
#' plot(cv_base, k = k_list)
#' beta <-  coef(cv_base, trim = FALSE, s = "lam.min")
#' k_opt <- cv_base$Ftrim$lam.min['df']
#' beta_C <- matrix(beta[1:(p*k_opt)], byrow = TRUE, nrow = p)
#' ## satisfies zero-sum constraints
#' cat("colSums:", colSums(beta_C))
#' Nonzero <- (1:p)[apply(beta_C, 1, function(x) max(abs(x)) >0)]
#' beta_curve <- splines::bs(sseq, df = k_opt, intercept = TRUE) %*% t(beta_C)
#' Curve_L2 <- colSums(beta_curve^2) - colSums(beta_curve[c(1, nrow(beta_curve)), ]^2) / 2
#' Curve_L2 <- sqrt(Curve_L2 * (Data$basis.info[2,1] - Data$basis.info[1,1]))
#' cutoff <- sum(Curve_L2) / p
#' Nonzero_cut <- (1:p)[which(Curve_L2 >= cutoff)]
#' y_hat <- predict(cv_base, Data$data$Comp, Data$data$Zc, s = "lam.min")
#' MSE <- sum((drop(Data$data$y) - y_hat)^2) / n_train
#' y_hat <- predict(cv_base, Test$data$Comp, Test$data$Zc, s = "lam.min")
#' PRE <- sum((drop(Test$data$y) - y_hat)^2) / n_test
#' base_result <- list(cv.result = cv_base, beta = beta,
#'                     Nonzero = c("Original" = Nonzero, "Cut" = Nonzero_cut),
#'                     MSE = MSE, PRE = PRE)
#' }
#'
#' @seealso
#' \code{\link{FuncompCGL}} and \code{\link{GIC.FuncompCGL}},
#' and \code{\link[=predict.cv.FuncompCGL]{predict}}, \code{\link[=coef.cv.FuncompCGL]{coef}}
#' and \code{\link[=plot.cv.FuncompCGL]{plot}} methods for \code{"cv.FuncompCGL"} object.
#'
#' @details
#' k-fold cross validation.
#'
#' @references
#' Sun, Z., Xu, W., Cong, X., Li G. and Chen K. (2020) \emph{Log-contrast regression with
#' functional compositional predictors: linking preterm infant's gut microbiome trajectories
#' to neurobehavioral outcome}, \href{https://arxiv.org/abs/1808.02403}{https://arxiv.org/abs/1808.02403}
#' \emph{Annals of Applied Statistics}
#'
#' @inherit FuncompCGL author
#'
#'
#' @export


cv.FuncompCGL <- function(y, X, Zc = NULL, lam = NULL, nlam = 100, k = 4:10, ref = NULL,
                          foldid, nfolds = 10, W = rep(1,times = p - length(ref)),
                          trim = 0, outer_maxiter = 1e+6, keep = FALSE, ...) {

  y <- drop(y)
  n <- length(y)

  if (missing(foldid)) {
      foldid <- sample(rep(seq(nfolds), length = n))
  } else {
      nfolds <- length(unique(foldid)) #max(foldid)
  }
  if (nfolds < 3) stop("nfolds must be bigger than 3; nfolds=10 recommended")


  object <- as.list(seq(length(k))) # list of data for different k
  names(object) <- k
  if(!is.null(lam) || length(k) == 1) {
    ## Case I
    if(dim(X)[1] == n ) p = ncol(X) / k else p <- ncol(X) - 2
    for(i in 1:length(k)) {
      object[[i]] <- FuncompCGL(y = y, X = X, Zc = Zc, lam = lam, nlam = nlam, ref = ref,
                                k = k[i], outer_maxiter = outer_maxiter, W = W, ...)
    }
  } else {
    ## Case II:
    ## find commom lambda
    p <- ncol(X) - 2
    this.call <- as.list(seq(length(k)))
    for(i in 1:length(k)) {
    ## Caculate integral Z and W matrix (if W is a functoin)
      object[[i]] <- FuncompCGL(y = y, X = X, Zc = Zc, nlam = 1, ref = ref,
                                k = k[i], outer_maxiter = 0, W = W, ...)
      this.call[[i]] <- object[[i]]$call
    }
    sseq <- object[[1]]$sseq

    lam0 <- max(sapply(object, "[[", "lam")) # Shared lam0 for different k
    for(i in 1:length(k)) {
        object[[i]] <- FuncompCGL(y = y, X = object[[i]]$Z, Zc = Zc,
                                  k = k[i], W = object[[i]]$W, ref = ref,
                                  lam = lam0, nlam = nlam,
                                  outer_maxiter = outer_maxiter, ...)
        object[[i]]$call <- this.call[[i]]
        object[[i]]$sseq <- sseq
    }
  }



  cv.nlam <- sapply(object, function(x) length(drop(x$lam))) # different stopping lambda sequence by dfmax/pfmax
  cv.nlam_id <- which.min(cv.nlam)
  cv.lam <- object[[cv.nlam_id]]$lam # shared lambda sequence for different k
  cv.nlam <- cv.nlam[cv.nlam_id]


  cvm <- matrix(NA, nrow = length(k), ncol = max(cv.nlam))
  rownames(cvm) <- paste0("df=", k)
  colnames(cvm) <- seq(cv.nlam)
  cvsd <- cvm

  if(trim > 0) {
    cvm.trim <- cvm
    cvsd.trim <- cvm
  }

  if(keep) preval <- array(NA, dim = c(n, cv.nlam, length(k)))

  ## Cross-validatoin via different folders
  for(l in 1:length(k)) {

    outlist <- as.list(seq(nfolds))
    for(i in 1:nfolds) {
      which <- foldid == i
      y_train <- y[!which]
      Z <- object[[l]]$Z
      Z_train <- Z[!which, , drop = FALSE]
      Zc_train <- Zc[!which, , drop = FALSE]
      outlist[[i]] <- FuncompCGL(y = y_train, X = Z_train, Zc = Zc_train, ref = ref,
                                 k = k[l], W = object[[l]]$W,
                                 lam = cv.lam, outer_maxiter = outer_maxiter, ...)
    }


    cvstuff <- cv.test(outlist, y, X = cbind(Z, Zc), foldid, lam = cv.lam, trim = trim, keep = keep)
    cvm[l, ] <- cvstuff$cvm
    cvsd[l, ] <- cvstuff$cvsd
    if(keep) preval[, , l] <- cvstuff$fit.preval

    if(trim > 0) {
      cvm.trim[l, ] <- cvstuff$cvmtrim
      cvsd.trim[l, ] <- cvstuff$cvsdtrim
    }

  }


  ## Select lambda
  Ftrim = list(cvm = cvm, cvsd = cvsd)
  Ftrim$cvup = Ftrim$cvm + Ftrim$cvsd
  Ftrim$cvlo = Ftrim$cvm - Ftrim$cvsd
  lammin <- ggetmin(lam = cv.lam, cvm = Ftrim$cvm, cvsd = Ftrim$cvsd, k_list = k)
  Ftrim <- c(Ftrim, lammin)


  if(trim > 0) {
    Ttrim <- list(cvm = cvm.trim, cvsd = cvsd.trim)
    Ttrim$cvup = Ttrim$cvm + Ttrim$cvsd
    Ttrim$cvlo = Ttrim$cvm - Ttrim$cvsd
    lammin <- ggetmin(lam = cv.lam, cvm = Ttrim$cvm, cvsd = Ttrim$cvsd, k_list = k)
    Ttrim <- c(Ttrim, lammin)
  } else {
    Ttrim <- NULL
  }

  result <- list(FuncompCGL.fit = object,
                 lam = cv.lam,
                 Ftrim = Ftrim,
                 Ttrim = Ttrim)
  if(keep) result <- c(result, list(fit.preval = preval, foldid = foldid))
  class(result) <- "cv.FuncompCGL"
  return(result)
}



#' @title
#' Cross-validation for compCL.
#'
#' @description
#' k-fold cross-validation for compCL; produce a plot and return
#' optimal values of \code{lam}.
#'
#' @usage
#' cv.compCL(y, Z, Zc = NULL, intercept = FALSE, lam = NULL,
#'           nfolds = 10, foldid, trim = 0, keep = FALSE, ...)
#'
#'
#' @param Z \code{z} matrix as in \code{compCL}.
#' @param Zc \code{Zc} matrix as in \code{compCL}. Default is \code{NULL}.
#' @param ... other arguments that can be passed to \code{compCL}.
#' @param intercept whether to include an intercept.
#'                  Default is \code{FALSE}.
#' @inheritParams cv.FuncompCGL
#'
#' @return an object of S3 class \code{"cv.compCL"} is returned, which is a list constaining:
#' \item{compCL.fit}{a fitted \code{\link{compCL}} object for the full data.}
#' \item{lam}{the sequence of \code{lam}.}
#' \item{Ftrim}{a list of cross-validation results without trimming:
#'                \itemize{
#'                \item \code{cvm} the mean cross-validated error -  a vector of
#'                      length \code{length(lam).}
#'                \item \code{cvsd} standard error of cvm.
#'                \item \code{cvupper} upper curve = \code{cvm+cvsd}.
#'                \item \code{cvlo} lower curve = \code{cvm-cvsd}.
#'                \item \code{lam.min} the optimal value of \code{lam} that gives minimum cross
#'                            validation error.
#'                \item \code{lam.1se} the largest value of lam such that the error is within 1
#'                      standard error of the minimum \code{cvm}.
#'                }
#'             }
#' \item{Ttrim}{a list of cross-validation result with \code{trim*100\%},
#'              The structure is the same as that for \code{Ftrim}.}
#' \item{foldid}{the values of \code{foldid}.}
#'
#' @details
#' cross-validation and fit full data with selected model.
#'
#' @examples
#' p = 30
#' n = 50
#' beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
#' beta = c(beta, rep(0, times = p - length(beta)))
#' Comp_data = comp_Model(n = n, p = p, beta = beta, intercept = FALSE)
#' cvm1 <- cv.compCL(y = Comp_data$y, Z = Comp_data$X.comp,
#'                   Zc = Comp_data$Zc, intercept = Comp_data$intercept)
#'
#' plot(cvm1)
#' coef(cvm1)
#' ## selection by "lam.min" criterion
#' which(abs(coef(cvm1, s = "lam.min")[1:p]) > 0)
#' ## selection by "lam.1se" criterion
#' which(abs(coef(cvm1, s= "lam.1se")[1:p]) > 0)
#'
#' Comp_data2 = comp_Model(n = 30, p = p, beta = Comp_data$beta, intercept = FALSE)
#' y_hat = predict(cvm1, Znew = Comp_data2$X.comp, Zcnew = Comp_data2$Zc)
#' plot(Comp_data2$y, y_hat,
#'      xlab = "Observed response", ylab = "Predicted response")
#'
#' @seealso
#' \code{\link{compCL}} and \code{\link{cv.compCL}},
#' and \code{\link[=coef.cv.compCL]{coef}},
#' \code{\link[=predict.cv.compCL]{predict}} and
#' \code{\link[=plot.cv.compCL]{plot}} methods for \code{"cv.compCL"} object.
#'
#' @inherit compCL references author
#'
#' @export


cv.compCL <- function(y, Z, Zc = NULL, intercept = FALSE,
                      lam = NULL,
                      nfolds = 10, foldid,
                      trim = 0, keep = FALSE,...) {
  y <- drop(y)
  n <- length(y)
  if (missing(foldid)) foldid <- sample(rep(seq(nfolds), length = n)) else nfolds <- max(foldid)
  if (nfolds < 3) stop("nfolds must be bigger than 3; nfolds=10 recommended")

  compCL.object <- compCL(y = y, Z = Z, Zc = Zc, intercept = intercept, lam = lam, ...)
  cv.lam <- compCL.object$lam # common lambda used

  outlist <- as.list(seq(nfolds))

  ## Now fit the nfold models and store them
  for (i in seq(nfolds)) {
    which <- foldid == i
    # suppress message from transforming into compositional data
    outlist[[i]] <- suppressMessages((compCL(y = y[!which], Z = Z[!which, , drop = FALSE],
                                             Zc = Zc[!which, , drop = FALSE],
                                             intercept = intercept,
                                             lam = cv.lam,  ...)))
  }

  newx <- cbind(compCL.object$Z_log, Zc)
  cvstuff <- cv.test(outlist, y, X = newx, foldid, lam = cv.lam, trim = trim, keep = keep)
  cvm <- drop(cvstuff$cvm)
  cvsd <- drop(cvstuff$cvsd)
  if(trim > 0) {
    cvm.trim <- drop(cvstuff$cvmtrim)
    cvsd.trim <- drop(cvstuff$cvsdtrim)
  }


  Ftrim = list(cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd, cvlo = cvm - cvsd)
  lam.min <- getmin(lam = cv.lam, cvm = cvm, cvsd = cvsd)
  Ftrim <- c(Ftrim, lam.min)

  if(trim > 0) {
    Ttrim <- list(cvm = cvm.trim, cvsd = cvsd.trim, cvupper = cvm.trim + cvsd.trim,
                  cvlo = cvm.trim-cvsd.trim)
    lam.min <- getmin(lam = cv.lam, cvm = cvm.trim, cvsd = cvsd.trim)
    Ttrim <- c(Ttrim, lam.min)
  } else {
    Ttrim <- NULL
  }

  result <- list(compCL.fit = compCL.object,
                 lam = cv.lam,
                 Ftrim = Ftrim,
                 Ttrim = Ttrim,
                 foldid = foldid)


  class(result) <- "cv.compCL"
  return(result)
}