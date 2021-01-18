#' @title
#' Compute information crieteria for the \code{FuncompCGL} model.
#'
#'
#' @description
#' Tune the grid values of the penalty parameter code{lam} and the degrees of freedom of
#' the basis function \code{k} in the \code{FuncompCGL} model by GIC, BIC, or AIC. This
#' function calculates the GIC, BIC, or AIC curve and returns the optimal values of
#' \code{lam} and \code{k}.
#'
#' @usage
#' GIC.FuncompCGL(y, X, Zc = NULL, lam = NULL, nlam = 100, k = 4:10, ref = NULL,
#'               intercept = TRUE, W = rep(1,times = p - length(ref)),
#'               type = c("GIC", "BIC", "AIC"),
#'               mu_ratio = 1.01, outer_maxiter = 1e+6, ...)
#'
#' @param k an integer vector specifying the degrees of freedom of the basis function.
#'
#' @param type a character string specifying which crieterion to use. The choices include
#'             \code{"GIC"} (default), \code{"BIC"}, and \code{"AIC"}.
#'
#' @param outer_maxiter maximum number of loops allowed for the augmented Lanrange method.
#'
#' @param \dots other arguments that could be passed to FuncompCL.
#'
#'
#' @inheritParams FuncompCGL
#'
#' @details
#' The \code{FuncompCGL} model estimation is conducted through minimizing the
#' linearly constrained group lasso criterion
#' \deqn{
#' \frac{1}{2n}\|y - 1_n\beta_0 - Z_c\beta_c - Z\beta\|_2^2 + \lambda \sum_{j=1}^{p} \|\beta_j\|_2,
#' s.t. \sum_{j=1}^{p} \beta_j = 0_k.}
#'
#' The tuning parameters can be selected by the generalized information crieterion (GIC),
#' \deqn{
#' GIC(\lambda,k) = \log{(\hat{\sigma}^2(\lambda,k))} +
#' (s(\lambda, k) - 1)k \log{(max(p*k+p_c+1, n))} \log{(\log{n})}/n
#' ,}
#' where
#' \eqn{\hat{\sigma}^2(\lambda,k) = \|y - 1_n\hat{\beta_0}(\lambda, k) -
#' Z_c\hat{\beta_c}(\lambda, k) - Z\hat{\beta}(\lambda, k) \|_{2}^{2}/n} with \eqn{\hat{\beta_0}(\lambda, k)},
#' \eqn{\hat{\beta_c}(\lambda, k)} and \eqn{\hat{\beta}(\lambda, k)} being the regularized estimators
#' of the regression coefficients, and \eqn{s(\lambda, k)} is the number of nonzero coefficient groups in
#' \eqn{\hat{\beta}(\lambda, k)}.
#'
#' @references
#' Sun, Z., Xu, W., Cong, X., Li G. and Chen K. (2020) \emph{Log-contrast regression with
#' functional compositional predictors: linking preterm infant's gut microbiome trajectories
#' to neurobehavioral outcome}, \href{https://arxiv.org/abs/1808.02403}{https://arxiv.org/abs/1808.02403}
#' \emph{Annals of Applied Statistics}.
#'
#' Fan, Y., and Tang, C. Y. (2013) \emph{Tuning parameter selection in high
#' dimensional penalized likelihood},
#' \href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssb.12001}{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssb.12001}
#' \emph{Journal of the Royal Statistical Society. Series B} \strong{75} 531-552.
#'
#' @return An object of S3 class \code{"GIC.FuncompCGL"} is returned, which is
#' a list containing:
#' \item{FuncompCGL.fit}{a list of length \code{length(k)},
#'                       with fitted \code{\link{FuncompCGL}}
#'                       objects of different degrees
#'                       of freedom of the basis function.}
#'
#' \item{lam}{the sequence of the penalty parameter \code{lam}.}
#'
#' \item{GIC}{a \code{k} by \code{length(lam)} matirx of GIC values.}
#'
#' \item{lam.min}{the optimal values of the degrees of freedom \code{k}
#'                and the penalty parameter \code{lam}.}
#'
#' \item{MSE}{a \code{k} by \code{length(lam)} matirx of mean squared errors.}
#'
#' @seealso
#' \code{\link{FuncompCGL}} and \code{\link{cv.FuncompCGL}},
#' and \code{\link[=predict.GIC.FuncompCGL]{predict}},
#' \code{\link[=coef.GIC.FuncompCGL]{coef}} and
#' \code{\link[=plot.GIC.FuncompCGL]{plot}} methods for \code{"GIC.FuncompCGL"}
#' object.
#'
#' @examples
#' \donttest{
#' df_beta = 5
#' p = 30
#' beta_C_true = matrix(0, nrow = p, ncol = df_beta)
#' beta_C_true[1, ] <- c(-0.5, -0.5, -0.5 , -1, -1)
#' beta_C_true[2, ] <- c(0.8, 0.8,  0.7,  0.6,  0.6)
#' beta_C_true[3, ] <- c(-0.8, -0.8 , 0.4 , 1 , 1)
#' beta_C_true[4, ] <- c(0.5, 0.5, -0.6  ,-0.6, -0.6)
#' n_train = 50
#' n_test = 30
#' k_list <- c(4,5)
#' Data <- Fcomp_Model(n = n_train, p = p, m = 0, intercept = TRUE,
#'                     SNR = 4, sigma = 3, rho_X = 0.2, rho_T = 0.5,
#'                     df_beta = df_beta, n_T = 20, obs_spar = 1, theta.add = FALSE,
#'                     beta_C = as.vector(t(beta_C_true)))
#' arg_list <- as.list(Data$call)[-1]
#' arg_list$n <- n_test
#' Test <- do.call(Fcomp_Model, arg_list)
#'
#' ## GIC_cgl: Constrained group lasso
#' GIC_cgl <- GIC.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
#'                           Zc = Data$data$Zc, intercept = Data$data$intercept,
#'                           k = k_list)
#' coef(GIC_cgl)
#' plot(GIC_cgl)
#' y_hat <- predict(GIC_cgl, Znew = Test$data$Comp, Zcnew = Test$data$Zc)
#' plot(Test$data$y, y_hat, xlab = "Observed response", ylab = "Predicted response")
#'
#' ## GIC_naive: ignoring the zero-sum constraints
#' ## set mu_raio = 0 to identifying without linear constraints,
#' ## no outer_loop for Lagrange augmented multiplier
#' GIC_naive <- GIC.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
#'                             Zc = Data$data$Zc, intercept = Data$data$intercept,
#'                             k = k_list, mu_ratio = 0)
#' coef(GIC_naive)
#' plot(GIC_naive)
#' y_hat <- predict(GIC_naive, Znew = Test$data$Comp, Zcnew = Test$data$Zc)
#' plot(Test$data$y, y_hat, xlab = "Observed response", ylab = "Predicted response")
#'
#' ## GIC_base: random select a component as reference
#' ## mu_ratio is set to 0 automatically once ref is set to a integer
#' ref <- sample(1:p, 1)
#' GIC_base <- GIC.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
#'                             Zc = Data$data$Zc, intercept = Data$data$intercept,
#'                             k = k_list, ref = ref)
#' coef(GIC_base)
#' plot(GIC_base)
#' y_hat <- predict(GIC_base, Znew = Test$data$Comp, Zcnew = Test$data$Zc)
#' plot(Test$data$y, y_hat, xlab = "Observed response", ylab = "Predicted response")
#' }
#'
#' @export


GIC.FuncompCGL <- function(y, X, Zc = NULL, lam = NULL, nlam = 100, k = 4:10, ref = NULL,
                           intercept = TRUE, W = rep(1,times = p - length(ref)),
                           type = c("GIC", "BIC", "AIC"),
                           mu_ratio = 1.01, outer_maxiter = 1e+6, ...) {
  y <- drop(y)
  n <- length(y)
  type <- match.arg(type)
  object <- as.list(seq(length(k)))
  names(object) <- k

  if(!is.null(lam) || length(k) == 1) {
    ## Case I
    if(dim(X)[1] == n ) p = ncol(X) / k else p <- ncol(X) - 2

    for(i in 1:length(k)){
      object[[i]] <- FuncompCGL(y = y, X = X, Zc = Zc, lam = lam, nlam = nlam, ref = ref,
                                k = k[i], outer_maxiter = outer_maxiter, W = W,
                                intercept = intercept, mu_ratio = mu_ratio, ...)
    }
  } else {
    ## Case II
    ## find commom lambda first
    p <- ncol(X) - 2
    this.call <- as.list(seq(length(k)))
    for(i in 1:length(k)){
      ## Caculate integral Z and W matrix (if W is a functoin)
      object[[i]] <- FuncompCGL(y = y, X = X, Zc = Zc, nlam = 1, ref = ref,
                                k = k[i], outer_maxiter = 0, W = W,
                                intercept = intercept, mu_ratio = mu_ratio, ...)
      this.call[[i]] <- object[[i]]$call
    }
    sseq <- object[[1]]$sseq

    lam0 <- max(sapply(object, "[[", "lam")) # Shared lam0 for different k
    ## Solution path for each df k
    for(i in 1:length(k)) {
      object[[i]] <- FuncompCGL(y = y, X = object[[i]]$Z, Zc = Zc,
                                k = k[i], W = object[[i]]$W, ref = ref,
                                lam = lam0, nlam = nlam,
                                outer_maxiter = outer_maxiter,
                                intercept = intercept, mu_ratio = mu_ratio, ...)
      object[[i]]$call <- this.call[[i]]
      object[[i]]$sseq <- sseq
    }
  }

  # lam_start = max of all lam_start
  # lam_end might diff beacuse of stopping criterion
  GIC.nlam <- sapply(object, function(x) length(drop(x$lam)))
  GIC.nlam_id <- which.min(GIC.nlam)
  GIC.lam <- object[[GIC.nlam_id]]$lam
  GIC.nlam <- GIC.nlam[GIC.nlam_id]


  ### >>> get GIC curves <<<
  GIC_curve <- matrix(NA, nrow = length(k), ncol = GIC.nlam)
  MSE <- GIC_curve
  pc = ifelse(is.null(Zc), 0, ncol(Zc)) + as.integer(intercept)
  for(i in seq(length(k))) {
    df = k[i]
    scaler = switch(type,
                    "GIC" = log(max(p*df+ pc, n)) * log(log(n)) / n,
                    "BIC" = log(n) / n,
                    "AIC" = 2 / n
                    )

    predmat <- predict.linear(object = object[[i]], newX = cbind(object[[i]]$Z, Zc))
    MSE[i, ] <- apply(predmat[, 1:GIC.nlam] - y, 2, function(x) mean(x^2))
    ## L: maximuize value of the likelihood fucntion for the estiamted model
    ## Gaussian model: -2*log(L) = n * log(MSE)
    GIC_curve[i, ] <- log(MSE[i, ])
    N_zero <- object[[i]]$df[1:GIC.nlam]


    if(mu_ratio != 0) {
      ## CGL or Baseline model, both with linear constraints
      GIC_curve[i, ] = GIC_curve[i, ] + (ifelse(N_zero == 0, 0, N_zero - 1) * df + pc) * scaler
    } else {
      ## naive model without linear constraints
      GIC_curve[i, ] = GIC_curve[i, ] + (N_zero * df + pc) * scaler
    }
  }
  rownames(GIC_curve) = paste(type, k, sep = "_")
  ## select the optimal grid valud of lambda and k which achieves minimum value of GIC curve
  lam.min <- ggetmin(lam=GIC.lam, cvm = GIC_curve, k_list = k)$lam.min
  result <- list(FuncompCGL.fit = object,
                 lam = GIC.lam,
                 GIC = GIC_curve,
                 lam.min = lam.min,
                 MSE = MSE)
  class(result) <- "GIC.FuncompCGL"
  return(result)
}



#' @title
#' Compute information crieteria for the \code{compCL} model.
#'
#' @description
#' Tune the penalty parameter code{lam} in the \code{compCGL} model by GIC, BIC, or AIC. This
#' function calculates the GIC, BIC, or AIC curve and returns the optimal value of
#' \code{lam}.
#'
#' @inheritParams compCL
#'
#' @param \dots other arguments that can be passed to compCL.
#'
#' @return an object of S3 class \code{GIC.compCL} is returned, which is a list:
#' \item{compCL.fit}{a fitted \code{\link{compCL}} object.}
#' \item{lam}{the sequence of \code{lam}.}
#' \item{GIC}{a vector of GIC value(s).}
#' \item{lam.min}{the \code{lam} value that minimizes \code{GIC}(\eqn{\lambda}).}
#'
#' @details
#' The model estimation is conducted through minimizing the following criterion:
#' \deqn{\frac{1}{2n}\|y-Z\beta\|_2^2 + \lambda\|\beta\|_1, s.t. \sum_{j=1}^{p} \beta_j = 0.}
#' The GIC is defined as:
#' \deqn{GIC(\lambda) = \log{\hat{\sigma}^2(\lambda)} +
#' (s(\lambda) -1) \log{(max(p, n))} * \log{(\log{n})} / n,}
#' where \eqn{\hat{\sigma}^2(\lambda) = \|y - Z\hat{\beta}(\lambda)\|_{2}^{2}/n},
#' \eqn{\hat{\beta}(\lambda)} is the regularized estimator,
#' and \eqn{s(\lambda)} is the number of nonzero coefficients in \eqn{\hat{\beta}(\lambda)}.
#' Because of the zero-sum constraint, the effective number of free parameters is
#' \eqn{s(\lambda) - 1} for \eqn{s(\lambda) \ge 2}.
#' The optimal \eqn{\lambda} is selected by minimizing \code{GIC}(\eqn{\lambda}).
#'
#'
#' @examples
#' p = 30
#' n = 50
#' beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
#' beta = c(beta, rep(0, times = p - length(beta)))
#' Comp_data = comp_Model(n = n, p = p, beta = beta, intercept = FALSE)
#' GICm1 <- GIC.compCL(y = Comp_data$y, Z = Comp_data$X.comp,
#'                     Zc = Comp_data$Zc, intercept = Comp_data$intercept)
#' coef(GICm1)
#' plot(GICm1)
#' test_data = comp_Model(n = 100, p = p, beta = Comp_data$beta, intercept = FALSE)
#' y_hat = predict(GICm1, Znew = test_data$X.comp, Zcnew = test_data$Zc)
#' plot(test_data$y, y_hat, xlab = "Observed value", ylab = "Predicted value")
#' abline(a = 0, b = 1, col = "red")
#'
#'
#' @references
#' Lin, W., Shi, P., Peng, R. and Li, H. (2014) \emph{Variable selection in
#' regression with compositional covariates},
#' \href{https://academic.oup.com/biomet/article/101/4/785/1775476}{https://academic.oup.com/biomet/article/101/4/785/1775476}.
#' \emph{Biometrika} \strong{101} 785-979
#'
#' Fan, Y., and Tang, C. Y. (2013) \emph{Tuning parameter selection in high
#' dimensional penalized likelihood},
#' \href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssb.12001}{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssb.12001}
#' \emph{Journal of the Royal Statistical Society. Series B} \strong{75} 531-552
#'
#' @seealso
#' \code{\link{compCL}} and \code{\link{cv.compCL}},
#' and \code{\link[=coef.GIC.compCL]{coef}}, \code{\link[=predict.GIC.compCL]{predict}} and
#' \code{\link[=plot.GIC.compCL]{plot}} methods for \code{"GIC.compCL"} object.
#'
#' @export


## TODO: add GIC and AIC feature
GIC.compCL <- function(y, Z, Zc = NULL, intercept = FALSE,
                       lam = NULL, ...) {

  this.call <- match.call()

  y <- drop(y)
  n <- length(y)

  p <- ncol(Z)
  pc <- ifelse(is.null(Zc), 0, dim(Zc)[2])
  pc <- pc + as.integer(intercept)


  compCL.object <- compCL(y = y, Z = Z, Zc = Zc, intercept = intercept, lam = lam, ...)
  lam <- compCL.object$lam

  ## >>> MSE <<<
  # predmat <- predict(compCL.object, newx = Z, newzc = Zc)
  newx <- cbind(compCL.object$Z_log, Zc)
  predmat <- predict.linear(compCL.object, newx, s = NULL)
  cvraw <- (y - predmat)^2
  MSE <- colSums(cvraw) / n
  ## <<< MSE >>>

  ## >>> GIC curve <<<
  scaler <- log(log(n)) * log(max(p + pc, n)) / n
  S <- apply(abs(compCL.object$beta[1:p, ]) > 0, 2, sum) # support
  GIC <- log(MSE) + scaler * (ifelse(S>=2, S-1, 0) + pc)
  # digits = 5
  # GIC <- round(GIC, digits = digits)
  ## <<< GIC curve >>>

  # >>> lambda selection <<<
  # GIC.min <- min(GIC[drop(compCL.object$df) > 0])
  GIC.min <- min(GIC)
  idmin <- GIC <= GIC.min
  # idmin[drop(compCL.object$df) < 2 ] <- FALSE
  lam.min <- max(lam[idmin])
  # <<< lambda selection >>>

  result <- list(compCL.fit = compCL.object,
                 lam = lam,
                 GIC = GIC,
                 lam.min = lam.min)

  class(result) <- "GIC.compCL"
  result$call <- this.call
  return(result)
}