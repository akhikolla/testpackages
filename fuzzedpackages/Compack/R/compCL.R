#' @title
#' Fit regularization path for log-contrast model of compositional data with lasso penalty.
#'
#' @description
#' Fit regression with compositional predictors via penalized \emph{log-contrast} model which was proposed by Lin et al. (2014) <doi:10.1093/biomet/asu031>.
#' The model estimation is conducted by minimizing a linearly constrained lasso criterion. The regularization paths are
#' computed at a grid of tuning parameter \code{lambda}.
#'
#'
#'
#'
#' @usage
#' compCL(y, Z, Zc = NULL, intercept = TRUE,
#'        lam = NULL, nlam = 100, lambda.factor = ifelse(n < p, 0.05, 0.001),
#'        pf = rep(1, times = p), dfmax = p, pfmax = min(dfmax * 1.5, p),
#'        u = 1, mu_ratio = 1.01, tol = 1e-10,
#'        inner_maxiter = 1e+4, inner_eps = 1e-6,
#'        outer_maxiter = 1e+08, outer_eps = 1e-8)
#'
#' @param Z a \eqn{n \times p}{n*p} design matrix of compositional data or categorical data.
#'          If \code{Z} is categorical data, i.e., row-sums of \code{Z} differ from 1, the program automatically transforms
#'          \code{Z} into compositional data by dividing each row by its sum.
#'          \code{Z} could NOT include entry of 0's.
#'
#' @param Zc a \eqn{n*p_c} design matrix of control variables (not penalized). Default is \code{NULL}.
#'
#' @param outer_maxiter,outer_eps \code{outer_maxiter} is the maximum number of loops allowed in the Augmented Lagrange method;
#'                                and \code{outer_eps} is the corresponding convergence tolerance.
#'
#' @param inner_maxiter,inner_eps \code{inner_maxiter} is the maximun number of loops allowed in the coordinate descent;
#'                                and \code{inner_eps} is the corresponding convergence tolerance.
#'
#' @param mu_ratio the increasing ratio, with value at least 1, for \code{u}. Default value is 1.01.
#'                 Inital values for scaled Lagrange multipliers are set as 0's.
#'                 If \code{mu_ratio} < 1, the program automatically set \code{u} as 0 and \code{outer_maxiter} as 1, indicating
#'                 that there is no linear constraints included.
#'
#' @param intercept Boolean, specifying whether to include an intercept.
#'                  Default is \code{FALSE}.
#'
#' @inheritParams classo
#'
#' @details
#' The \emph{log-contrast} regression model with compositional predictors is expressed as
#' \deqn{y = Z\beta + e, s.t. \sum_{j=1}^{p}\beta_j=0,}
#' where \eqn{Z} is the n-by-p design matrix of log-transforemd compositional data,
#' \eqn{\beta} is the p-vector of regression cofficients,
#' and \eqn{e} is an n-vector of random errors.
#' If zero(s) exists in the original compositional data, user should pre-process these zero(s).
#' \cr
#' To enable variable selection, we conduct model estimation via linearly constrained lasso
#' \deqn{
#' argmin_{\beta}(\frac{1}{2n}\|y-Z\beta\|_2^2 + \lambda\|\beta\|_1), s.t. \sum_{j=1}^{p}\beta_j= 0.
#' }
#'
#' @return An object with S3 calss \code{"compCL"} is a list containing:
#' \item{beta}{a matrix of coefficients for \eqn{p+p_c+1} rows.
#'             If \code{intercept=FALSE}, then the last row of \code{beta} is set to 0's.}
#' \item{lam}{the sequence of \code{lam} values used.}
#' \item{df}{the number of non-zero \eqn{\beta_p}'s in estimated coefficients for \code{Z} at each value of \code{lam}.}
#' \item{npass}{total iterations.}
#' \item{error}{error messages. If 0, no error occurs.}
#' \item{call}{the call that produces this object.}
#' \item{dim}{dimension of the coefficient matrix \code{beta}.}
#'
#' @author
#' Zhe Sun and Kun Chen
#'
#' @references
#' Lin, W., Shi, P., Peng, R. and Li, H. (2014) \emph{Variable selection in regression with compositional covariates},
#' \href{https://academic.oup.com/biomet/article/101/4/785/1775476}{https://academic.oup.com/biomet/article/101/4/785/1775476}.
#' \emph{Biometrika} \strong{101} 785-979
#'
#' @seealso
#' \code{\link[=coef.compCL]{coef}}, \code{\link[=predict.compCL]{predict}},
#' \code{\link[=print.compCL]{print}} and \code{\link[=plot.compCL]{plot}} methods
#' for \code{"compCL"} object
#' and \code{\link{cv.compCL}} and \code{\link{GIC.compCL}}.
#'
#' @examples
#' p = 30
#' n = 50
#' beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
#' beta = c(beta, rep(0, times = p - length(beta)))
#' Comp_data = comp_Model(n = n, p = p, beta = beta, intercept = FALSE)
#' m1 <- compCL(y = Comp_data$y, Z = Comp_data$X.comp,
#'              Zc = Comp_data$Zc, intercept = Comp_data$intercept)
#' print(m1)
#' plot(m1)
#' beta = coef(m1)
#' Test_data = comp_Model(n = 30, p = p, beta = Comp_data$beta, intercept = FALSE)
#' predmat = predict(m1, Znew = Test_data$X.comp, Zcnew = Test_data$Zc)
#'
#' @export

## TODO: add feature for baseline method

compCL <- function(y, Z, Zc = NULL, intercept = TRUE,
                   lam = NULL, nlam = 100, lambda.factor = ifelse(n < p, 0.05, 0.001),
                   pf = rep(1, times = p), dfmax = p, pfmax = min(dfmax * 1.5, p),
                   u = 1, mu_ratio = 1.01, tol =1e-10 ,
                   inner_maxiter = 1e+4, inner_eps = 1e-6,
                   outer_maxiter = 1e+08, outer_eps = 1e-8) {

  if(!is.null(lam) && TRUE %in% (lam < 0)) stop("User provided lambda must be positive vector")

  this.call <- match.call()

  y <- drop(y)
  n <- length(y)
  Z <- proc.comp(Z)
  p <- ncol(Z)
  Znames <- colnames(Z)
  if (is.null(Znames)) Znames <- paste0("Z", seq(p))

  inter <- as.integer(intercept)

  # Zc <- as.matrix(cbind(Zc, rep(inter, n)))
  if(is.null(Zc)) Zc <- matrix(inter, nrow = n) else Zc <- cbind2(as.matrix(Zc), inter)

  if(inter == 1) {
    Zc_proj <- solve(crossprod(Zc)) %*% t(Zc)
  } else {
    if(dim(Zc)[2] == 1) Zc_proj <- t(Zc) else Zc_proj <- rbind(solve(crossprod(Zc[, -ncol(Zc)])) %*% t(Zc[, -ncol(Zc)]), rep(inter, n))
  }

  beta_c <- as.vector(Zc_proj %*% y)
  beta.ini <- c(rep(0, times = p), beta_c)

  m <- ncol(Zc)
  if(m > 1) {
    Zcnames <- colnames(Zc)[-ncol(Zc)]
    if (is.null(Zcnames)) Zcnames <- paste0("Zc", seq(m-1))
  } else {
    Zcnames <- NULL
  }

  if(is.null(lam)) {
    lam0 <- t(Z) %*% (y - Zc %*% beta_c) / n
    lam0 <- max(abs(lam0) / pf)
    lam <- exp(seq(from=log(lam0), to=log(lambda.factor*lam0), length=nlam))
  } else {
    if(length(lam) == 1 && nlam > 1) {
      lam <- exp(seq(from=log(lam), to=log(lambda.factor*lam), length=nlam))
    } else {
      lam <- sort(lam, decreasing = TRUE)
    }
  }

  pfmax <- as.integer(pfmax)
  dfmax <- as.integer(dfmax)
  inner_maxiter <- as.integer(inner_maxiter)
  outer_maxiter <- as.integer(outer_maxiter)
  inner_eps <- as.double(inner_eps)
  outer_eps <- as.double(outer_eps)
  tol <- as.double(tol)
  u <- as.double(u)
  mu_ratio <- as.double(mu_ratio)

  output <- ALM_CD_comp(y = y, Z = Z, Zc = Zc, Zc_proj = Zc_proj, beta = beta.ini,
                        lambda = lam, pf = pf, dfmax = dfmax, pfmax = pfmax,
                        inner_eps = inner_eps, outer_eps = outer_eps,
                        inner_maxiter = inner_maxiter, outer_maxiter = outer_maxiter,
                        u_ini = u, mu_ratio = mu_ratio, tol = tol

  )

  output$call <- this.call
  class(output) <- "compCL"
  output$dim <- dim(output$beta)
  output$lam <- drop(output$lam)
  output$Z_log <- Z
  dimnames(output$beta) <- list(c(Znames, Zcnames, "Intercept"), paste0("L", seq(along = output$lam)) )
  return(output)

}



