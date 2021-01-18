######################################################################
## These functions are minor modifications or directly copied from the
## glmnet package:
##        Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
## Regularization Paths for Generalized Linear Models via Coordinate
#   Descent.
##        Journal of Statistical Software, 33(1), 1-22.
##        URL http://www.jstatsoft.org/v33/i01/.
## The reason they are copied here is because they are internal functions
## and hence are not exported into the global environment.
#######################################################################

#######################################################################
## roxygen2
## references and author sections are passed from coef.compCL and coef.FuncompCL
#######################################################################

#### >>> coef method <<< ####

#' @title
#' Extract estimated coefficients from a \code{"FuncompCGL"} object.
#'
#' @description
#' get the coefficients at the requested values for \code{lam}
#' from a fitted \code{\link{FuncompCGL}} object.
#'
#' @param object fitted \code{\link{FuncompCGL}} object.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which coefficients
#'          are requested. Default (\code{NULL}) is the entire sequence used to
#'          create the model.
#'
#' @param \dots Not used.
#'
#' @details
#' \code{s} is a vector of lambda values at which the coefficients are requested. If \code{s} is not in the
#' \code{lam} sequence used for fitting the model, the \code{coef} function will use linear
#' interpolation, so the function should be used with caution.
#'
#' @return
#' The coefficients at the requested tuning parameter values in \code{s}.
#'
#' @author
#' Zhe Sun and Kun Chen
#'
#' @references
#' Sun, Z., Xu, W., Cong, X., Li G. and Chen K. (2020) \emph{Log-contrast regression with
#' functional compositional predictors: linking preterm infant's gut microbiome trajectories
#' to neurobehavioral outcome}, \href{https://arxiv.org/abs/1808.02403}{https://arxiv.org/abs/1808.02403}
#' \emph{Annals of Applied Statistics}
#'
#' @seealso
#' \code{\link{FuncompCGL}}, and \code{\link[=predict.FuncompCGL]{predict}},
#' \code{\link[=plot.FuncompCGL]{plot}} and \code{\link[=print.FuncompCGL]{print}}
#' methods for \code{"FuncompCGL"} object.
#'
#' @examples
#' df_beta = 5
#' p = 30
#' beta_C_true = matrix(0, nrow = p, ncol = df_beta)
#' beta_C_true[1, ] <- c(-0.5, -0.5, -0.5 , -1, -1)
#' beta_C_true[2, ] <- c(0.8, 0.8,  0.7,  0.6,  0.6)
#' beta_C_true[3, ] <- c(-0.8, -0.8 , 0.4 , 1 , 1)
#' beta_C_true[4, ] <- c(0.5, 0.5, -0.6  ,-0.6, -0.6)
#' Data <- Fcomp_Model(n = 50, p = p, m = 2, intercept = TRUE,
#'                     SNR = 2, sigma = 2, rho_X = 0, rho_T = 0.5, df_beta = df_beta,
#'                     n_T = 20, obs_spar = 1, theta.add = FALSE,
#'                     beta_C = as.vector(t(beta_C_true)))
#'
#' m1 <- FuncompCGL(y = Data$data$y, X = Data$data$Comp , Zc = Data$data$Zc,
#'                  intercept = Data$data$intercept, k = df_beta)
#' coef(m1)
#' coef(m1, s = c(0.5, 0.1, 0.01))
#'
#' @export


coef.FuncompCGL <- function(object, s = NULL, ...) {
  beta <- object$beta
  if (!is.null(s)) {
    lam <- object$lam
    lamlist <- point.interp(lam, s)
    if(length(s) == 1)
    {
      beta = beta[, lamlist$left, drop = FALSE] *  (1 - lamlist$frac) +
        beta[, lamlist$right, drop = FALSE] * lamlist$frac
    } else {
      beta = beta[, lamlist$left, drop = FALSE] %*% diag(1 - lamlist$frac) +
        beta[, lamlist$right, drop = FALSE] %*% diag(lamlist$frac)
    }
    rownames(seq(s))
  }

  if(is.numeric(object$ref)) {
    beta <- apply(beta, 2, function(x, k, p1, ref)
      c(x[ifelse(ref==1, 0, 1):((ref-1)*k)], -colSums(matrix(x[1:p1], byrow = TRUE, ncol = k)), x[((ref-1)*k+1):(length(x))]),
      ref = object$ref, k = as.list(object$call)$k, p1 = ncol(object$Z) )
  }

  return(beta)
}


#' @title
#' Extract estiamted coefficients from a \code{"cv.FuncompCGL"} object.
#'
#' @description
#' This function gets the coefficients from a cross-validated
#' \code{FuncompCGL} model, using the stored \code{"FuncompCGL.fit"} object,
#' and the optimal grid values of the penalty parameter \code{lam} and the degrees of freedom \code{k}.
#'
#' @param object fitted \code{\link{cv.FuncompCGL}} object.
#'
#' @param trim logical; whether to use the trimmed result. Default is \code{FALSE}.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which coefficients are requested.
#'          \itemize{
#'          \item \code{s="lam.min"}(default), grid value of
#'          \code{lam} and \code{k} stored in the \code{"cv.FuncompCGL"} object
#'          such that the minimum cross-validation error is achieved.
#'          \item \code{s="lam.1se"}, grid value of
#'          \code{lam} and \code{k} stored on the \code{"cv.FuncompCGL"} object
#'          such that the 1 standard error above the miminum cross-validation error is achieved.
#'          \item If \code{s} is numeric, it is taken as the value(s) of \code{lam} to be used. In this
#'          case, \code{k} must be provided.
#'          \item If \code{s = NULL}, the whole sequence of \code{lam} stored in the \code{cv.FuncompCGL}
#'                object is used.
#'          }
#'
#' @param k value(s) of the degrees of freedom of the basis function at which coefficents are requested.
#'          \code{k} can be \code{NULL} (default) or integer(s).
#'          \itemize{
#'          \item \code{k = NULL}, \code{s} must be either \code{"lam.min"} or \code{"lam.1se"}.
#'          \item if \code{k} is an integer(s), it is taken as the value of \code{k} to be used and
#'                it must be one(s) of these in the \code{"cv.FuncompCGL"} object.
#'          }
#'
#' @param \dots not used.
#'
#' @return
#' The coefficients at the requested values of \code{s} and \code{k}.
#' If \code{k} is a vector, a list of coefficient matrices is returned.
#'
#' @seealso
#' \code{\link{cv.FuncompCGL}} and \code{\link{FuncompCGL}}, and
#' \code{\link[=predict.cv.FuncompCGL]{predict}} and
#' \code{\link[=plot.cv.FuncompCGL]{plot}} methods for \code{"cv.FuncompCGL"} object.
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
#' Data <- Fcomp_Model(n = 50, p = p, m = 0, intercept = TRUE,
#'                     SNR = 4, sigma = 3, rho_X = 0, rho_T = 0.6, df_beta = df_beta,
#'                     n_T = 20, obs_spar = 1, theta.add = FALSE,
#'                     beta_C = as.vector(t(beta_C_true)))
#'
#' cv_m1 <-  cv.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
#'                         Zc = Data$data$Zc, intercept = Data$data$intercept,
#'                         k = c(4,5), nfolds = 5, nlam = 50,
#'                         keep = TRUE)
#' coef(cv_m1)
#' coef(cv_m1, s = "lam.1se")
#' coef(cv_m1, s = c(0.5, 0.1, 0.05), k = c(4,5))
#' coef(cv_m1, s = NULL, k = c(4,5))
#' }
#'
#' @inherit coef.FuncompCGL details references author
#'
#' @export

coef.cv.FuncompCGL <- function(object, trim = FALSE, s = c("lam.min", "lam.1se"), k = NULL, ...) {
  trim <- ifelse(trim, "Ttrim", "Ftrim")

  if (is.numeric(s) || is.null(s) ) {
    lam_se <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    lam_se <- object[[trim]][[s]]["lam"]
    k_se <- object[[trim]][[s]]["df"]
  } else stop("Invalid form for s")

  if(is.numeric(k)){
    if(FALSE %in% (k %in% as.integer(names(object$FuncompCGL.fit)))) stop("K is out of range")
    k_se <- k
  }
  k_se <- as.character(k_se)
  if(length(k_se) > 1) {
    beta <- as.list(k_se)
    for(i in seq(length(k_se))) {
      cv.fit <- object$FuncompCGL.fit[[k_se[i]]]
      beta[[i]] <- coef(cv.fit, s = lam_se)
    }
    names(beta) <- paste("k", k_se, sep = "=")
  } else {
    cv.fit <- object$FuncompCGL.fit[[k_se]]
    beta <- coef(cv.fit, s = lam_se)
  }

  return(beta)
}



#' @title
#' Extract model estimated coefficients from a \code{"GIC.FuncompCGL"} object.
#'
#' @description
#' This function gets coefficients from a \code{"GIC.FuncompCGL"} object,
#' using the stored \code{"FuncompCGL.fit"} object, and the optimal values of
#' \code{lam} and \code{k}.
#'
#'
#' @param object fitted \code{\link{GIC.FuncompCGL}} object.
#'
#'
#' @param s value(s) of the regularization parameter \code{lam} at which coefficients are requested.
#'          \itemize{
#'          \item \code{s="lam.min"} (default), grid value of \code{lam} and \code{k} stored in
#'                \code{"GIC.FuncompCGL"} object such that the minimun value of GIC is achieved.
#'          \item If \code{s} is numeric, it is taken as the value(s) of \code{lam} to be used.
#'                In this case, k must be provided.
#'          \item If \code{s = NULL}, used the whole sequence of \code{lam} stored in the \code{GIC.FuncompCGL}
#'                object.
#'          }
#'
#' @param k value(s) of degrees of freedom of the basis function at which coefficents are requested.
#'          \code{k} can be \code{NULL} (default) or integer(s).
#'          \itemize{
#'          \item \code{k = NULL}, \code{s} must be \code{"lam.min"}.
#'          \item if \code{k} is integer(s), it is taken as the value of \code{k} to be used
#'                and it must be one(s) of these in \code{"GIC.FuncompCGL"} model.
#'          }
#'
#' @param \dots not used.
#'
#' @seealso
#' \code{\link{GIC.FuncompCGL}} and \code{\link{FuncompCGL}}, and
#' \code{\link[=predict.GIC.FuncompCGL]{predict}} and
#' \code{\link[=plot.GIC.FuncompCGL]{plot}} methods for \code{"GIC.FuncompCGL"} object.
#'
#' @examples
#' df_beta = 5
#' p = 30
#' beta_C_true = matrix(0, nrow = p, ncol = df_beta)
#' beta_C_true[1, ] <- c(-0.5, -0.5, -0.5 , -1, -1)
#' beta_C_true[2, ] <- c(0.8, 0.8,  0.7,  0.6,  0.6)
#' beta_C_true[3, ] <- c(-0.8, -0.8 , 0.4 , 1 , 1)
#' beta_C_true[4, ] <- c(0.5, 0.5, -0.6  ,-0.6, -0.6)
#' n = 50
#' k_list <- c(4,5)
#' Data <- Fcomp_Model(n = n, p = p, m = 0, intercept = TRUE,
#'                     SNR = 4, sigma = 3, rho_X = 0.6, rho_T = 0,
#'                     df_beta = df_beta, n_T = 20, obs_spar = 1, theta.add = FALSE,
#'                     beta_C = as.vector(t(beta_C_true)))
#'
#' GIC_m1 <-  GIC.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
#'                           Zc = Data$data$Zc, intercept = Data$data$intercept,
#'                           k = k_list)
#' coef(GIC_m1)
#' coef(GIC_m1, s = c(0.05, 0.01), k = c(4,5))
#' coef(GIC_m1, s = NULL, k = c(4,5))
#'
#' @inherit coef.FuncompCGL details return references author
#'
#' @export

coef.GIC.FuncompCGL <- function(object, s ="lam.min", k = NULL, ...) {

  if (is.numeric(s) || is.null(s) ) {
    lam_se <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    lam_se <- object[[s]]['lam']
    k_se <- object[[s]]['df']
  } else stop("Invalid form for s")

  if(is.numeric(k)){
    if(FALSE %in% (k %in% as.integer(names(object$FuncompCGL.fit)))) stop("K is out of range")
    k_se <- k
  }
  k_se <- as.character(k_se)
  if(length(k_se) > 1) {
    beta <- as.list(k_se)
    for(i in seq(length(k_se))) {
      GIC.fit <- object$FuncompCGL.fit[[k_se[i]]]
      beta[[i]] <- coef(GIC.fit, s = lam_se)
    }
    names(beta) <- paste("k", k_se, sep = "=")
  } else {
    GIC.fit <- object$FuncompCGL.fit[[k_se]]
    beta <- coef(GIC.fit, s = lam_se)
  }
  return(beta)
}



#' @title
#' extracts model estimated coefficients from a \code{"compCL"} object.
#'
#' @description
#' gets the coefficients at the requested values for \code{lam} from a fitted
#' \code{"\link{compCL}"} object.
#'
#' @param object fitted \code{"\link{compCL}"} object.
#'
#' @inheritParams coef.FuncompCGL
#'
#' @author
#' Zhe Sun and Kun Chen
#'
#' @references
#' Lin, W., Shi, P., Peng, R. and Li, H. (2014) \emph{Variable selection in
#' regression with compositional covariates},
#' \href{https://academic.oup.com/biomet/article/101/4/785/1775476}{https://academic.oup.com/biomet/article/101/4/785/1775476}.
#' \emph{Biometrika} \strong{101} 785-979.
#'
#' @seealso
#' \code{\link{compCL}} and \code{\link[=predict.compCL]{predict}},
#' \code{\link[=plot.compCL]{plot}} and \code{\link[=print.compCL]{print}} methods
#' for \code{"compCL"} object.
#'
#' @examples
#' Comp_data = comp_Model(n = 50, p = 30)
#' Comp_fit = compCL(y = Comp_data$y, Z = Comp_data$X.comp, Zc = Comp_data$Zc,
#'                   intercept = Comp_data$intercept)
#' coef(Comp_fit)
#' coef(Comp_fit, s = Comp_fit$lam[50])
#' coef(Comp_fit, s = c(1, 0.5, 0.1))
#'
#' @inherit coef.FuncompCGL details return
#'
#' @export


## TODO: adjust after adding baseline feature in compCL

coef.compCL <- function(object, s = NULL, ...) {
  beta <- object$beta
  if (!is.null(s)) {
    lam <- object$lam
    lamlist <- point.interp(lam, s)
    if(length(s) == 1) {
      beta = beta[, lamlist$left, drop=FALSE] *  (1 - lamlist$frac) +
        beta[, lamlist$right, drop=FALSE] * lamlist$frac
    } else {
      beta = beta[, lamlist$left, drop=FALSE] %*% diag(1 - lamlist$frac) +
        beta[, lamlist$right, drop=FALSE] %*% diag(lamlist$frac)
    }
    colnames(beta) <- paste(seq(along = s))
  }
  return(beta)
}



#' @title
#' Extract estimated coefficients from a \code{"cv.compCL"} object.
#'
#' @description
#' This function gets coefficients from a \code{compCL} object, using the
#' stored \code{"compCL.fit"} object.
#'
#' @param object fitted \code{"\link{cv.compCL}"} object.
#'
#' @param trim whether to use the trimmed result. Default is FASLE.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which coefficients are requested.
#'          \itemize{
#'          \item \code{s="lam.min"} (default) stored in the \code{cv.compCL} object,
#'                which gives value of \code{lam} that achieves the minimum
#'                cross-vadilation error.
#'          \item \code{s="lambda.min"} which gives the largest value of \code{lam}
#'                 such that 1 standard error above the minimum of the cross-validation errors is achieved.
#'          \item If \code{s} is numeric, it is taken as the value(s) of \code{lam} to
#'                be used.
#'          \item If \code{s = NULL}, the whole sequence of \code{lam} stored in the
#'                \code{cv.compCGL} object is used.
#'          }
#'
#' @param \dots not used.
#'
#' @examples
#' p = 30
#' n = 50
#' beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
#' beta = c(beta, rep(0, times = p - length(beta)))
#' Comp_data = comp_Model(n = n, p = p, beta = beta, intercept = FALSE)
#' cvm1 <- cv.compCL(y = Comp_data$y, Z = Comp_data$X.comp,
#'                   Zc = Comp_data$Zc, intercept = Comp_data$intercept)
#' coef(cvm1)
#' coef(cvm1, s = NULL)
#' coef(cvm1, s = c(1, 0.5, 0.1))
#'
#' @seealso
#' \code{\link{cv.compCL}} and \code{\link{compCL}}, and
#' \code{\link[=predict.cv.compCL]{predict}} and \code{\link[=plot.cv.compCL]{plot}} methods
#' for \code{"cv.compCL"} object.
#'
#' @inherit coef.compCL details return author references
#'
#' @export

coef.cv.compCL <- function(object, trim = FALSE, s = c("lam.min", "lam.1se" ),...) {
  trim <- ifelse(trim, "Ttrim", "Ftrim")
  object_use <- object[[trim]]
  if (is.numeric(s) || is.null(s) ) {
    lam <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    lam <- object_use[[s]]
  } else {
    stop("Invalid form for s")
  }
  beta = coef(object$compCL.fit, s = lam, ...)
  return(beta)
}



#' @title
#' Extracts estimated coefficients from a \code{"GIC.compCL"} object.
#'
#' @description
#' This function gets coefficients from a \code{compCL} object, using
#' the stored \code{"compCL.fit"} object.
#'
#' @param object fitted \code{"\link{GIC.compCL}"} object.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which coefficients are requested.
#'          \itemize{
#'          \item \code{s="lam.min"} (default) stored in the \code{GIC.compCL} object,
#'                which gives value of \code{lam} that achieves the minimum value of GIC.
#'          \item If \code{s} is numeric, it is taken as the value(s) of \code{lam} to be used.
#'          \item If \code{s = NULL}, the whole sequence of \code{lam} stored
#'                in the \code{GIC.compCGL} object is used.
#'          }
#'
#' @param \dots not used.
#'
#' @seealso
#' \code{\link{GIC.compCL}} and \code{\link{compCL}}, and
#' \code{\link[=predict.GIC.compCL]{predict}}, and
#' \code{\link[=plot.GIC.compCL]{plot}} methods for \code{"GIC.compCL"} object.
#'
#' @examples
#' p = 30
#' n = 50
#' beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
#' beta = c(beta, rep(0, times = p - length(beta)))
#' Comp_data = comp_Model(n = n, p = p, beta = beta, intercept = FALSE)
#' GICm1 <- GIC.compCL(y = Comp_data$y, Z = Comp_data$X.comp, Zc = Comp_data$Zc,
#'                     intercept = Comp_data$intercept)
#' coef(GICm1)
#' coef(GICm1, s = c(1, 0.5, 0.1))
#'
#' @inherit coef.compCL details return author references
#'
#' @export


coef.GIC.compCL <- function(object, s = "lam.min",...){
  if (is.numeric(s) || is.null(s) ) {
    lam <- s
  } else if (s == "lam.min") {
    lam <- object[[s]]
  } else {
    stop("Invalid form for s")
  }
  beta = beta = coef(object$compCL.fit, s = lam, ...)
  return(beta)
}

#### <<< coef method >>> ####




#### >>> print method <<< ####

#' @title
#' Print a \code{"FuncompCGL"} object.
#'
#' @description
#' print the number of nonzero coefficient curves for the functional compositional
#' predictors at each \code{lam} along the FuncompCGL path.
#'
#' @param x fitted \code{\link{FuncompCGL}} object.
#'
#' @param digits significant digits in printout.
#'
#' @param \dots not used.
#'
#' @return
#' a two-column matrix; the first column \code{DF} gives the number of nonzero
#' coefficients and
#' the second column \code{Lam} gives the correspondint \code{lam} values.
#'
#' @seealso
#' \code{\link{FuncompCGL}}, and \code{\link[=coef.FuncompCGL]{coef}},
#' \code{\link[=predict.FuncompCGL]{predict}} and
#' \code{\link[=plot.FuncompCGL]{plot}} methods for \code{"FuncompCGL"} object.
#'
#' @examples
#' df_beta = 5
#' p = 30
#' beta_C_true = matrix(0, nrow = p, ncol = df_beta)
#' beta_C_true[1, ] <- c(-0.5, -0.5, -0.5 , -1, -1)
#' beta_C_true[2, ] <- c(0.8, 0.8,  0.7,  0.6,  0.6)
#' beta_C_true[3, ] <- c(-0.8, -0.8 , 0.4 , 1 , 1)
#' beta_C_true[4, ] <- c(0.5, 0.5, -0.6  ,-0.6, -0.6)
#' Data <- Fcomp_Model(n = 50, p = p, m = 2, intercept = TRUE,
#'                     SNR = 2, sigma = 2, rho_X = 0, rho_T = 0.5, df_beta = df_beta,
#'                     n_T = 20, obs_spar = 1, theta.add = FALSE,
#'                     beta_C = as.vector(t(beta_C_true)))
#' m1 <- FuncompCGL(y = Data$data$y, X = Data$data$Comp ,
#'                  Zc = Data$data$Zc, intercept = Data$data$intercept,
#'                  k = df_beta)
#' print(m1)
#'
#' @inherit coef.FuncompCGL references author
#'
#' @export

print.FuncompCGL <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  show <- cbind(Df = x$df, Lam = signif(x$lam, digits))
  colnames(show) <- c("DF", "Lam")
  rownames(show) <- seq(nrow(show))
  print(show)
}

#' @title
#' Print a \code{"compCL"} object.
#'
#' @description
#' print the number of nonzero coefficients for the compositional varaibles at each step along the compCL path.
#'
#' @param x fitted \code{"\link{compCL}"} object.
#'
#' @inheritParams print.FuncompCGL
#'
#' @return
#' a two-column matrix; the first column \code{DF} gives the number of nonzero
#' coefficients for the compositional predictors and
#' the second column \code{Lam} gives the corresponding \code{lam} values.
#'
#' @seealso
#' \code{\link{compCL}} and \code{\link[=coef.compCL]{coef}},
#' \code{\link[=predict.compCL]{predict}} and \code{\link[=plot.compCL]{plot}} methods
#' for \code{"compCL"} object.
#'
#' @examples
#' Comp_data = comp_Model(n = 50, p = 30)
#' Comp_fit = compCL(y = Comp_data$y, Z = Comp_data$X.comp,
#'                   Zc = Comp_data$Zc, intercept = Comp_data$intercept)
#' print(Comp_fit)
#'
#' @inherit coef.compCL references author
#'
#' @export


print.compCL <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  show <- cbind(Df = x$df, Lam = signif(x$lam, digits))
  colnames(show) <- c("Df", "Lam")
  rownames(show) <- paste0("L", seq(nrow(show)))
  print(show)
  #print(cbind(Df = object$df, Lam = signif(object$lam, digits)))
}
#### <<< print method >>> ####




#### >>> predict method <<< ####

#' @title
#' Make prediction from a \code{"FuncompCGL"} object.
#'
#' @description
#' Make prediction based on a fitted \code{\link{FuncompCGL}} object.
#'
#' @usage
#' \method{predict}{FuncompCGL}(object, Znew, Zcnew = NULL, s = NULL,
#'         T.name = "TIME", ID.name = "Subject_ID",
#'         Trange, interval, insert, basis_fun, degree, method, sseq,
#'         ...)
#'
#' @param object fitted \code{\link{FuncompCGL}} object.
#'
#' @param Znew data frame or matrix \code{X} as in \code{FuncompCGL} with new
#'             functional compositional data at which prediction is to be made.
#'
#' @param Zcnew matrix \code{Zc} as in \code{FuncompCGL} with new values of
#'              time-invariate covariates at which prediction is to be made.
#'              Default is \code{NULL}.
#'
#' @param s value(s) of the penalty parameter \code{lam} at which predictions are requested.
#'          Default is the entire sequence used to fit the model.
#'
#' @param Trange,interval,insert,basis_fun,degree,method the same as those in \code{FuncompCGL}.
#'
#' @param sseq full set of potential time points of observations;
#'             used for interpolation when \code{insert = "X"}
#'             or \code{insert = "basis"}.
#'
#' @param \dots not used.
#'
#' @inheritParams FuncompCGL
#'
#' @details
#' \code{s} is the vector at which predictions are requested. If \code{s} is not in the \code{lam}
#' sequence used for fitting the model, the \code{predict} function uses linear interpolation.
#' \cr
#' If the data frame \code{X} is provided in \code{FuncompCGL} mode, the integral
#' for new data \code{newx} is taken the same as that in the fitted
#' \code{FuncompCGL} model. This means that the parameters \code{degree},
#' \code{basis_fun}, \code{insert}, \code{method}, \code{inteval},
#' \code{Trange}, and \code{K} are exactly the same as these in the provided
#' \code{object}. If \code{insert="X"} or \code{"basis"}, \code{sseq} is the
#' sorted sequence of all the observed time points in fitting \code{FuncompCGL} model and
#' all the observed time points in \code{newx}. Then interpolation is
#' conducted on \code{sseq}. If matrix \code{X} after integral is provided in
#' the \code{FuncompCGL} object, these parameters are required.
#'
#' @return
#' predicted values at the requested value(s) for \code{s}.
#'
#' @seealso
#' \code{\link{FuncompCGL}}, and \code{\link[=coef.FuncompCGL]{coef}},
#' \code{\link[=plot.FuncompCGL]{plot}} and \code{\link[=print.FuncompCGL]{print}}
#' methods for \code{"FuncompCGL"} object.
#'
#' @examples
#' p = 30
#' n_train = 50
#' n_test = 30
#' df_beta = 5
#' beta_C_true = matrix(0, nrow = p, ncol = df_beta)
#' beta_C_true[1, ] <- c(-0.5, -0.5, -0.5 , -1, -1)
#' beta_C_true[2, ] <- c(0.8, 0.8,  0.7,  0.6,  0.6)
#' beta_C_true[3, ] <- c(-0.8, -0.8 , 0.4 , 1 , 1)
#' beta_C_true[4, ] <- c(0.5, 0.5, -0.6  ,-0.6, -0.6)
#' Data <- Fcomp_Model(n = n_train, p = p, m = 0, intercept = TRUE,
#'                     SNR = 2, sigma = 2,
#'                     rho_X = 0, rho_T = 0.5, df_beta = df_beta,
#'                     n_T = 20, obs_spar = 1, theta.add = c(3,4,5),
#'                     beta_C = as.vector(t(beta_C_true)))
#' m1 <- FuncompCGL(y = Data$data$y, X = Data$data$Comp , Zc = Data$data$Zc,
#'                  intercept = Data$data$intercept, k = df_beta)
#' arg_list <- as.list(Data$call)[-1]
#' arg_list$n <- n_test
#' TEST <- do.call(Fcomp_Model, arg_list)
#' predmat <- predict(m1, Znew = TEST$data$Comp, Zcnew = TEST$data$Zc)
#' predmat <- predict(m1, Znew = TEST$data$Comp, Zcnew = TEST$data$Zc, s = c(0.5, 0.1, 0.05))
#'
#' @inherit coef.FuncompCGL references author
#'
#' @export
#'
#' @importFrom methods cbind2


predict.FuncompCGL <- function(object, Znew, Zcnew = NULL, s = NULL,
                               T.name = "TIME", ID.name = "Subject_ID",
                               Trange, interval, insert, basis_fun, degree, method, sseq,
                               ...)
{
  newx = Znew
  newZc = Zcnew

  beta = coef(object)
  if(is.vector(newZc)){
    newZc <- matrix(newZc, nrow = 1)
  } else if (is.data.frame(newZc)) {
    newZc <- as.matrix(newZc)
  }
  pc = ifelse(is.null(dim(newZc)), 0, ncol(newZc))
  p1 = nrow(beta)

  if(ncol(newx) == (p1 - pc - 1)){
    ## case I: already take integral and ready to go
    Z <- as.matrix(newx)
    p1 <- ncol(Z)
  } else {
    ## case II: take integral first
    digits = 10
    list_param = as.list(object$call)

    newx = as.data.frame(newx)
    if( !(T.name %in% names(newx)) || !(ID.name %in% names(newx))) stop("Please provide ID.name and T.name")
    newx[, T.name] = round(newx[, T.name], digits = digits)
    newx[, ID.name] <- factor(newx[, ID.name], levels = unique(newx[, ID.name])) # In case that newx is not represented in numerical order of Subject_ID
    X.names <- colnames(newx)[!colnames(newx) %in% c(T.name, ID.name)]
    newx[, X.names] = proc.comp(newx[, X.names])
    p = length(X.names)
    k_choice = (p1 - 1 - pc) / p
    if(k_choice %% 1 != 0) stop("fitted data and predicting data with different dimensions")

    # >>> integral domain <<<
    Trange_choice = c(list_param$Trange[1], list_param$Trange[2])
    if(is.null(Trange_choice)) {
      if(missing(Trange)) stop("\"Trange\" shoulde be provided") else Trange_choice = Trange
    }
    filter_id <- newx[, T.name] >= Trange_choice[1]
    filter_id <- newx[filter_id, T.name] <= Trange_choice[2]
    newx <- newx[filter_id, ]

    interval_choice = list_param$interval
    if(is.null(interval_choice)) {
      if(missing(interval)) {
        stop("\"interval\" shoulde be provided")
      } else {
        interval_choice = match.arg(interval, choices = c("Original", "Standard"))
      }
    }

    switch(interval_choice,
           "Standard" = {
             ## mapping time sequence on to [0,1]
             newx[, T.name] = (newx[, T.name] - min(Trange_choice)) / diff(Trange_choice)
             interval_choice = c(0, 1)
             Trange_choice = c(0, 1)
           },
           "Original" = {
             interval_choice = Trange_choice
           })
    # <<< integral domain >>>

    # >>> integral interpolation: sseq <<<
    insert_choice = list_param$insert
    if(is.null(insert_choice)) {
      if(missing(insert)) {
        stop("\"insert\" shoulde be provided")
      } else {
        insert_choice = match.arg(insert, choices = c("FALSE", "X", "basis"))
      }
    }

    sseq_choice <- object$sseq
    if(is.null(sseq_choice)) {
      if(missing(sseq)) stop("\"sseq\" shoulde be provided") else sseq_choice = sort(unique(sseq))
    }
    sseq_choice <- sort(unique(c(round(Trange_choice, digits = digits), as.vector(newx[, T.name]), sseq_choice)))


    # if(insert != "FALSE") sseq_choice <- round(seq(from = interval_choice[1], to = interval_choice[2],
    #                                            by = min(diff(sseq)_choice)/10 # by = length(sseq_choice) * 2
    #                                            ),
    #                                            digits = digits)
    # <<< integral interpolation: sseq >>>

    # >>> integral basis <<<
    basis_fun_choice = list_param$basis_fun
    if(is.null(basis_fun_choice)) {
      if(missing(basis_fun)) {
        stop("\"basis_fun\" shoulde be provided")
      } else {
        basis_fun_choice = match.arg(basis_fun, choices = c("bs", "OBasis", "fourier"))
      }
    }

    degree_choice <- list_param$degree
    if(is.null(degree_choice)) {
      if(missing(degree)) stop("\"degree\" shoulde be provided") else degree_choice = degree
    }
    basis <- switch(basis_fun_choice,
                    "bs" = bs(x = sseq_choice, df = k_choice, degree = degree_choice,
                              Boundary.knots = interval_choice, intercept = TRUE),
                    "fourier" = eval.basis(sseq_choice,
                                           basisobj = create.fourier.basis(rangeval = interval_choice, nbasis = k_choice )),
                    "OBasis" = {
                      nknots <- k_choice - (degree_choice + 1)
                      if(nknots > 0) {
                        knots <- ((1:nknots) / (nknots + 1)) * diff(interval_choice) + interval_choice[1]
                      } else  knots <- NULL

                      evaluate(OBasis(expand.knots(c(interval_choice[1], knots, interval_choice[2])),
                                      order = degree_choice + 1),
                               sseq_choice)
                    })
    # <<< integral basis >>>

    # >>> integral <<< #
    D <- split(newx[, c(T.name, X.names)], newx[, ID.name])
    n <- length(D)
    method_choice = list_param$method
    if(is.null(method_choice)) {
      if(missing(method)) {
        stop("\"method\" shoulde be provided")
      } else {
        method_choice = match.arg(method, choices = c("trapezoidal", "step"))
      }
    }
    newZ <- matrix(NA, nrow = n, ncol = p1 - pc - 1)
    for(i in 1:n) newZ[i, ] <- ITG(D[[i]], basis, sseq_choice, T.name, interval_choice, insert_choice, method_choice)$integral
    # <<< integral >>> #
  }

  # >>> prediction <<< #
  newX = cbind(newZ, newZc)
  if(!is.null(s)) {
    beta <- lamtobeta(lambda = object$lam, beta = beta, s = s)
  }
  fitting <- as.matrix(cbind2(newX, 1)) %*% beta
  # <<< prediction >>> #

  return(fitting)
}




#' @title
#' Make predictions based on a \code{"cv.FuncompCGL"} object.
#'
#' @description
#' This function makes prediction based on a \code{cv.FuncompCGL}
#' object, using the stored \code{"FuncompCGL.fit"} object and the optimal values of the
#' regularization parameter \code{lam} and the degrees of freedom \code{k}.
#'
#' @usage
#' \method{predict}{cv.FuncompCGL}(object, Znew, Zcnew = NULL,
#'         s = c("lam.1se", "lam.min"), k = NULL, trim = FALSE, ...)
#'
#' @param Znew data frame or matrix \code{X} as in \code{FuncompCGL} with new
#'             functional compositional data at which prediction is to be made.
#'
#' @param Zcnew matrix \code{Zc} as in \code{FuncompCGL} with new values of
#'              time-invariate covariates at which prediction is to be made.
#'              Default is \code{NULL}.
#'
#' @param ... Other arguments passed to \code{predict.FuncompCGL}
#'
#' @inheritParams coef.cv.FuncompCGL
#'
#' @seealso
#' \code{\link{cv.FuncompCGL}} and \code{\link{FuncompCGL}}, and
#' \code{\link[=coef.cv.FuncompCGL]{coef}} and
#' \code{\link[=plot.cv.FuncompCGL]{plot}} methods for \code{"cv.FuncompCGL"} object.
#'
#' @return
#' The prediction values at the requested value(s) for \code{s} and \code{k}.
#' If \code{k} is a vector, a list of prediction matrix is returned,
#' otherwise a prediction matrix is returned.
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
#' Data <- Fcomp_Model(n = n_train, p = p, m = 0, intercept = TRUE,
#'                     SNR = 4, sigma = 3, rho_X = 0, rho_T = 0.6, df_beta = df_beta,
#'                     n_T = 20, obs_spar = 1, theta.add = FALSE,
#'                     beta_C = as.vector(t(beta_C_true)))
#' arg_list <- as.list(Data$call)[-1]
#' arg_list$n <- n_test
#' Test <- do.call(Fcomp_Model, arg_list)
#' k_list = c(4,5)
#' cv_m1 <-  cv.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
#'                         Zc = Data$data$Zc, intercept = Data$data$intercept,
#'                         k = k_list, nfolds = 5)
#' y_hat = predict(cv_m1, Znew = Test$data$Comp, Zcnew = Test$data$Zc)
#' predict(cv_m1, Znew = Test$data$Comp, Zcnew = Test$data$Zc, s = "lam.1se")
#' predict(cv_m1, Znew = Test$data$Comp, Zcnew = Test$data$Zc, s = c(0.5, 0.1, 0.05), k = k_list)
#' plot(Test$data$y, y_hat, xlab = "Observed Response", ylab = "Predicted Response")
#' }
#'
#' @inherit predict.FuncompCGL details references author
#'
#' @export

predict.cv.FuncompCGL <- function(object, Znew, Zcnew = NULL,
                                  s = c("lam.1se", "lam.min"), k = NULL,
                                  trim = FALSE, ...) {

  trim <- ifelse(trim, "Ttrim", "Ftrim")

  if (is.numeric(s) || is.null(s) ) {
    lam_se <- s
    if(is.null(k)) {
      stop("Need to provide df k")
    } else {
      if(is.numeric(k)){
        if(FALSE %in% (k %in% as.integer(names(object$FuncompCGL.fit)))) stop("element in \"K\" is out of range")
        k_se <- k
      }
      k_se <- as.character(k_se)
    }
  } else if (is.character(s)) {
    s <- match.arg(s)
    lam_se <- object[[trim]][[s]]['lam']
    k_se <- object[[trim]][[s]]['df']
  } else stop("Invalid form for s")

  if(length(k_se) > 1) {
    predmat <- as.list(k_se)
    for(i in seq(length(k_se))) {
      predmat[[i]] <- predict(object$FuncompCGL.fit[[ which(names(object$FuncompCGL.fit) == k_se[i])]],
                              Znew, Zcnew, s = lam_se, ...)
    }
    names(predmat) <- paste("k", k_se, sep = "=")
  } else {
    predmat <- predict(object$FuncompCGL.fit[[which(names(object$FuncompCGL.fit) == k_se)]],
                       Znew, Zcnew, s = lam_se, ...)
  }

  return(predmat)
}



#' @title
#' Make predictions based on a \code{"GIC.FuncompCGL"} object.
#'
#' @description
#' This function makes prediction based on a \code{"GIC.FuncompCGL"} object, using the
#' stored \code{"FuncompCGL.fit"} object and the optimal values of
#' the regularization parameter \code{lam} and the degrees of freedom \code{k}.
#'
#' @usage
#' \method{predict}{GIC.FuncompCGL}(object, Znew, Zcnew = NULL,
#'         s = "lam.min", k = NULL, ...)
#'
#' @param Znew data frame or matrix \code{X} as in \code{FuncompCGL} with new
#'             functional compositional data at which prediction is to be made.
#'
#' @param Zcnew matrix \code{Zc} as in \code{FuncompCGL} with new values of
#'              time-invariate covariates at which prediction is to be made.
#'              Default is \code{NULL}.
#'
#' @param ... Other arguments passed to \code{predict.FuncompCGL}
#'
#' @inheritParams coef.GIC.FuncompCGL
#'
#' @seealso
#' \code{\link{GIC.FuncompCGL}} and \code{\link{FuncompCGL}}, and
#' \code{\link[=coef.GIC.FuncompCGL]{coef}} and
#' \code{\link[=plot.GIC.FuncompCGL]{plot}} methods for \code{"GIC.FuncompCGL"} object.
#'
#' @examples
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
#'                     SNR = 4, sigma = 3, rho_X = 0.6, rho_T = 0,
#'                     df_beta = df_beta, n_T = 20, obs_spar = 1, theta.add = FALSE,
#'                     beta_C = as.vector(t(beta_C_true)))
#' arg_list <- as.list(Data$call)[-1]
#' arg_list$n <- n_test
#' Test <- do.call(Fcomp_Model, arg_list)
#' GIC_m1 <-  GIC.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
#'                           Zc = Data$data$Zc, intercept = Data$data$intercept,
#'                           k = k_list)
#' y_hat <- predict(GIC_m1, Znew = Test$data$Comp, Zcnew = Test$data$Zc)
#' predict(GIC_m1, Znew = Test$data$Comp, Zcnew = Test$data$Zc, s = NULL, k = k_list)
#' plot(Test$data$y, y_hat, xlab = "Observed response", ylab = "Predicted response")
#'
#' @inherit predict.cv.FuncompCGL details return references author
#'
#' @export

predict.GIC.FuncompCGL <- function(object, Znew, Zcnew = NULL,
                                   s = "lam.min", k = NULL, ...) {

  if (is.numeric(s) || is.null(s) ){
    lam_se <- s
    if(is.null(k)) {
      stop("Need to provide df k")
    } else {
      if(is.numeric(k)){
        if(FALSE %in% (k %in% as.integer(names(object$FuncompCGL.fit)))) stop("element in \"K\" is out of range")
        k_se <- k
      }
      k_se <- as.character(k_se)
    }
  } else if (is.character(s)) {
    s <- match.arg(s)
    lam_se <- object$lam.min['lam']
    k_se <- object$lam.min['df']
  } else stop("Invalid form for s")

  if(length(k_se) > 1) {
    predmat <- as.list(k_se)
    for(i in seq(length(k_se))) {
      predmat[[i]] <- predict(object$FuncompCGL.fit[[which(names(object$FuncompCGL.fit) == k_se[i])]],
                              Znew, Zcnew, s = lam_se, ...)
    }
    names(predmat) <- paste("k", k_se, sep = "=")
  } else {
    predmat <- predict(object$FuncompCGL.fit[[ which(names(object$FuncompCGL.fit) == k_se)]],
                       Znew, Zcnew, s = lam_se, ...)
  }
  return(predmat)
}


#' @title
#' Make predictions based on a \code{"compCL"} object.
#'
#' @description
#' Make predictions based on a fitted \code{"\link{compCL}"} object.
#'
#' @param object fitted \code{"\link{compCL}"} object.
#'
#' @param Znew \code{z} matrix as in \code{compCL} with new compositional data
#'          or categorical data.
#' @param Zcnew \code{Zc} matrix as in \code{compCL} with new data for other
#'              covariates. Default is \code{NULL}
#'
#' @param s value(s) of the penalty parameter \code{lam} at which predictions are required.
#'          Default is the entire sequence used in the fitted object.
#'
#' @param \dots not used.
#'
#' @details
#' \code{s} is the vector at which predictions are requested. If \code{s} is not in the lambda
#' sequence used for fitting the model, the \code{predict} function uses linear interpolation.
#'
#' @return
#' predicted values at the requested values of \code{s}.
#'
#' @seealso
#' \code{\link{compCL}} and \code{\link[=coef.compCL]{coef}},
#' \code{\link[=print.compCL]{predict}} and \code{\link[=plot.compCL]{plot}} methods
#' for \code{"compCL"} object.
#'
#' @examples
#' Comp_data = comp_Model(n = 50, p = 30)
#' Comp_data2 = comp_Model(n = 30, p = 30, beta = Comp_data$beta)
#' m1 = compCL(y = Comp_data$y, Z = Comp_data$X.comp,
#'             Zc = Comp_data$Zc, intercept = Comp_data$intercept)
#' predict(m1, Znew = Comp_data2$X.comp, Zcnew = Comp_data2$Zc)
#' predict(m1, Znew = Comp_data2$X.comp, Zcnew = Comp_data2$Zc, s = c(1, 0.5, 0.1))
#'
#' @inherit coef.compCL references author
#'
#' @export

predict.compCL <- function(object, Znew, Zcnew = NULL,  s = NULL, ...) {
  Znew <- proc.comp(Znew) # log is taken
  if( !is.null(Zcnew) && is.null(dim(Zcnew)) ) {
    #if only one new data point
    newZc = matrix(Zcnew, nrow = 1)
  }
  predict.linear(object, cbind(Znew, Zcnew), s)
}



#' @title
#' Make predictions based on a \code{"cv.compCL"} object.
#'
#' @description
#' This function makes prediction based on a cross-validated \code{compCL} model,
#' using the stored \code{compCL.fit} object.
#'
#' @usage
#' \method{predict}{cv.compCL}(object, Znew, Zcnew = NULL, s = c("lam.min", "lam.1se" ),
#'         trim = FALSE, ...)
#'
#' @param object fitted \code{"cv.compCL"} model.
#'
#' @param s specify the \code{lam} at which prediction(s) is requested.
#'          \itemize{
#'            \item \code{s = "lam.min"} (default), value of \code{lam} that obtains
#'                  the minimun value of cross-validation error.
#'            \item \code{s = "lam.1se"} value of \code{lam} that obtains 1 standard
#'                  error above the miminum of the cross-validation errors.
#'            \item if \code{s} is numeric, it is taken as the value(s) of lam to be used.
#'            \item if \code{s = NULL}, uses the whole sequence of \code{lam} stored in the
#'                  \code{"cv.compCL"} object.
#'          }
#'
#' @param trim Whether to use the trimmed result. Default is FASLE.
#'
#' @inheritParams predict.compCL
#'
#' @seealso
#' \code{\link{cv.compCL}} and \code{\link{compCL}},
#' and \code{\link[=coef.cv.compCL]{coef}} and \code{\link[=plot.cv.compCL]{plot}} methods
#' for \code{"cv.compCL"} object.
#'
#' @examples
#' p = 30
#' n = 50
#' beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
#' beta = c( beta, rep(0, times = p - length(beta)) )
#' Comp_data = comp_Model(n = n, p = p, beta = beta, intercept = FALSE)
#' test_data = comp_Model(n = 30, p = p, beta = beta, intercept = FALSE)
#' cvm1 <- cv.compCL(y = Comp_data$y, Z = Comp_data$X.comp,
#'                   Zc = Comp_data$Zc, intercept = Comp_data$intercept)
#' y_hat = predict(cvm1, Znew = test_data$X.comp, Zcnew = test_data$Zc)
#' predmat = predict(cvm1, Znew = test_data$X.comp, Zcnew = test_data$Zc, s = NULL)
#' plot(test_data$y, y_hat, xlab = "Observed response", ylab = "Predicted response")
#' abline(a = 0, b = 1, col = "red")
#'
#' @inherit predict.compCL details return references author
#'
#' @export

predict.cv.compCL <- function(object, Znew, Zcnew = NULL, s = c("lam.min", "lam.1se" ),
                              trim = FALSE, ...){
  if (is.numeric(s) || is.null(s) ) {
    lam <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    trim <- ifelse(trim, "Ttrim", "Ftrim")
    lam <- object[[trim]][[s]]
  } else {
    stop("Invalid form for s")
  }
  predict(object$compCL.fit, Znew, Zcnew, s = lam, ...)
}


#' @title
#' Make predictions based on a \code{"GIC.compCL"} object.
#'
#' @description
#' This function makes prediction based on a \code{"GIC.compCL"} model,
#' using the stored \code{"compCL.fit"} object and the optimal value of \code{lambda}.
#'
#' @param object fitted \code{"GIC.compCL"} model.
#'
#' @param s specify the \code{lam} at which prediction(s) is requested.
#'          \itemize{
#'            \item \code{s = "lam.min"} (default), \code{lam} that obtains the minimun value of GIC values.
#'            \item if \code{s} is numeric, it is taken as the value(s) of lam to be used.
#'            \item if \code{s = NULL}, uses the whole sequence of \code{lam} stored in the
#'                  \code{"GIC.compCL"} object.
#'          }
#'
#' @inheritParams predict.cv.compCL
#'
#' @seealso
#' \code{\link{GIC.compCL}} and \code{\link{compCL}}, and
#' \code{\link[=coef.GIC.compCL]{coef}} and
#' \code{\link[=plot.GIC.compCL]{plot}} methods for \code{"GIC.compCL"}.
#'
#' @examples
#' p = 30
#' n = 50
#' beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
#' beta = c(beta, rep(0, times = p - length(beta)))
#' Comp_data = comp_Model(n = n, p = p, beta = beta, intercept = FALSE)
#' test_data = comp_Model(n = 100, p = p, beta = beta, intercept = FALSE)
#' GICm1 <- GIC.compCL(y = Comp_data$y, Z = Comp_data$X.comp,
#'                     Zc = Comp_data$Zc, intercept = Comp_data$intercept)
#' y_hat = predict(GICm1, Znew = test_data$X.comp, Zcnew = test_data$Zc)
#' predmat = predict(GICm1, Znew = test_data$X.comp, Zcnew = test_data$Zc, s = c(1, 0.5, 1))
#' plot(test_data$y, y_hat, xlab = "Observed value", ylab = "Predicted value")
#' abline(a = 0, b = 1, col = "red")
#'
#' @inherit predict.compCL details return references author
#'
#' @export


predict.GIC.compCL <- function(object, Znew, Zcnew = NULL, s = "lam.min", ...){
  if (is.numeric(s) || is.null(s) ) {
    lam <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    lam <- object$lam.min
  } else {
    stop("Invalid form for s")
  }

  predict(object$compCL.fit, Znew, Zcnew, s = lam, ...)
}
#### <<< predict method >>> ####




#### >>> plot method <<< ####


#' @title
#' Plot solution paths from a \code{"FuncompCGL"} object.
#'
#' @description
#' Produce a coefficient profile plot of the coefficient paths for a fitted
#' \code{"FuncompCGL"} object.
#'
#' @param x fitted \code{"FuncompCGL"} object.
#' @param ylab what is the on Y-axis,
#'             \code{"L2"} (default) plots against the L2-norm of each group of coefficients,
#'             \code{"L1"} against L1-norm.
#' @param xlab what is on the X-axis,
#'             \code{"log"} plots against \code{log(lambda)} (default),
#'             \code{"-log"} against \code{-log(lambda)}, and \code{"lambda"} against \code{lambda}.
#'
#' @param \dots other graphical parameters.
#'
#' @details
#' A solution path plot is produced.
#'
#' @return
#' No return value. Side effect is a base R plot.
#'
#' @seealso
#' \code{\link{FuncompCGL}}, and \code{\link[=predict.FuncompCGL]{predict}},
#' \code{\link[=coef.FuncompCGL]{coef}} and \code{\link[=print.FuncompCGL]{print}}
#' methods for \code{"FuncompCGL"} object.
#'
#' @examples
#' df_beta = 5
#' p = 30
#' beta_C_true = matrix(0, nrow = p, ncol = df_beta)
#' beta_C_true[1, ] <- c(-0.5, -0.5, -0.5 , -1, -1)
#' beta_C_true[2, ] <- c(0.8, 0.8,  0.7,  0.6,  0.6)
#' beta_C_true[3, ] <- c(-0.8, -0.8 , 0.4 , 1 , 1)
#' beta_C_true[4, ] <- c(0.5, 0.5, -0.6  ,-0.6, -0.6)
#' Data <- Fcomp_Model(n = 50, p = p, m = 0, intercept = TRUE,
#'                     SNR = 4, sigma = 3, rho_X = 0, rho_T = 0.6, df_beta = df_beta,
#'                     n_T = 20, obs_spar = 1, theta.add = FALSE,
#'                     beta_C = as.vector(t(beta_C_true)))
#' m1 <- FuncompCGL(y = Data$data$y, X = Data$data$Comp, Zc = Data$data$Zc,
#'                  intercept = Data$data$intercept, k = df_beta, tol = 1e-10)
#' plot(m1)
#' plot(m1, ylab = "L1", xlab = "-log")
#'
#' @inherit coef.FuncompCGL author references
#'
#' @export

## TODO: add "coef" and "L2_curve" features on ylab

plot.FuncompCGL <- function(x,
                            ylab = c("L2", "L1"),
                            xlab = c("log", "-log", "lambda"),
                            ...) {
  obj <- x
  ylab = match.arg(ylab)
  xlab = match.arg(xlab)

  obj$call <- as.list(obj$call)
  k <- obj$call[['k']]
  p1 <- ncol(obj$Z)
  p <- p1 / k

  ## implicitly converting baseline result into full-long
  beta <- coef(obj)
  ## groups that ever been non-zero
  include_which <- sort(unique(unlist(apply(beta, 2, Nzero, p = p, k = k))))
  group <- matrix(1:p1, nrow = k)
  beta <- beta[group[, include_which], ]

  switch(ylab,
         "L1" = {
           yvalue = apply(beta, 2, function(x, k) {
             value = vector()
             for(l in 1:(length(x)/k)) value[l] <- sum(abs(x[(l*k-k+1):(l*k)]))
             return(value)
           }, k = k)
           ylab = "L1 norm"},
         "L2" = {
           yvalue = apply(beta, 2, function(x, k) {
             value = vector()
             for(l in 1:(length(x)/k)) value[l] <- sqrt(sum((x[(l*k-k+1):(l*k)])^2))
             return(value)
           }, k = k)
           ylab = "L2 norm"
         }
  )

  switch(xlab,
         "lambda" = {
           xlab = "Lambda"
           xvalue = drop(obj$lam)
         },
         "log" = {
           xlab = "Log(Lambda)"
           xvalue = log(drop(obj$lam))
           },
          "-log" = {
            xlab = "Log(Lambda)"
            xvalue = -log(drop(obj$lam))
          }
         )


  dotlist=list(...) #list(...)
  type=dotlist$type
  if(is.null(type))
    matplot(xvalue,t(yvalue),lty=1,xlab=xlab,ylab=ylab,type="l", ylim = c(-0.1, max(yvalue)), ...
    )  else matplot(xvalue,t(yvalue),lty=1,xlab=xlab,ylab=ylab
                    , ylim = c(-0.1, max(yvalue)),...
    )

  #text(x = min(xvalue), y = yvalue[,dim(yvalue)[2]], labels = paste(include_which))
  # cvraw <- (drop(y) - cv.m1$fit.preval[, , 1])^2
  # N <- length(y) - apply(is.na(cv.m1$fit.preval[, , 1]), 2, sum)
  # cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
  # cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 1))

  at.txt = apply(yvalue, 1, function(x) which(x>0)[1])
  at.index <- rbind(at.txt, sort(at.txt))
  #matplot(xvalue,t(yvalue),lty=1,xlab=xlab,ylab=ylab,type="l", ylim = c(-0.4, max(yvalue)))
  cex = 1
  text(x = xvalue[at.txt[order(at.txt)[-((1:floor(length(at.txt)/2))*2)]]],
       y = -0.05, labels = paste(include_which[order(at.txt)][-((1:floor(length(at.txt)/2))*2)]),
       col = rep(c("black", "grey"), length.out = length(include_which) - length((1:floor(length(at.txt)/2))*2)),
       cex = cex)
  text(x = xvalue[at.txt[order(at.txt)[((1:floor(length(at.txt)/2))*2)]]],
       y = -0.1, labels = paste(include_which[order(at.txt)][(1:floor(length(at.txt)/2))*2]),
       col = rep(c("black", "grey"), length.out = length((1:floor(length(at.txt)/2))*2)),
       cex = cex)

}


#' @title
#' Plot the cross-validation curve produced by \code{"cv.FuncompCGL"}.
#'
#' @description
#' Plot the cross-validation curve with its upper and lower standard
#' deviation curves.
#'
#' @param x fitted \code{"cv.FuncompCGL"} model.
#'
#' @param trim logical; whether to use the trimmed result.
#'        Default is \code{FALSE}.
#'
#' @param k a vector or character string
#'               \itemize{
#'               \item if character string, either \code{"lam.1se"} or
#'                     \code{"lam.min"}.
#'               \item if it is an integer vector, specify the set of degrees of freedom
#'                     \code{k} to plot.
#'               \item if it is missing (default), cross-validation curves for \code{k}
#'                     that are associated with \code{lambda.min} (blue)
#'                     and \code{lambda.1se} (red) are plotted.
#'               }
#'
#' @inheritParams plot.FuncompCGL
#'
#' @details
#' A cross-validation curve is produced.
#'
#' @inherit plot.FuncompCGL return
#'
#' @seealso
#' \code{\link{cv.FuncompCGL}} and \code{\link{FuncompCGL}}, and
#' \code{\link[=predict.cv.FuncompCGL]{predict}} and
#' \code{\link[=coef.cv.FuncompCGL]{coef}} methods for \code{"cv.FuncompCGL"} object.
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
#' Data <- Fcomp_Model(n = 50, p = p, m = 0, intercept = TRUE,
#'                     SNR = 4, sigma = 3, rho_X = 0, rho_T = 0.6, df_beta = df_beta,
#'                     n_T = 20, obs_spar = 1, theta.add = FALSE,
#'                     beta_C = as.vector(t(beta_C_true)))
#' k_list <- 4:5
#' cv_m1 <-  cv.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
#'                         Zc = Data$data$Zc, intercept = Data$data$intercept,
#'                         k = k_list, nfolds = 5, keep = TRUE)
#' plot(cv_m1)
#' plot(cv_m1, xlab = "-log", k = k_list)
#' }
#'
#' @inherit coef.FuncompCGL author references
#'
#' @export
#'
#' @import graphics


plot.cv.FuncompCGL <- function(x, xlab = c("log", "-log", "lambda"),
                               trim = FALSE, k, ...) {
  cvobj = x
  xlab <- match.arg(xlab)
  trim <- ifelse(trim, "Ttrim", "Ftrim")
  cvobj_use <- cvobj[[trim]]

  # >>> x-axis matching <<< #
  switch(xlab,
         "log" = {
           xlab = "Log(Lambda)"
           xvalue = log(drop(cvobj$lam))
         },
         "-log" = {
           xlab = "-Log(Lambda)"
           xvalue = -log(drop(cvobj$lam))
         },
         "lambda" = {
           xlab = "Lambda"
           xvalue = drop(cvobj$lam)
         })
  # <<< x-axis matching >>> #

  # >>> k value(s) <<< #
  k_all <- as.numeric(names(cvobj$FuncompCGL.fit))
  rule_min <- cvobj_use$lam.min
  rule_1se <- cvobj_use$lam.1se
  if (missing(k)) {
    k <- c(rule_min['df'], rule_1se['df'])
  } else if (is.numeric(k)) {
    k <- k[k %in% k_all]
    if(length(k) == 0) stop("Elements in \"k\" MUST be one(s) of these used in model fitting.")
  } else if (is.character(k)) {
    k <- match.arg(k, choices = c("lam.min", "lam.1se"))
    k <- switch(k,
                "lam.1se" = rule_1se['df'],
                "lam.min" = rule_min['df'])
  }
  emp1 <- which(rule_min['df'] == k)
  emp2 <- which(rule_1se['df'] == k)
  normal <- which(!(k %in% c(rule_min['df'], rule_1se['df'])))
  k <- c(sort(k[normal]), k[emp2], k[emp1])
  k_id <- match(k, k_all)
  # <<< k value(s) >>> #

  # >>> pre-process <<< #
  N <- apply(!is.na(cvobj_use$cvm[k_id, , drop = FALSE]), 1, function(x) max(which(x)))
  N_id <- seq(max(N))
  xvalue = xvalue[N_id]
  cvobj_use <- lapply(cvobj_use[c("cvm", "cvup", "cvlo")], function(x, slice) x[, slice], slice = N_id)
  ## red for lam.min; blue for lam.1se
  bar_col <- c(rep("lightgrey", length(normal)), rep("blue", length(emp2)), rep("red", length(emp1)))
  bar_width <- 0.01
  cuv_col <- c(rep("limegreen", length(normal)), rep("blue", length(emp2)), rep("red", length(emp1)))
  cuv_pch <- c(rep(18, length(normal)), rep(20, length(emp2)), rep(20, length(emp1)))
  # <<< pre-process >>> #

  # >>> plot <<< #
  plot.args = list(xlab = xlab, ylab = "MEAN-Squared Error", type = "n", 0,
                   xlim = range(xvalue),
                   ylim = range(cvobj_use$cvup[k_id, ], cvobj_use$cvlo[k_id, ], na.rm = TRUE))
  new.args = list(...)
  if(length(new.args)) plot.args[names(new.args)] = new.args
  do.call("plot", plot.args)

  ## base layer: error bars
  for(l in 1:length(k)) {
    error.bars(xvalue, cvobj_use$cvup[k_id[l], ], cvobj_use$cvlo[k_id[l], ],
              width = bar_width, col = bar_col[l])
  }

  ## top layer: error curve
  for(l in 1:length(k)) {
    points(xvalue, cvobj_use$cvm[k_id[l], ], pch = cuv_pch[l], col = cuv_col[l])
  }

  # <<< plot >>> #

  # >>> label DF along path: axis = 3 <<< #
  if(length(emp2)) {
    labels <- drop(cvobj$FuncompCGL.fit[[as.character(rule_1se['df'])]]$df)[N_id]
    axis(side = 3, tick = FALSE, col.axis = "blue", line = -0.4, # pos = plot.args$ylim[2],
         at = xvalue, labels = as.character(labels))
    abline(v = switch(xlab,
                      "Lambda" = rule_1se["lam"],
                      "Log(Lambda)" = log(rule_1se["lam"]),
                      "-Log(Lambda)" = -log(rule_1se["lam"])),
           lty = 3, col = "blue")
  }

  if(length(emp1)) {
    labels <- drop(cvobj$FuncompCGL.fit[[as.character(rule_min['df'])]]$df)[N_id]
    axis(side = 3, tick = FALSE, line = 0.5, col.axis = "red",
         at = xvalue, labels = as.character(labels))
    abline(v = switch(xlab,
                      "Lambda" = rule_min["lam"],
                      "Log(Lambda)" = log(rule_min["lam"]),
                      "-Log(Lambda)" = -log(rule_min["lam"])),
           lty = 3, col = "red")
  }
  # <<< label DF along path: axis = 3 >>> #

}


#' @title
#' Plot the GIC curve produced by \code{"GIC.FuncompCGL"} object.
#'
#' @description
#' Plot the GIC curve as a function of the \code{lam} values used for different
#' degree of freedom \code{k}.
#'
#' @param x fitted \code{"\link{GIC.FuncompCGL}"} object.
#'
#' @param k value(s) of the degrees of freedom at which GIC cuvre(s) are plotted.
#'          \itemize{
#'          \item if missing (default), GIC curve for \code{k} that is associated
#'                with \code{"lam.min"} (RED) stored on \code{x} is plotted.
#'          \item if it is an integer vector, specify what set of degrees of freedom
#'                to plot.
#'          }
#'
#' @inheritParams plot.cv.FuncompCGL
#'
#' @details
#' A GIC curve is produced.
#'
#' @inherit plot.FuncompCGL return
#'
#' @examples
#' df_beta = 5
#' p = 30
#' beta_C_true = matrix(0, nrow = p, ncol = df_beta)
#' beta_C_true[1, ] <- c(-0.5, -0.5, -0.5 , -1, -1)
#' beta_C_true[2, ] <- c(0.8, 0.8,  0.7,  0.6,  0.6)
#' beta_C_true[3, ] <- c(-0.8, -0.8 , 0.4 , 1 , 1)
#' beta_C_true[4, ] <- c(0.5, 0.5, -0.6  ,-0.6, -0.6)
#' Data <- Fcomp_Model(n = 50, p = p, m = 0, intercept = TRUE,
#'                     SNR = 4, sigma = 3, rho_X = 0.6, rho_T = 0,
#'                     df_beta = df_beta, n_T = 20, obs_spar = 1, theta.add = FALSE,
#'                     beta_C = as.vector(t(beta_C_true)))
#'
#' k_list <- c(4,5)
#' GIC_m1 <-  GIC.FuncompCGL(y = Data$data$y, X = Data$data$Comp,
#'                           Zc = Data$data$Zc, intercept = Data$data$intercept,
#'                           k = k_list)
#' plot(GIC_m1)
#' plot(GIC_m1, xlab = "-log", k = k_list)
#'
#' @seealso
#' \code{\link{GIC.FuncompCGL}} and \code{\link{FuncompCGL}}, and
#' \code{\link[=predict.GIC.FuncompCGL]{predict}} and
#' \code{\link[=coef.GIC.FuncompCGL]{coef}} methods for \code{"GIC.FuncompCGL"} object.
#'
#' @inherit coef.FuncompCGL author references
#'
#' @export
#'
#' @importFrom grDevices rainbow
#' @importFrom utils tail


plot.GIC.FuncompCGL <- function(x, xlab = c("log", "-log", "lambda"), k, ...) {
  ## >>> x-axis <<<
  xlab <- match.arg(xlab)
  switch(xlab,
         "lambda" = {
           xlab = "Lambda"
           xvalue = drop(x$lam)
           legend_pos = "topright"
           line_value = x$lam.min['lam']
         },
         "log" = {
           xlab = "Log(Lambda)"
           xvalue = log(drop(x$lam))
           legend_pos = "topright"
           line_value = log(x$lam.min['lam'])
         },
         "-log" = {
           xlab = "-Log(Lambda)"
           xvalue = -log(drop(x$lam))
           legend_pos = "topleft"
           line_value = -log(x$lam.min['lam'])
         })
  ## <<< x-axis >>>

  ## >>> K value(s) <<<
  k_all <- as.numeric(names(x$FuncompCGL.fit))
  k_opt <- x$lam.min['df']
  if(missing(k)) k <- k_opt
  if(is.numeric(k)) {
    k <- k[k %in% k_all]
    if(length(k) == 0) stop("Elements in \"k\" MUST be one(s) of these used in model fitting.")
  }
  emp <- which(k_opt == k)
  normal <- which(k_opt != k)
  k <- c(sort(k[normal]), k[emp]) # k[c(normal, emp)]
  k_id <-  match(k, k_all)
  ## <<< K value(s) >>>

  # >>> pre-process <<< #
  N <- apply(!is.na(x$GIC[k_id, , drop = FALSE]), 1, function(x) max(which(x)))
  N_id <- seq(max(N))
  xvalue <- xvalue[N_id]
  Yvalue <- x$GIC[, N_id]
  ## reserve pch = 1 for k_opt
  pch = c(0, seq(length(k))[-1])
  ## reserve vcol = "red" ("#FF0000FF") for k_opt
  col = rainbow(length(k)+1)[-1]
  ## reverse order to let k_opt, if included, list as the first
  col <- c(switch(length(emp) == 0 + 1, "#FF0000FF", NULL), head(col, length(normal)))
  pch <- c(switch(length(emp) == 0 + 1, 1, NULL), head(pch, length(normal)))
  text <- paste("k", rev(k), sep = "=")
  ylab = unlist(strsplit(rownames(x$GIC)[1], split = "_"))[1]
  ylab = paste(ylab, "curve(s)")
  # <<< pre-process >>> #

  ## >>> plot <<<
  plot.args = list(xlab = xlab, ylab = ylab, type = "n", 0,
                   xlim = range(xvalue), ylim = range(Yvalue[k_id, ], na.rm = TRUE))
  new.args = list(...)
  if(length(new.args)) plot.args[names(new.args)] = new.args
  do.call("plot", plot.args)
  for(l in 1:length(k)) points(xvalue, Yvalue[k_id[l], ], pch = pch[l], col = col[l])
  if(length(emp)) {
    labels = drop(x$FuncompCGL.fit[[as.character(k_opt)]]$df)[N_id]
    axis(side = 3, tick = FALSE, col.axis = "red", line = 0.2,
         at = xvalue, labels = as.character(labels))
    abline(v = line_value, lty = 3, col = "red")
  }

  legend(legend_pos, inset = 0.02, cex = 0.8,
         box.lty = 2, box.lwd = 0.1, box.col = "grey",
         legend = text, col = col, pch = pch)
  # <<<  plot >>> #
}



#' @title
#' Plot solution paths from a \code{"compCL"} object.
#'
#' @description
#' Produce a coefficient profile plot from a fitted
#' \code{"\link{compCL}"} object.
#'
#' @param x fitted \code{"\link{compCL}"} model.
#' @param xlab what is on the X-axis. \code{"lam"} plots against the log-lambda sequence
#'             (default) and \code{"norm"} against the L1-norm of the coefficients.
#' @param label if TRUE, label the curve with the variable sequence numbers. Default is \code{FALSE}.
#' @inheritParams plot.FuncompCGL
#'
#'
#' @details
#' A coefficient profile plot for the compositional predictors is produced.
#'
#' @return
#' No return value. Side effect is a base R plot.
#'
#' @seealso
#' \code{\link{compCL}} and \code{\link[=print.compCL]{print}},
#' \code{\link[=predict.compCL]{predict}} and \code{\link[=coef.compCL]{coef}} methods
#' for \code{"compCL"} object.
#'
#' @examples
#' Comp_data = comp_Model(n = 50, p = 30)
#' Comp_fit = compCL(y = Comp_data$y, Z = Comp_data$X.comp, Zc = Comp_data$Zc,
#'                   intercept = Comp_data$intercept)
#' plot(Comp_fit)
#' plot(Comp_fit, xlab = "norm", label = TRUE)
#'
#' @inherit coef.compCL author references
#'
#' @export


plot.compCL <- function(x, xlab=c("lam", "norm"), label = FALSE, ...) {

  xlab <-  match.arg(xlab)
  p <- ncol(x$Z_log)
  beta <- x$beta[1:p, ]
  lambda <- x$lam
  df <- x$df

  # nr = nrow(beta)
  # if(nr == 1) {
  #   index_which <- ifelse(any(abs(beta) > 0), 1, NULL)
  # } else {
  #   index_which <- apply(beta, 1, function(x) ifelse(any(abs(x) > 0), TRUE, FALSE))
  #   index_which <- seq(nrow(beta))[index_which]
  # }

  index_which <- apply(beta, 1, function(x) ifelse(any(abs(x) > 0), TRUE, FALSE))
  index_which <- seq(nrow(beta))[index_which]
  nwhich = length(index_which)

  if(nwhich == 0) {
    warning("No plot produced since all coefficients are zero")
    return()
  }

  beta = as.matrix(beta[index_which, , drop=FALSE])
  name_ylab = "Coefficients"
  switch(xlab,
         "norm" = {
           index_xlab = apply(abs(beta), 2, sum)
           name_xlab = "L1 Norm"
           approx.f = 1
         },
         "lam" = {
           index_xlab = log(lambda)
           name_xlab = "log Lambda"
           approx.f = 0
         })

  dotlist = list(...)
  type = dotlist$type
  if(is.null(type)) {
    matplot(index_xlab, t(beta), lty=1, xlab = name_xlab, ylab = name_ylab, type = "l", ...)
  } else {
    matplot(index_xlab, t(beta), lty=1, xlab = name_xlab, ylab = name_ylab, ...)
  }

  atdf=pretty(index_xlab)
  ## compute df by interpolating to df at next smaller lambda thanks to R
  ## package glment and Yunyang Qian
  prettydf=approx(x=index_xlab,y=df,xout=atdf,rule=2,method="constant",f=approx.f)$y

  axis(3, at = atdf, labels = prettydf, tcl = NA, col.axis = "red", tick = FALSE)
  # mtext("DF", side=3, line=3, cex.lab=2,las=0, col="red")
  if(label){
    switch(xlab,
           "norm" = {
             xpos = max(index_xlab)
             pos = 4
           },
           "lam" = {
             xpos = min(index_xlab)
             pos = 2
           })
    pos_xlab = rep(xpos, nwhich)
    pos_ylab = beta[, ncol(beta)]
    text(pos_xlab, pos_ylab, paste(index_which), cex = 0.8, pos = pos)
  }
}



#' @title
#' Plot the cross-validation curve produced by \code{"cv.compCL"} object.
#'
#' @description
#' Plot the cross-validation curve with its upper and lower standard deviation
#' curves.
#'
#' @param x fitted \code{"cv.compCL"} object.
#'
#' @inheritParams plot.cv.FuncompCGL
#'
#' @details
#' A cross-validation curve is produced.
#'
#' @inherit plot.compCL return
#'
#' @seealso
#' \code{\link{cv.compCL}} and \code{\link{compCL}},
#' and \code{\link[=coef.cv.compCL]{coef}} and \code{\link[=plot.cv.compCL]{plot}} methods
#' for \code{"cv.compCL"} object.
#'
#' @examples
#' p = 30
#' n = 50
#' beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
#' beta = c(beta, rep(0, times = p - length(beta)))
#' Comp_data = comp_Model(n = n, p = p, intercept = FALSE)
#' cvm1 <- cv.compCL(y = Comp_data$y, Z = Comp_data$X.comp,
#'                   Zc = Comp_data$Zc, intercept = Comp_data$intercept)
#' plot(cvm1)
#' plot(cvm1, xlab = "-log")
#'
#' @inherit coef.compCL author references
#'
#' @export


plot.cv.compCL <- function(x, xlab = c("log", "-log", "lambda"), trim = FALSE, ...) {
  cvobj <- x
  xlab <- match.arg(xlab)
  trim <- ifelse(trim, "Ttrim", "Ftrim")
  cvobj_use <- cvobj[[trim]]
  switch(xlab,
         "log" = {
           xlab = "Log(Lambda)"
           xvalue = log(drop(cvobj$lam))
         },
         "-log" = {
           xlab = "-Log(Lambda)"
           xvalue = -log(drop(cvobj$lam))
         },
         "lambda" = {
           xlab = "Lambda"
           xvalue = drop(cvobj$lam)
         })


  plot.args = list(x = xvalue,y = cvobj_use$cvm,
                   ylim = range(cvobj_use$cvupper,cvobj_use$cvlo),
                   xlab = xlab,
                   ylab = "MEAN-Squared Error",
                   type = "n")
  new.args=list(...)
  if(length(new.args)) plot.args[names(new.args)] = new.args
  do.call("plot",plot.args)
  error.bars(xvalue, cvobj_use$cvupper, cvobj_use$cvlo, width = 0.01, col="darkgrey")
  points(xvalue, cvobj_use$cvm, pch=20, col = "limegreen")
  axis(side = 3,at = xvalue,labels = paste(cvobj$compCL.fit$df),
       tick = FALSE, line=0, col.axis = "orange")
  # mtext("DF", side=3, line=3, cex.lab=2,las=0, col="orange")

  abline(v = switch(xlab,
                    "Log(Lambda)" = log(cvobj_use$lam.1se),
                    "-Log(Lambda)" = -log(cvobj_use$lam.1se),
                    "Lambda" = cvobj_use$lam.1se),
         lty=3, col = "blue")

  abline(v=switch(xlab,
                  "Log(Lambda)" = log(cvobj_use$lam.min),
                  "-Log(Lambda)" = -log(cvobj_use$lam.min),
                  "Lambda" = cvobj_use$lam.min),
         lty=3, col = "red")

}


#' @title
#' Plot the GIC curve produced by \code{"GIC.compCL"} object.
#'
#' @description
#' Plot the CIC curve as a function of the \code{lam} values.
#'
#' @param x fitted \code{"GIC.compCL"} object.
#'
#' @inheritParams plot.cv.compCL
#'
#' @details
#' A GIC curve is produced.
#'
#' @inherit plot.compCL return
#'
#' @seealso
#' \code{\link{GIC.compCL}} and \code{\link{compCL}}, and
#' \code{\link[=predict.GIC.compCL]{predict}} and
#' \code{\link[=coef.GIC.compCL]{coef}} methods for \code{"GIC.compCL"} object.
#'
#' @examples
#' p = 30
#' n = 50
#' beta = c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2)
#' beta = c(beta, rep(0, times = p - length(beta)))
#' Comp_data = comp_Model(n = n, p = p, beta = beta, intercept = FALSE)
#' GICm1 <- GIC.compCL(y = Comp_data$y, Z = Comp_data$X.comp,
#'                     Zc = Comp_data$Zc, intercept = Comp_data$intercept)
#' plot(GICm1)
#' plot(GICm1, xlab = "-log")
#'
#' @inherit coef.compCL author references
#'
#' @export


plot.GIC.compCL <- function(x, xlab = c("log", "-log", "lambda"), ...) {
  xlab <- match.arg(xlab)
  switch(xlab,
         "lambda" = {
           xlab = "Lambda"
           xvalue = drop(x$lam)
         },
         "log" = {
           xlab = "Log(Lambda)"
           xvalue = log(drop(x$lam))
         },
         "-log" = {
           xlab = "-Log(Lambda)"
           xvalue = -log(drop(x$lam))
         })


  plot.args = list(x = xvalue,y = x$GIC,
                   ylim = range(x$GIC),
                   xlab=xlab,
                   ylab="GIC",
                   type="n")
  new.args=list(...)
  if(length(new.args)) plot.args[names(new.args)]=new.args
  do.call("plot",plot.args)
  points(xvalue, x$GIC)
  axis(side = 3, at = xvalue, labels = paste(x$compCL.fit$df),
       tick = FALSE, line = 0, col.axis = "orange")
  # mtext("DF", side = 3, line = 3, cex.lab = 2, las = 0, col = "orange")

  abline(v=switch(xlab,
                  "Lambda" = x$lam.min,
                  "Log(Lambda)" = log(x$lam.min),
                  "-Log(Lambda)" = -log(x$lam.min)
  ), lty=3, col = "red")

}
#### <<< plot method >>> ####