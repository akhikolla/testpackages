#' @title
#' Fit regularization paths of sparse log-contrast regression with functional
#' compositional predictors.
#'
#' @description
#' Fit the penalized \emph{log-contrast} regression with functional compositional predictors
#' proposed by Zhe et al. (2020) <arXiv:1808.02403>. The model estimation is
#' conducted by minimizing a linearly constrained group lasso criterion. The regularization
#' paths are computed for the group lasso penalty at grid values of the regularization
#' parameter \code{lam} and the degree of freedom of the basis function \code{K}.
#'
#' @usage
#' FuncompCGL(y, X, Zc = NULL, intercept = TRUE, ref = NULL,
#'            k, degree = 3, basis_fun = c("bs", "OBasis", "fourier"),
#'            insert = c("FALSE", "X", "basis"), method = c("trapezoidal", "step"),
#'            interval = c("Original", "Standard"), Trange,
#'            T.name = "TIME", ID.name = "Subject_ID",
#'            W = rep(1,times = p - length(ref)),
#'            dfmax = p - length(ref), pfmax = min(dfmax * 1.5, p - length(ref)),
#'            lam = NULL, nlam = 100, lambda.factor = ifelse(n < p1, 0.05, 0.001),
#'            tol = 1e-8, mu_ratio = 1.01,
#'            outer_maxiter = 1e+6, outer_eps = 1e-8,
#'            inner_maxiter = 1e+4, inner_eps = 1e-8)
#
#' @param y response vector with length n.
#'
#' @param X data frame or matrix.
#'       \itemize{
#'       \item If \code{nrow(X)} > \eqn{n},
#'             \code{X} should be a data frame or matrix of the functional compositional
#'             predictors with \eqn{p} columns for the values of the composition components,
#'             a column indicating subject ID and a column of observation times.
#'             Order of Subject ID should be the SAME as that of \code{y}.
#'             Zero entry is not allowed.
#'       \item If \code{nrow(X)[1]}=\eqn{n},
#'             \code{X} is considered as after taken integration, a
#'             \code{n*(k*p - length(ref))} matrix.
#'       }
#'
#' @param T.name,ID.name a character string specifying names of the time variable and the Subject
#'                       ID variable in \code{X}.
#'                       This is only needed when X is a data frame or matrix of the functional
#'                       compositional predictors.
#'                       Default are \code{"TIME"} and \code{"Subject_ID"}.
#'
#' @param Zc a \eqn{n\times p_{c}}{n*p_c} design matrix of unpenalized variables.
#'           Default is NULL.
#'
#' @param k an integer, degrees of freedom of the basis function.
#' @param ref reference level (baseline), either an integer between \eqn{[1,p]} or \code{NULL}.
#'            Default value is \code{NULL}.
#'            \itemize{
#'            \item If \code{ref} is set to be an integer between \code{[1,p]}, the group lasso penalized \emph{log-contrast} model (with log-ratios) is fitted
#'            with the \code{ref}-th component chosed as baseline.
#'            \item If \code{ref} is set to be \code{NULL}, the linearly constrained group lasso penalized \emph{log-contrast} model is fitted.
#'            }
#'
#' @param degree degrees of freedom of the basis function. Default value is 3.
#'
#' @param basis_fun method to generate basis:
#'        \itemize{
#'        \item \code{"bs"}(Default) B-splines. See fucntion \code{\link[splines]{bs}}.
#'        \item \code{"OBasis"} Orthoganal B-splines. See function \code{\link[orthogonalsplinebasis]{OBasis}}
#'              and package \pkg{orthogonalsplinebasis}.
#'        \item \code{"fourier"} Fourier basis. See fucntion \code{\link[fda]{create.fourier.basis}}
#'              and package \pkg{fda}.
#'        }
#'
#' @param insert a character string sepcifying method to perform functional interpolation.
#'               \itemize{
#'               \item \code{"FALSE"}(Default) no interpolation.
#'               \item \code{"X"} linear interpolation of functional compositional
#'                     data along the time grid.
#'               \item \code{"basis"} the functional compositional data is interplolated
#'                     as a step function along the time grid.
#'               }
#'               If \code{insert} = \code{"X"} or \code{"basis"}, interplolation is conducted
#'               on \code{sseq}, where \code{sseq} is the sorted sequence of all the observed time points.
#                Denser time sequence is generated,
#                equally spaced by \code{min(diff(sseq))/20)},
#                where \code{sseq} is the sorted sequence of all the observed time points.
#'
#' @param method a character string sepcifying method used to approximate integral.
#'               \itemize{
#'               \item \code{"trapezoidal"}(Default) Sum up areas under the trapezoids.
#                     See \code{\link{ITG_trap}}.
#'               \item \code{"step"} Sum up area under the rectangles.
#                     See \code{\link{ITG_step}}.
#'               }
#'
#' @param interval a character string sepcifying the domain of the integral.
#'        \itemize{
#'          \item \code{"Original"}(Default) On the original time scale,
#'                \code{interval} = \code{range(Time)}.
#'          \item \code{"Standard"} Time points are mapped onto \eqn{[0,1]},
#'                \code{interval} = \code{c(0,1)}.
#'        }
#'
#' @param W a vector of length p (the total number of groups),
#'          or a matrix with dimension \eqn{p_{1} \times p_{1}}{p1*p1}, where
#'          \code{p1=(p - length(ref)) * k},
#'          or character specifying the function used to calculate weight matrix for each group.
#'          \itemize{
#'          \item a vector of penalization weights for the groups of coefficients. A zero weight implies no shrinkage.
#'          \item a diagonal matrix with positive diagonal elements.
#'          \item if character string of function name or an object of type \code{function} to compute the weights.
#'          }
#'
#'
#' @param Trange range of time points
#'
#' @param outer_maxiter,outer_eps \code{outer_maxiter} is the maximum number of loops allowed for the augmented Lanrange method;
#'                                and \code{outer_eps} is the corresponding convergence tolerance.
#'
#' @param inner_maxiter,inner_eps \code{inner_maxiter} is the maximum number of loops allowed for blockwise-GMD;
#'                                and \code{inner_eps} is the corresponding convergence tolerance.
#'
#' @inheritParams cglasso
#'
#' @return An object with S3 class \code{"FuncompCGL"}, which is a list containing:
#' \item{Z}{the integral matrix for the functional compositional predictors with
#'          dimension \eqn{n \times (pk)}{n*(pk)}.}
#' \item{lam}{the sequence of \code{lam} values.}
#' \item{df}{the number of non-zero groups in the estimated coefficients for
#'           the functional compositional predictors at each value of \code{lam}.}
#' \item{beta}{a matrix of coefficients with \code{length(lam)} columns and
#'             \eqn{p_{1}+p_{c}+1}{p_1+p_c+1} rows, where \code{p_1=p*k}.
#'             The first \eqn{p_{1}}{p_1} rows are the estimated values for
#'             the coefficients for the functional compositional preditors, and the
#'             last row is for the intercept. If \code{intercept = FALSE}, the last row is 0's.}
#' \item{dim}{dimension of the coefficient matrix.}
#' \item{sseq}{sequence of the time points.}
#' \item{call}{the call that produces this object.}
#'
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
#' print(m1)
#' plot(m1)
#' beta <- coef(m1)
#' arg_list <- as.list(Data$call)[-1]
#' arg_list$n <- 30
#' TEST <- do.call(Fcomp_Model, arg_list)
#' y_hat <- predict(m1, Znew = TEST$data$Comp, Zcnew = TEST$data$Zc)
#' plot(y_hat[, floor(length(m1$lam)/2)], TEST$data$y,
#'      ylab = "Observed Response", xlab = "Predicted Response")
#'
#' beta <- coef(m1, s = m1$lam[20])
#' beta_C <- matrix(beta[1:(p*df_beta)], nrow = p, byrow = TRUE)
#' colSums(beta_C)
#' Non.zero <- (1:p)[apply(beta_C, 1, function(x) max(abs(x)) > 0)]
#' Non.zero
#'
#' @seealso
#' \code{\link{cv.FuncompCGL}} and \code{\link{GIC.FuncompCGL}}, and
#' \code{\link[=predict.FuncompCGL]{predict}}, \code{\link[=coef.FuncompCGL]{coef}},
#' \code{\link[=plot.FuncompCGL]{plot}} and \code{\link[=print.FuncompCGL]{print}}
#' methods for \code{"FuncompCGL"} object.
#'
#'
#' @references
#' Sun, Z., Xu, W., Cong, X., Li G. and Chen K. (2020) \emph{Log-contrast regression with
#' functional compositional predictors: linking preterm infant's gut microbiome trajectories
#' to neurobehavioral outcome}, \href{https://arxiv.org/abs/1808.02403}{https://arxiv.org/abs/1808.02403}
#' \emph{Annals of Applied Statistics}.
#'
#' Yang, Y. and Zou, H. (2015) \emph{A fast unified algorithm for computing
#' group-lasso penalized learning problems}, \href{https://link.springer.com/article/10.1007/s11222-014-9498-5}{
#' https://link.springer.com/article/10.1007/s11222-014-9498-5}
#' \emph{Statistics and Computing} \strong{25(6)} 1129-1141.
#'
#' Aitchison, J. and Bacon-Shone, J. (1984) \emph{Log-contrast models for experiments with mixtures},
#  \href{https://www.jstor.org/stable/pdf/2336249.pdf}{https://www.jstor.org/stable/pdf/2336249.pdf}
#' \emph{Biometrika} \strong{71} 323-330.
#'
#' @author
#' Zhe Sun and Kun Chen
#'
#' @details
#' The \emph{functional log-contrast regression model}
#' for compositional predictors is defined as
#' \deqn{
#' y = 1_n\beta_0 + Z_c\beta_c + \int_T Z(t)\beta(t)dt + e, s.t. (1_p)^T \beta(t)=0  \forall t \in T,
#' }
#' where \eqn{\beta_0} is the intercept,
#' \eqn{\beta_c} is the regression coefficient vector with length \eqn{p_c} corresponding to the control variables,
#' \eqn{\beta(t)} is the functional regression coefficient vector with length \eqn{p} as a funtion of \eqn{t}
#' and \eqn{e} is the random error vector with zero mean with length \eqn{n}.
#' Moreover, \code{Z(t)} is the log-transformed functional compostional data.
#' If zero(s) exists in the original functional compositional data, user should pre-process these zero(s).
#' For example, if count data provided, user could replace 0's with 0.5.
#' \cr
#' After adopting a truncated basis expansion approach to re-express \eqn{\beta(t)}
#' \deqn{\beta(t) = B \Phi(t),}
#' where \eqn{B} is a p-by-k unkown but fixed coefficient matrix, and \eqn{\Phi(t)} consists of basis with
#' degree of freedom \eqn{k}. We could write \emph{functional log-contrast regression model} as
#' \deqn{
#' y = 1_n\beta_0 + Z_c\beta_c + Z\beta + e, s.t. \sum_{j=1}^{p}\beta_j=0_k,
#' }
#' where \eqn{Z} is a n-by-pk matrix corresponding to the integral, \eqn{\beta=vec(B^T)} is a
#' pk-vector with every each k-subvector corresponding to the coefficient vector for the \eqn{j}-th
#' compositional component.
#' \cr
#' To enable variable selection, \code{FuncompCGL} model is estimated via linearly constrained group lasso,
#' \deqn{
#' argmin_{\beta_0, \beta_c, \beta}(\frac{1}{2n}\|y - 1_n\beta_0 - Z_c\beta_c - Z\beta\|_2^2
#' + \lambda \sum_{j=1}^{p} \|\beta_j\|_2), s.t. \sum_{j=1}^{p} \beta_j = 0_k.
#' }
#'
#' @export
#'
#' @import stats
#' @importFrom splines bs
#' @importFrom fda eval.basis create.fourier.basis
#' @importFrom orthogonalsplinebasis evaluate OBasis expand.knots


FuncompCGL <- function(y, X, Zc = NULL, intercept = TRUE, ref = NULL,
                       k, degree = 3, basis_fun = c("bs", "OBasis", "fourier"),
                       insert = c("FALSE", "X", "basis"), method = c("trapezoidal", "step"),
                       interval = c("Original", "Standard"), Trange,
                       T.name = "TIME", ID.name = "Subject_ID",
                       W = rep(1,times = p - length(ref)),
                       dfmax = p - length(ref), pfmax = min(dfmax * 1.5, p - length(ref)),
                       lam = NULL, nlam = 100, lambda.factor = ifelse(n < p1, 0.05, 0.001),
                       tol = 1e-8, mu_ratio = 1.01,
                       outer_maxiter = 1e+6, outer_eps = 1e-8,
                       inner_maxiter = 1e+4, inner_eps = 1e-8) {

  n <- length(y)
  this.call <- match.call()
  basis_fun <- match.arg(basis_fun)
  method <- match.arg(method)
  insert <- match.arg(insert)
  interval <- match.arg(interval)
  interval_choice <- interval

  if(missing(k)) stop("Value(s) of k must be provided.")
  if(basis_fun %in% c("bs", "OBasis") && k <= degree) stop("df k is too small, must be larger than provided degree.")
  if(basis_fun == "fourier" && k%%2 == 0) stop("df k must be an odd integer")
  if(!is.null(lam) && TRUE %in% (lam < 0)) stop("User provided lambda must be positive vector")

  digits = 10

  if(dim(X)[1] == n) {
    # Case I: already take integral and ready to go
    Z <-  as.matrix(X)
    p1 <- ncol(Z)
    p <- p1 / k + length(ref)
    if(p %% 1 != 0) stop("Wrong dimension of k.")
    if(is.numeric(ref)) {
      if( !(ref %in% 1:p) ) stop("Reference level is out of range")
      mu_ratio = 0 # outer_maxiter will be set to 1 automatically
    }
    if(!missing(Trange)) Trange_choice <- Trange else Trange_choice <- NULL
    sseq <- NULL
  } else {
    ## Case II: take integral
    X <- as.data.frame(X)
    X.names <- colnames(X)[!colnames(X) %in% c(T.name, ID.name)]
    p <- length(X.names)

    ## transform X into percentage and take log
    if( any(X[, X.names] == 0) ) stop("There is entry with value 0")
    X[, X.names] <- log(X[, X.names] / rowSums(X[, X.names]))
    # X[, X.names] <- proc.comp(as.matrix(X[, X.names]))

    if(is.numeric(ref)) {
      if( !(ref %in% 1:p) ) stop("Reference level is out of range")
      mu_ratio = 0 # outer_maxiter will be set to 1 automatically
    }

    # Time range provided in case that sample do not cover the entire range
    if(missing(Trange)) {
      Trange <- range(X[, T.name])
    } else {
      filter_id <- X[, T.name] >= Trange[1]
      filter_id <- X[filter_id, T.name] <= Trange[2]
      X <- X[filter_id, ]
    }
    if(nrow(X) == 0) stop("Trange is out of range of observed data points.")
    Trange_choice <- Trange

    switch(interval,
           "Standard" = {
              ## mapping time sequence on to [0,1]
              X[, T.name] = (X[, T.name] - min(Trange)) / diff(Trange)
              interval = c(0, 1)
              Trange = c(0, 1)
           },
           "Original" = {
              interval = Trange
           })


    X[, T.name] = round(X[, T.name], digits = digits)
    ## In case that X is not represented in numerical order of Subject_ID
    X[, ID.name] <- factor(X[, ID.name], levels = unique(X[, ID.name]))

    ## Take log-ratio w.r.t. reference
    X[, X.names] <- apply(X[, X.names], 2, function(x, ref) x - ref,
                          ref = if(is.null(ref)) 0 else X[, X.names[ref], drop = TRUE])

    ## If ref is NULL, take the p+1-th component off (nothing);
    ## Ortherwise take the ref-th component off
    D <- split(X[, c(T.name, X.names[-ifelse(is.null(ref), p + 1, ref)])], X[, ID.name])
    p1 <- (p - length(ref)) * k




    ## Generate discrete basis matrix
    sseq <- sort(unique(c(round(Trange, digits = digits), as.vector(X[, T.name]))))

    # if(insert != "FALSE") sseq <- round(seq(from = interval[1], to = interval[2],
    #                                         by = min(diff(sseq))/10 # by = length(sseq) * 2
    #                                         ),
    #                                     digits = digits)




    basis <- switch(basis_fun,
                    "bs" = bs(x = sseq, df = k, degree = degree,
                              Boundary.knots = interval, intercept = TRUE),
                    "fourier" = eval.basis(sseq,
                                           basisobj = create.fourier.basis(rangeval = interval, nbasis = k )),
                    "OBasis" = {
                      ## Generate knots equally
                      nknots <- k - (degree + 1)
                      if(nknots > 0) {
                        knots <- ((1:nknots) / (nknots + 1)) * diff(interval) + interval[1]
                      } else  knots <- NULL

                      evaluate(OBasis(expand.knots(c(interval[1], knots, interval[2])),
                                                   order = degree + 1),
                               sseq)
                      }
                    )

    Z <- matrix(NA, nrow = n, ncol = p1)
    for(i in 1:n) Z[i, ] <- ITG(D[[i]], basis, sseq, T.name, interval, insert, method)$integral

  }



  group.index <- matrix(1:p1, nrow = k)
  if( is.vector(W) ) {
    if(length(W) != (p1 / k) ) stop("W shoulde be a vector of length p=", p1/ k)
  } else if( is.matrix(W) ) {
    if(any(dim(W) - c(p1, p1) != 0)) stop('Wrong dimensions of Weights matrix')
  } else if( is.function(W) || is.character(W) ){
    if(is.character(W)) W <- get(W)
    W_fun <- W
    W <- diag(p1)
    for(i in 1:(p1/k)) W[group.index[, i], group.index[, i]] <- W_fun(Z[, group.index[, i]])
  }

  output <- cglasso(y = y, Z = Z, Zc = Zc, k = k, W = W, intercept = intercept,
                    mu_ratio = mu_ratio,
                    lam = lam, nlam = nlam, lambda.factor = lambda.factor,
                    dfmax = dfmax, pfmax = pfmax,
                    tol = tol,
                    outer_maxiter = outer_maxiter, outer_eps = outer_eps,
                    inner_maxiter = inner_maxiter, inner_eps = inner_eps)
  output$Z <- Z
  output$W <- W
  output$call <- as.list(this.call)
  output$call['basis_fun'] = basis_fun
  output$call['method'] = method
  output$call['insert'] = insert
  output$call['interval'] = interval_choice
  output$call['degree'] = degree
  output$call$Trange = Trange_choice # c(Trange_choice[1], Trange_choice[2])
  output$call[['k']] = k
  output$call <- as.call(output$call)
  output$ref <- ref
  output$sseq <- sseq
  if(is.numeric(ref)) output$df[output$df != 0] <- output$df[output$df != 0] + 1
  class(output) <- "FuncompCGL"
  output$dim <- dim(output$beta)
  return(output)
}