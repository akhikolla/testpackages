
tpscov <- function (r, m = 2, d = 1)
{
  r <- as.matrix(r)
  num.row <- nrow(r)
  num.col <- ncol(r)
  r <- as.vector(r)
  nzi <- (1:length(r))[r != 0]
  ans <- rep(0, length(r))
  if ((d + 1)%%2 != 0)
    ans[nzi] <- (abs(r[nzi]))^(2 * m - d) * log(abs(r[nzi]))
  else ans[nzi] <- (abs(r[nzi]))^(2 * m - d)
  if (num.col > 1)
    ans <- matrix(ans, num.row, num.col)
  return(ans)
}

#' Creates a design matrix from a bivariate smoothing algorithm
#'
#' \code{create_bivariate_design} accepts two numeric vectors of equal length as inputs. From these inputs, a bivariate smoothing design matrix is produced using thin plate splines. 
#'
#' @export
#' @param X1 numeric vector for first variable
#' @param X2 numeric vector for second variable
#' @param num_knots optional: number of knots
#' @param knots optional: matrix of knot locations for bivariate smoothing
#' @references Ruppert, David, Matt P. Wand, and Raymond J. Carroll. \emph{Semiparametric Regression}. No. 12. Cambridge university press, 2003. Section 13.5
#' @references Matt Wand (2018). SemiPar: Semiparametric Regression. R package version 1.0-4.2.
#' @return list containing the design matrix \code{Z} and matrix of \code{knots}
#' @examples
#' x1 <- rnorm(100)
#' x2 <- rnorm(100)
#' res <- create_bivariate_design(x1, x2)
#' res$knots
#' dim(res$Z)
create_bivariate_design <- function(X1, X2, num_knots = NULL, knots = NULL) {

  X_bivariate <- cbind(X1, X2)

  if (is.null(knots)) {
    if (is.null(num_knots))
      num_knots <- max(10, min(50, round(nrow(X_bivariate)/4)))
    knots <- cluster::clara(X_bivariate, num_knots)$medoids
  } else {
    num_knots <- nrow(knots)
  }

  dist_mat <- matrix(0, num_knots, num_knots)

  # euclidean distance between knots
  dist_mat[lower.tri(dist_mat)] <- dist(as.matrix(knots))
  dist_mat <- dist_mat + t(dist_mat)

  Omega <- tpscov(dist_mat, m = 2, d = 2)
  x.knots.diff1 <- outer(X_bivariate[, 1], knots[, 1], "-")
  x.knots.diff2 <- outer(X_bivariate[, 2], knots[, 2], "-")
  x.knots.dists <- sqrt(x.knots.diff1^2 + x.knots.diff2^2)
  prelim.Z <- tpscov(x.knots.dists, m = 2, d = 2)

  # square root of omega
  sva <- svd(Omega)
  sqrt.Omega <- t(sva$v %*% (t(sva$u) * sqrt(sva$d)))
  trans.mat <- list(sqrt.Omega)

  # random effects design matrix
  retval <- list(Z = t(solve(sqrt.Omega, t(prelim.Z))),
                 knots = knots)
  return(retval)

}

# helper function for dataframe case
create_bivariate_design_dat <- function(X, num_knots = NULL, knots = NULL) {
  create_bivariate_design(X1=X[, 1], X2=X[, 2], num_knots=num_knots, knots=knots)
}

select_knots <- function(x, num_knots=NULL, ...) {
  if (is.null(num_knots)) {
    num_knots <- pmin(0.25* length(unique(x)), 25)
  }
  quantile(unique(x), seq(0,1,length=(num_knots+2))[-c(1,(num_knots+2))], ...)
}

#' Creates design matrices for univariate and bivariate applications
#' 
#' \code{np} accepts one or two numeric vectors of equal length as inputs. From these inputs, univariate or bivariate smoothing design matrices are produced. Currently available basis functions are truncated polynomials and thin plate splines.  
#' When bivariate smoothing is selected, \code{np} calls \code{\link{create_bivariate_design}}. 
#' @export
#' @param x1 numeric vector
#' @param x2 optional vector for bivariate non-parametric function
#' @param num_knots optional number of knots
#' @param knots optional numeric vector of knots
#' @param basis character vector for basis function.  \code{tps} for thin-plate spline and \code{trunc.poly} for truncated polynomial
#' @param degree for truncated polynomial basis function
#' @references Ruppert, David, Matt P. Wand, and Raymond J. Carroll. \emph{Semiparametric Regression}. No. 12. Cambridge university press, 2003. Section 5.6.
#' @references Matt Wand (2018). SemiPar: Semiparametric Regression. R package version 1.0-4.2.
#' @return list with the following elements:
#' \itemize{
#' \item \code{X} parametric design matrix
#' \item \code{Z} non-parametric design matrix
#' \item \code{knots} numeric vector of knots for the model
#' \item \code{Xnms} names of parameters passed to `np`
#' \item \code{basis} selected basis function
#' \item \code{degree} degree for truncated polynomial basis function
#' }
#' @examples
#' x1 <- rnorm(100)
#' res <- np(x1, num_knots=10, basis="trunc.poly", degree=2)
#' res
np <- function(x1, x2=NULL, num_knots=NULL, knots=NULL, basis="tps", degree=3) {

  # get names passed to x1 x2
  arg1 <- deparse(substitute(x1))
  arg2 <- NULL
  
  # set basis to tps if two variables
  if (!is.null(x2)) {
    arg2 <- deparse(substitute(x2))
    basis <- "tps"
  }
  
  # Fixed design
  # select knots if not specified
  if (is.null(x2)) {
    
    if (is.null(knots)) {
      knots <- select_knots(x=x1, num_knots=num_knots, na.rm=TRUE)
    }

    # univariate
    # truncated polynomial basis
    if (basis == "trunc.poly") {
      npow <- 1:degree
      X_tp <- matrix(x1^(rep(npow, each=length(x1))),ncol=degree,byrow=FALSE)
      colnames(X_tp) <- rep("", ncol(X_tp))

      Z <- outer(x1, knots, "-")
      Z[Z < 0] <- 0
      Z <- abs(Z)^degree
      Z_tp <- Z
    } else if (basis == "tps") {
      X_tp <- matrix(x1)
      
      # thin-plated spline
      Z_K<-(abs(outer(x1, knots,"-")))^degree
      OMEGA_all<-(abs(outer(knots, knots, "-")))^degree
      svd.OMEGA_all<-svd(OMEGA_all)
      sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*%
                          (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d)))
      Z_tp <-t(solve(sqrt.OMEGA_all,t(Z_K)))
    } else {
      stop("basis must by trunc.poly or tps")
    }
  } else {
    # bivariate
    X_tp <- cbind(x1, x2)
    biv <- create_bivariate_design(x1, x2, num_knots)
    Z_tp <- biv$Z
    knots <- biv$knots
  }
  
  return(list(X=X_tp,
              Z=Z_tp,
              knots=knots,
              Xnms = c(arg1, arg2),
              basis=basis,
              degree=degree))
}

