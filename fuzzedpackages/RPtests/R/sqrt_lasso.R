#' Square-root Lasso regression
#'
#' Fits a linear model to potentially high-dimensional data using the
#' square-root Lasso, also known as the scaled Lasso. The Lasso path is computed
#' using the \pkg{glmnet} package.
#'
#' @param x Input matrix of dimension nobs by nvars; each row is an observation
#'   vector.
#' @param y Response variable; shoud be a numeric vector.
#' @param lam0 Tuning parameter for the square-root / scaled Lasso. If left
#'   blank (recommended) this is chosen using the method of Sun & Zhang (2013)
#'   implemented in the \pkg{scalreg} package.
#' @param exclude Indices of variables to be excluded from the model; default is
#'   none.
#' @param output_all In addition to the vector of coefficients, if \code{TRUE},
#'   also outputs the intercept, an estimate of the noise standard deviation,
#'   and the output of \code{\link[glmnet]{glmnet}}.
#' @param ... Additional arguments to be passed to \code{\link[glmnet]{glmnet}}.
#' @details First the Lasso path is computed using \code{glmnet} from
#'   \pkg{glmnet}. Next the particular point on the path corresponding to the
#'   square-root Lasso solution is found. As the path is only computed on a grid
#'   of points, the square-root solution is approximate.
#' @return Either an estimated vector of regression coefficients with nvars
#'   components or, if \code{output_all} is \code{true}, a list with components
#'   \describe{
#'     \item{\code{beta}}{the vector of regression coefficents}
#'     \item{\code{a0}}{an intercept term}
#'     \item{\code{sigma_hat}}{an estimate of
#'       the noise standard deviation; this is calculated as square-root of the
#'       average residual sums of squares}
#'     \item{\code{glm_obj}}{the fitted \code{glmnet} object, an S3 class ``\code{glmnet}"}
#'   }
#' @references
#'   A. Belloni, V. Chernozhukov, and L. Wang. (2011)
#'   \emph{Square-root lasso: pivotal recovery of sparse signals via conic
#'   programming. Biometrika, 98(4):791-806.}
#'   \url{http://biomet.oxfordjournals.org/content/98/4/791.refs} T. Sun and
#'   C.-H. Zhang. (2012) \emph{Scaled sparse linear regression. Biometrika,
#'   99(4):879-898.}
#'   \url{http://biomet.oxfordjournals.org/content/early/2012/09/24/biomet.ass043.short}
#'    T. Sun and C.-H. Zhang. (2013) \emph{Sparse matrix inversion with scaled
#'   lasso. The Journal of Machine Learning Research, 14(1):3385-3418.}
#'   \url{www.jmlr.org/papers/volume14/sun13a/sun13a.pdf}
#' @seealso \code{\link[glmnet]{glmnet}} and \code{\link[scalreg]{scalreg}}.
#' @examples
#' x <- matrix(rnorm(100*250), 100, 250)
#' y <- x[, 1] + x[, 2] + rnorm(100)
#' out <- sqrt_lasso(x, y)
#' @importFrom glmnet glmnet predict.glmnet
#' @export
sqrt_lasso <- function(x, y, lam0=NULL, exclude = integer(0), output_all = FALSE, ...) {
  # Do this first so error checking can take place
  out <- glmnet(x, y, exclude=exclude, ...)

  n <- nrow(x)
  p <- ncol(x) - length(exclude)

  # lam0 copied from scalreg
  if (is.null(lam0)) {
    if (p == 1) {
      L = 0.5
    } else {
      L = 0.1
      Lold = 0
      while (abs(L - Lold) > 0.001) {
        k = (L^4 + 2 * L^2)
        Lold = L
        L = -qnorm(min(k/p, 0.99))
        L = (L + Lold)/2
      }
    }
    lam0 = sqrt(2/n) * L
  }

  resids <- y - predict.glmnet(out, newx=x)
  full_MSE <- colMeans(resids^2)
  index_sel <- which.min(abs(full_MSE*lam0^2 / out$lambda^2 - 1))
  if (index_sel == length(full_MSE)) warning("Smallest lambda chosen by sqrt_lasso")
  if (output_all) {
    return(list("beta"=as.numeric(out$beta[, index_sel]),
                "a0"=out$a0[index_sel],
                "sigma_hat"=sqrt(full_MSE[index_sel]),
                "glmnet_obj"=out))
  }
  return(as.numeric(out$beta[, index_sel]))
}
