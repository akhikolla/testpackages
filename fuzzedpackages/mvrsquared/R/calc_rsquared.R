
#' @title Calculate R-Squared.
#' @description
#'     Calculate R-Squared for univariate or multivariate outcomes.
#'
#' @param y The true outcome. This must be a numeric vector, numeric matrix, or
#'          coercible to a sparse matrix of class \code{dgCMatrix}. See 'Details'
#'          below for more information.
#' @param yhat The predicted outcome or a list of two matrices whose dot product
#'             makes the predicted outcome. See 'Details' below for more information.
#' @param ybar Numeric scalar or vector; the mean of \code{y}. Useful for parallel
#'             computation in batches.
#' @param return_ss_only Logical. Do you want to forego calculating R-squared and
#'                       only return the sums of squares?
#' @param threads Integer number of threads for parallelism; defaults to 1.
#' @return
#'     If \code{return_ss_only = FALSE}, \code{calc_rsqured} returns a numeric
#'     scalar R-squared. If \code{return_ss_only = TRUE}, \code{calc_rsqured}
#'     returns a vector; the first element is the error sum of squares (SSE) and
#'     the second element is the total sum of squares (SST). R-squared may then
#'     be calculated as \code{1 - SSE / SST}.
#'
#' @details
#'     There is some flexibility in what you can pass as \code{y} and \code{yhat}.
#'     In general, \code{y} can be a numeric vector, numeric matrix, a sparse
#'     matrix of class \code{dgCMatrix} from the \code{\link[Matrix]{Matrix}} package,
#'     or any object that can be coerced into a \code{dgCMatrix}.
#'
#'     \code{yhat} can be a numeric vector, numeric matrix, or a list of two
#'     matrices whose dot product has the same dimensionality as \code{y}. If
#'     \code{yhat} is a list of two matrices you may optionally name them \code{x}
#'     and \code{w} indicating the order of multiplication (\code{x} left
#'     multiplies \code{w}). If unnamed or ambiguously named, then it is assumed
#'     that \code{yhat[[1]]} left multiplies \code{yhat[[2]]}.
#'
#' @note
#'     On some Linux systems, setting \code{threads} greater than 1 for parallelism
#'     may introduce some imprecision in the calculation. As of this writing, the
#'     cause is still under investigation. In the meantime setting \code{threads = 1}
#'     should fix the issue.
#'
#'     Setting \code{return_ss_only} to \code{TRUE} is useful for parallel or
#'     distributed computing for large data sets, particularly when \code{y} is
#'     a large matrix. However if you do parallel execution you MUST pre-calculate
#'     'ybar' and pass it to the function. If you do not, SST will be calculated
#'     based on means of each batch independently. The resulting r-squared will
#'     be incorrect.
#'
#'     See example below for parallel computation with \code{\link[furrr]{future_map}}
#'     from the \code{furr} package.
#' @examples
#'
#' # standard r-squared with y and yhat as vectors
#' f <- stats::lm(mpg ~ cyl + disp + hp + wt, data = datasets::mtcars)
#'
#' y <- f$model$mpg
#'
#' yhat <- f$fitted.values
#'
#' calc_rsquared(y = y, yhat = yhat)
#'
#' # standard r-squared with y as a matrix and yhat containing 'x' and linear coefficients
#' s <- summary(f)
#'
#' x <- cbind(1, as.matrix(f$model[, -1]))
#'
#' w <- matrix(s$coefficients[, 1], ncol = 1)
#'
#' calc_rsquared(y = matrix(y, ncol = 1), yhat = list(x, w))
#'
#' # multivariate r-squared with y and yhat as matrices
#' calc_rsquared(y = cbind(y, y), yhat = cbind(yhat, yhat))
#'
#' # multivariate r-squared with yhat as a linear reconstruction of two matrices
#' calc_rsquared(y = cbind(y, y), yhat = list(x, cbind(w,w)))
#' @export
calc_rsquared <- function(
  y,
  yhat,
  ybar = NULL,
  return_ss_only = FALSE,
  threads = 1
) {

  ### Preliminary input checking ----
  if (! is.logical(return_ss_only) || is.na(return_ss_only))
    stop("'return_ss_only' must be logical TRUE/FALSE.")

  # is ybar null but return_ss_only not null?
  if (is.null(ybar) & return_ss_only)
    warning("'return_ss_only' is TRUE but 'ybar' is NULL. If you are calculating
            in batches, you may get a misleading result. If you are not, then
            disregard this warning.")

  ### Main part of function ----

  Y <- handle_y(y)

  Yhat <- handle_yhat(yhat, dim(Y))

  if (! is.null(ybar))
    if (length(ybar) != ncol(Y))
      stop("'ybar' is the wrong size. If 'ybar' is not NULL, length(ybar) must be
           the same as ncol(y) if 'y' is a matrix or length(ybar) must 1 if
           'y' is a vector.")

  if (is.null(ybar))
    ybar <- Matrix::colMeans(Y)

  result <- calc_sum_squares_latent(
    Y = Y,
    X = Yhat[[1]],
    W = Yhat[[2]],
    ybar = ybar,
    threads = threads
  )

  if (return_ss_only) {

    return(c(sse = result[1], sst = result[2]))

  } else {

    out <- 1 - result[1] / result[2]

    names(out) <- NULL

    return(out)

  }
}


# internal function checks and formats y
handle_y <- function(y) {

  # exception handling
  if ((! is.vector(y) && ! methods::canCoerce(y, "dgCMatrix")) || is.list(y))
    stop("'y' must be a numeric vector, numeric matrix, or be coercable to sparse
         matrix class 'dgCMatrix'. See 'help(calc_rsquared)' for more information")

  # get it formatted
  if (is.vector(y) && ! is.list(y)) {
    if (! is.numeric(y))
      stop("'y' does not appear to be numeric.")

    Y <- matrix(y, ncol = 1)

    Y <- methods::as(Y, "dgCMatrix")

  } else {

    if (! "dgCMatrix" %in% class(y)) {

      Y <- methods::as(y, "dgCMatrix")

    } else {

      Y <- y

    }

  }

  Y

}

# internal function checks and formats yhat
handle_yhat <- function(yhat, dim_y) {

  # exception handling
  if (! is.vector(yhat) && ! is.list(yhat) && ! is.matrix(yhat))
    stop("'yhat' must be a numeric vector, a numeric matrix, or a list of two
         numeric matrices whose dot product makes the predicted value of 'y'.")

  # does yhat have the right dimensions?
  if (is.vector(yhat) && ! is.list(yhat)) {

    if (dim_y[2] != 1)
      stop("'yhat' does not appear to have enough columns. 'yhat' must have the
           same dimensions as 'y' or 'yhat' must be a list of two numeric
           matrices whose dot product has the same dimensions as 'y'")

    if (dim_y[1] != length(yhat))
      stop("'yhat' does not appear to have enough observations. 'yhat' must have the
           same dimensions as 'y' or 'yhat' must be a list of two numeric
           matrices whose dot product has the same dimensions as 'y'")

  } else if (is.matrix(yhat)) {

    if (sum(dim(yhat) == dim_y) != 2) {

      msg <- paste("'y' and 'yhat' do not have the same dimensionality.",
                   "dim(y) is", paste(dim_y, collapse = " "),
                   "but dim(yhat) is", paste(dim(yhat), collapse = " "))

      stop(msg)

    }

  } else if (is.list(yhat)) {

    if (length(yhat) < 2)
      stop("'yhat' must be a numeric vector, a numeric matrix, or a list of two
         numeric matrices whose dot product makes the predicted value of 'y'.")

    if (sum(sapply(yhat[1:2], is.matrix)) != 2)
      stop("'yhat' appears to be a list but its first two elements are not matrices.")

    if (length(yhat) > 2)
      message("'yhat' is a list with more than two elements. Only the first two
              will be used.")

    # below allows you to pass named elements to the list that will get figured out
    if (sum(names(yhat) %in% c("x", "w")) > 0) {

      if (sum(names(yhat) %in% c("x", "w")) == 2) {

        yhat <- yhat[c("x", "w")]

      } else if (sum(names(yhat) %in% c("x", "w")) == 1) {

        warning("Found only one element of 'yhat' named 'x' or 'w'. Ignoring
                names of 'yhat' and using the first element as 'x' and the
                second elment as 'w'. You can explicitly name the elements of
                of 'yhat' as 'x' for the data matrix and 'w' for the weights.")

        yhat <- yhat[1:2]

      } else if (sum(names(yhat) %in% c("x", "w")) > 2) {
        warning("Found only multiple elements of 'yhat' named 'x' or 'w'. Ignoring
                names of 'yhat' and using the first element as 'x' and the
                second elment as 'w'. You can explicitly name the elements of
                of 'yhat' as 'x' for the data matrix and 'w' for the weights.")

        yhat <- yhat[1:2]

      }

    }

    if (sum(c(nrow(yhat[[1]]), ncol(yhat[[2]])) == dim_y) != 2 ||
        ncol(yhat[[1]]) != nrow(yhat[[2]]))
      stop("Dimensions of matrices in 'yhat' are not compatible. Either these
           matrices cannot be multiplied together or the result of their dot
           product does not have the same dimensions as 'y'.")



  }

  # I am explicitly not checking dimnames. I figure that it's too pedantic and
  # makes this code harder to read and troubleshoot. Best to make a note in the
  # docs or vignette(s) (or both) to encourage users to make sure they line up.
  # CAVEAT EMPTOR!

  # format yhat
  if (is.vector(yhat) && ! is.list(yhat)) {

    x <- matrix(yhat, ncol = 1)

    w <- matrix(1, ncol = 1, nrow = 1)

  } else if (is.matrix(yhat)) {

    x <- yhat

    w <- diag(ncol(yhat))

  } else { # must be a list. If it wasn't, would've gotten an error above

    x <- yhat[[1]]

    w <- yhat[[2]]

  }

  list(x = x, w = w)

}
