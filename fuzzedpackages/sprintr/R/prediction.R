# Calculate prediction from a \code{sprinter} object.
#
# Here I deliberately not make this as a S3 generic of \code{sprinter} object. Because there seems to be a discrepensy:
# (1) If I make it as S3 generic, then Roxygen2 will not export it, but just register it as a method of the S3 object. However, this behavious will make it hard if I want to explicitly use it somewhere.
# (2) If I manually export it by calling @export predict.sprinter, then this function is not registered as a S3 generic, which will also make a note in CRAN.
predict_sprinter <- function(object, newdata, ...) {
  # input check
  stopifnot(ncol(newdata) == object$p)
  idx <- object$idx[, , drop = FALSE]
  # selected indices for main effects
  idxm <- idx[idx[, 1] == 0, 2]
  # selected index pairs for interactions
  idxi <- idx[idx[, 1] != 0, , drop = FALSE]

  # need to standardize the main effects to construct interactions
  xm <- myscale(newdata)
  xint <- xm[, idxi[, 1]] * xm[, idxi[, 2]]

  fitted <- as.matrix(cbind(newdata[, idxm], xint) %*% object$coef)
  colnames(fitted) <- NULL
  fitted <- t(object$a0 + t(fitted))
  return(fitted)
}

#' Calculate prediction from a \code{cv.sprinter} object.
#'
#' @param object a fitted \code{cv.sprinter} object.
#' @param newdata a design matrix of all the \code{p} main effects of some new observations of which predictions are to be made.
#' @param ... additional argument (not used here, only for S3 generic/method consistency)
#' @return The prediction of \code{newdata} by the cv.sprinter fit \code{object}.
#' @examples
#' n <- 100
#' p <- 200
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] + 2 * x[, 2] - 3 * x[, 1] * x[, 2] + rnorm(n)
#' mod <- cv.sprinter(x = x, y = y)
#' fitted <- predict(mod, newdata = x)
#'
#' @export
predict.cv.sprinter <- function(object, newdata, ...) {
  # input check
  stopifnot(ncol(newdata) == object$p)
  idx <- object$compact[, 1:2, drop = FALSE]
  # selected indices for main effects
  idxm <- idx[idx[, 1] == 0, 2]
  # selected index pairs for interactions
  idxi <- idx[idx[, 1] != 0, , drop = FALSE]

  # need to standardize the main effects to construct interactions
  xm <- myscale(newdata)
  xint <- xm[, idxi[, 1]] * xm[, idxi[, 2]]

  fitted <- as.numeric(object$a0 + cbind(newdata[, idxm], xint) %*% object$compact[, 3])
  return(fitted)
}
