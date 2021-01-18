#' Running sprinter with cross-validation
#'
#' The main cross-validation function to select the best sprinter fit for a path of tuning parameters.
#'
#' @param x An \code{n} by \code{p} design matrix of main effects. Each row is an observation of \code{p} main effects.
#' @param y A response vector of size \code{n}.
#' @param num_keep Number of candidate interactions to keep in Step 2. If \code{num_keep} is not specified (as default), it will be set to \code{[n / log n]}.
#' @param square Indicator of whether squared effects should be fitted in Step 1. Default to be FALSE.
#' @param lambda A user specified list of tuning parameter. Default to be NULL, and the program will compute its own \code{lambda} path based on \code{nlam} and \code{lam_min_ratio}.
#' @param nlam The number of \code{lambda} values. Default value is \code{100}.
#' @param lam_min_ratio The ratio of the smallest and the largest values in \code{lambda}. The largest value in \code{lambda} is usually the smallest value for which all coefficients are set to zero. Default to be \code{1e-2} in the \code{n} < \code{p} setting.
#' @param nfold Number of folds in cross-validation. Default value is 5. If each fold gets too view observation, a warning is thrown and the minimal \code{nfold = 3} is used.
#' @param foldid A vector of length \code{n} representing which fold each observation belongs to. Default to be \code{NULL}, and the program will generate its own randomly.
#' @return An object of S3 class "\code{sprinter}".
#'  \describe{
#'   \item{\code{n}}{The sample size.}
#'   \item{\code{p}}{The number of main effects.}
#'   \item{\code{a0}}{estimate of intercept corresponding to the CV-selected model.}
#'   \item{\code{compact}}{A compact representation of the selected variables. \code{compact} has three columns, with the first two columns representing the indices of a selected variable (main effects with first index = 0), and the last column representing the estimate of coefficients.}
#'   \item{\code{fit}}{The whole \code{glmnet} fit object in Step 3.}
#'   \item{\code{fitted}}{fitted value of response corresponding to the CV-selected model.}
#'   \item{\code{lambda}}{The sequence of \code{lambda} values used.}
#'   \item{\code{cvm}}{The averaged estimated prediction error on the test sets over K folds.}
#'   \item{\code{cvsd}}{The standard error of the estimated prediction error on the test sets over K folds.}
#'   \item{\code{foldid}}{Fold assignment. A vector of length \code{n}.}
#'   \item{\code{ibest}}{The index in \code{lambda} that is chosen by CV.}
#'   \item{\code{call}}{Function call.}
#'  }
#' @seealso
#'   \code{\link{predict.cv.sprinter}}
#' @examples
#' n <- 100
#' p <- 200
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5] + rnorm(n)
#' mod <- cv.sprinter(x = x, y = y)
#'
#' @import glmnet
#' @export
cv.sprinter <- function(x, y, num_keep = NULL, square = FALSE,
                        lambda = NULL, nlam = 100,
                        lam_min_ratio = ifelse(nrow(x) < ncol(x), 0.01, 1e-04),
                        nfold = 5, foldid = NULL){
  n <- nrow(x)
  p <- ncol(x)
  stopifnot(n == length(y))

  # by default, num_keep is set to be [n / log(n)]
  if(is.null(num_keep)){
    num_keep = ceiling(n / log(n))
  }
  else{
    num_keep = min(ceiling(num_keep), ifelse(square, p*(p - 1) / 2, p * (p - 1)/ 2 + p))
  }

  # first fit the sprinter using all data with all lambdas
  fit <- sprinter(x = x, y = y, num_keep = num_keep, square = square, lambda = lambda, nlam = nlam, lam_min_ratio = lam_min_ratio)
#  cat("cv initial finished!", fill = TRUE)

  colnames(fit$coef) <- NULL
  rownames(fit$coef) <- NULL
  if(is.null(lambda)){
    # if lambda is not provided
    lambda <- fit$lambda
  }
  nlam <- length(lambda)

  # use cross-validation to select the best lambda
  if (is.null(foldid)){
    # foldid is a vector of values between 1 and nfold
    # identifying what fold each observation is in.
    # If supplied, nfold can be missing.
    foldid <- sample(seq(nfold), size = n, replace = TRUE)
  }

  # mse of lasso estimate of coef
  mat_mse <- matrix(NA, nrow = nlam, ncol = nfold)
  for (i in seq(nfold)){
    # train on all but i-th fold
    id_tr <- (foldid != i)
    id_te <- (foldid == i)
    # training/testing data set
    # standardize the training data
    x_tr <- myscale(x[id_tr, ])
    # use the scaling for training data
    x_te <- myscale(x[id_te, ],
                  center = attr(x = x_tr, which = "scaled:center"),
                  scale = attr(x = x_tr, which = "scaled:scale"))
    y_tr <- y[id_tr]
    y_te <- y[id_te]

    # get the fit using lasso on training data
    # and fit/obj on the test data
    fit_tr <- sprinter(x = x_tr, y = y_tr, num_keep = num_keep, square = square, lambda = lambda, nlam = nlam, lam_min_ratio = lam_min_ratio)
    pred_te <- predict_sprinter(object = fit_tr, newdata = x_te)
    mat_mse[, i] <- sqrt(as.numeric(colMeans((y_te - pred_te)^2)))
#    cat(paste("cv: fold ", i, " finished!"), fill = TRUE)
  }

  # extract information from CV
  # the mean cross-validated error, a vector of length nlam
  cvm <- rowMeans(mat_mse)
  cvsd <- apply(mat_mse, 1, stats::sd)
  # the index of best lambda
  ibest <- which.min(cvm)

  idx <- fit$idx
  coef <- fit$coef[, ibest]
  compact <- cbind(idx[which(coef != 0), , drop = FALSE], coef[coef != 0])
  colnames(compact) <- c("index_1", "index_2", "coefficient")

  # finally return the best lambda
  out <- list(n = n,
              p = p,
              a0 = fit$a0[ibest],
              compact = compact,
              fit = fit,
              fitted = fit$fitted[, ibest],
              lambda = lambda,
              cvm = cvm,
              cvsd = cvsd,
              foldid = foldid,
              ibest = ibest,
              call = match.call())
  class(out) <- "cv.sprinter"
  return(out)
}
