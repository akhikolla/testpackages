#' Perform nfolds-cross validation
#'
#' Select tuning parameter by cross validation according to the likelihood on testing data, with and without refitting.
#'
#' @param x A n-by-p sample matrix, each row is an observation of the p-dim random vector.
#' @param w Logical. Should we use weighted version of the penalty or not? If \code{TRUE}, we use general weight. If \code{FALSE}, use unweighted penalty. Default is \code{FALSE}.
#' @param lasso Logical. Should we use l1 penalty instead of hierarchical group lasso penalty? Note that by using l1 penalty, we lose the banded structure in the resulting estimate. And when using l1 penalty, the becomes CSCS (Convex Sparse Cholesky Selection) introduced in Khare et al. (2016). Default value for \code{lasso} is \code{FALSE}.
#' @param lamlist A list of non-negative tuning parameters \code{lambda}.
#' @param nlam If lamlist is not provided, create a lamlist with length \code{nulam}. Default is 60.
#' @param flmin If lamlist is not provided, create a lamlist with ratio of the smallest and largest lambda in the list equal to \code{flmin}. Default is 0.01.
#' @param folds Folds used in cross-validation
#' @param nfolds If folds are not provided, create folds of size \code{nfolds}.
#'
#' @return A list object containing \describe{
#' \item{errs_fit: }{A \code{nlam}-by-\code{nfolds} matrix of negative Gaussian log-likelihood values on the CV test data sets. \code{errs[i,j]} is negative Gaussian log-likelihood values incurred in using \code{lamlist[i]} on fold \code{j}}.
#' \item{errs_refit: }{A \code{nlam}-by-\code{nfolds} matrix of negative Gaussian log-likelihood values of the refitting.}
#' \item{folds: }{Folds used in cross validation.}
#' \item{lamlist: }{\code{lambda} grid used in cross validation.}
#' \item{ibest_fit: }{index of \code{lamlist} minimizing CV negative Gaussian log-likelihood.}
#' \item{ibest_refit: }{index of \code{lamlist} minimizing refitting CV negative Gaussian log-likelihood.}
#' \item{i1se_fit: }{Selected value of \code{lambda} using the one-standard-error rule.}
#' \item{i1se_refit: }{Selected value of \code{lambda} of the refitting process using the one-standard-error rule.}
#' \item{L_fit: }{Estimate of L corresponding to \code{ibest_fit}.}
#' \item{L_refit: }{Refitted estimate of L corresponding to \code{ibest_refit}.}
#' }
#'
#' @examples
#' set.seed(123)
#' p <- 50
#' n <- 50
#' true <- varband_gen(p = p, block = 5)
#' x <- sample_gen(L = true, n = n)
#' res_cv <- varband_cv(x = x, w = FALSE, nlam = 40, flmin = 0.03)
#' @export
#'
#' @seealso \code{\link{varband}} \code{\link{varband_path}}
varband_cv <- function(x, w = FALSE, lasso = FALSE, lamlist = NULL, nlam = 60, flmin = 1e-2, folds = NULL, nfolds = 5) {
  n <- nrow(x)
  p <- ncol(x)

  S <- crossprod(scale(x, center=TRUE, scale=FALSE)) / n
  if(is.null(folds))
    folds <- makefolds(n, nfolds = nfolds)
  nfolds <- length(folds)

  if (is.null(lamlist)) {
    lam_max <- lammax(S = S)
    lamlist <- pathGen(nlam = nlam, lam_max = lam_max,
                       flmin = flmin, S = S)
  } else {
    nlam <- length(lamlist)
  }

  errs_fit <- matrix(NA, nlam, nfolds)
  errs_refit <- matrix(NA, nlam, nfolds)

  # error function is the negative log Gaussian likelihood
  for (i in seq(nfolds)) {
    # train on all but i-th fold:
    x_tr <- x[-folds[[i]],]
    meanx <- colMeans(x_tr)
    x_tr <- scale(x_tr, center = meanx, scale = FALSE)
    S_tr <- crossprod(x_tr) / (dim(x_tr)[1])

    path_fit <- varband_path(S = S_tr, w = w, lasso = lasso,
                             lamlist = lamlist)$path
    path_refit <- refit_path(S = S_tr, path = path_fit)

    # evaluate this on left-out fold:
    x_te <- x[folds[[i]], ]
    x_te <- scale(x_te, center = meanx, scale = FALSE)
    S_te <- crossprod(x_te) / (dim(x_te)[1])

    for (j in seq(nlam)) {
      errs_fit[j, i] <- likelihood(crossprod(path_fit[, , j]), S_te)
      errs_refit[j, i] <- likelihood(crossprod(path_refit[, , j]), S_te)
    }
  }

  m_fit <- rowMeans(errs_fit)
  se_fit <- apply(errs_fit, 1, sd) / sqrt(nfolds)
  m_refit <- rowMeans(errs_refit)
  se_refit <- apply(errs_refit, 1, sd) / sqrt(nfolds)
  ibest_fit <- which.min(m_fit)
  i1se_fit <- min(which(m_fit < m_fit[ibest_fit] + se_fit[ibest_fit]))
  ibest_refit <- which.min(m_refit)
  i1se_refit <- min(which(m_refit < m_refit[ibest_refit] + se_refit[ibest_refit]))


  fit_cv <- varband(S = S, lambda = lamlist[ibest_fit], init = path_fit[, , ibest_fit], w = w, lasso = lasso)
  refit_cv <- varband(S = S, lambda = lamlist[ibest_refit], init = path_refit[, , ibest_refit], w = w, lasso = lasso)
  refit_cv <- refit_matrix(S = S, mat = refit_cv)

  return(list(errs_fit = errs_fit, errs_refit = errs_refit,
              folds = folds, lamlist = lamlist,
              ibest_fit = ibest_fit, ibest_refit = ibest_refit,
              i1se_fit = i1se_fit, i1se_refit = i1se_refit,
              L_fit = fit_cv, L_refit = refit_cv))
}

makefolds <- function(n, nfolds) {
  # divides the indices 1:n into nfolds random folds of about the same size.
  nn <- round(n / nfolds)
  sizes <- rep(nn, nfolds)
  sizes[nfolds] <- sizes[nfolds] + n - nn * nfolds
  b <- c(0, cumsum(sizes))
  ii <- sample(n)
  folds <- list()
  for (i in seq(nfolds))
    folds[[i]] <- ii[seq(b[i] + 1, b[i + 1])]
  folds
}

likelihood <- function (Omega, S){
  # Calculate the negative log-Gaussian likelihood with
  # precision matrix Omega and sample covariance S
  return(-determinant(Omega, logarithm = TRUE)$modulus[1] + sum(S*Omega))
}
