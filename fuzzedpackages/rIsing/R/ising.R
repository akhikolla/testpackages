#' @title High-Dimensional Ising Model Selection
#'
#' @description Ising Model selection using L1-regularized logistic regression and extended BIC.
#'
#' @param X The design matrix.
#' @param gamma (non-negative double) Parameter for the extended BIC (default 0.5). Higher gamma encourages sparsity. See references for more details.
#' @param min_sd (non-negative double) Columns of \code{X} with standard deviation less than this value will be excluded from the graph.
#' @param nlambda (positive integer) The number of parameters in the regularization path (default 50). A longer regularization path will likely yield more accurate results, but will take more time to run.
#' @param lambda.min.ratio (non-negative double) The ratio \code{min(lambda) / max(lambda)} (default \code{1e-3}).
#' @param symmetrize The method used to symmetrize the output adjacency matrix. Must be one of "min", "max", "mean" (default), or FALSE. "min" and "max" correspond to the Wainwright min/max, respectively
#' (see reference 1). "mean" corresponds to the coefficient-wise mean of the output adjacency matrix and its transpose. If FALSE, the output matrix is not symmetrized.
#'
#' @return A list containing the estimated adjacency matrix (\code{Theta}) and the optimal regularization parameter for each node (\code{lambda}), as selected by extended BIC.
#'
#' @references
#' \enumerate{
#'   \item Ravikumar, P., Wainwright, M. J. and Lafferty, J. D. (2010). High-dimensional Ising model selection using L1-regularized logistic regression. https://arxiv.org/pdf/1010.0311v1
#'   \item Barber, R.F., Drton, M. (2015). High-dimensional Ising model selection with Bayesian information criteria. https://arxiv.org/pdf/1403.3374v2
#' }
#'
#' @examples
#'
#' \dontrun{
#' # simulate a dataset using IsingSampler
#' library(IsingSampler)
#' n = 1e3
#' p = 10
#' Theta <- matrix(sample(c(-0.5,0,0.5), replace = TRUE, size = p*p), nrow = p, ncol = p)
#' Theta <- Theta + t(Theta) # adjacency matrix must be symmetric
#' diag(Theta) <- 0
#' X <- unname(as.matrix(IsingSampler(n, graph = Theta, thresholds = 0, method = "direct") ))
#' m1 <- ising(X, symmetrize = "mean", gamma = 0.5, nlambda = 50)
#'
#' # Visualize output using igraph
#' library(igraph)
#' ig <- graph_from_adjacency_matrix(m1$Theta, "undirected", weighted = TRUE, diag = FALSE)
#' plot.igraph(ig, vertex.color = "skyblue")
#' }
#' @export
ising <- function(X, gamma = 0.5, min_sd = 0, nlambda = 50, lambda.min.ratio = 1e-3, symmetrize = "mean") {
  x <- unique(as.vector(X))
  if (length(x) > 2) stop("X must contain binary data.")
  if (gamma < 0) stop("gamma must be >= 0.")
  if (min_sd < 0) stop("min_sd must be >= 0.")
  if (nlambda < 0) stop("nlambda must be > 0.")
  if (lambda.min.ratio < 0) stop("lambda.min.ratio must be >=0.")
  if (!(symmetrize %in% c("min", "max", "mean", FALSE))) stop("symmetrize must be one of 'min', 'max', 'mean', or FALSE.")
  if (!all(x %in% c(0,1))) {
    X <- apply(X, 2, function(x) c(0,1)[as.factor(x)])
  }
  X <- as.matrix(X)

  if (is.null(colnames(X)))
    colnames(X) <- seq(ncol(X))
  X_ <- X[complete.cases(X),]
  lrs <- logreg_setup(X_, rep(1,ncol(X_)), TRUE, FALSE, nlambda, as.integer(lambda.min.ratio))

  # include only columns with sd > min_sd
  randomNodes <- which(lrs$sds > min_sd)
  if (length(randomNodes) < ncol(X)) {
    print(paste("Following nodes have sd <= min_sd and are not included in graph:",
                paste(colnames(X_)[-randomNodes], collapse = ", ") ))
  }
  X_ <- X_[,randomNodes]
  Xs <- lrs$Xs[,randomNodes]
  means <- lrs$means[randomNodes]
  sds <- lrs$sds[randomNodes]
  n <- nrow(X_); p <- ncol(X_)

  Theta_hat <- matrix(nrow=p,ncol=p)
  colnames(Theta_hat) <- rownames(Theta_hat) <- colnames(X_)
  lambda_vec <- rep(NA,p)
  for (i in 1:p) {
    lambda <- regpath_ising(Xs[,-i], X_[,i], nlambda, lambda.min.ratio)
    fit <- logreg_ising(Xs[,-i], X_[,i], means[-i], sds[-i], lambda)
    wmat <- fit$wmat
    logliks <- fit$logliks
    nvars <- colSums(wmat[-1,] != 0)
    BIC <- -2 * n * logliks + nvars*(log(n) + 2*gamma*log(p))
    m_ind <- which.min(BIC)

    Theta_hat[i, -i] <- wmat[-1,m_ind]
    lambda_vec[i] = lambda[m_ind]
  }

  diag(Theta_hat) <- 0
  #symmetrize Theta_hat
  if (symmetrize != FALSE) {
    Theta_hat <- symmetrize_mat(Theta_hat, symmetrize)
  }

  return (list(Theta = Theta_hat, lambda = lambda_vec))
}

symmetrize_mat <- function(X, symmetrize) {
  tX <- t(X)
  if (symmetrize == "min") {
    return (ifelse(abs(X) < abs(tX), X, tX))
  } else if (symmetrize == "max") {
    return (ifelse(abs(X) > abs(tX), X, tX))
  } else if (symmetrize == "mean") {
    return (0.5 * (X + tX))
  } else stop("Invalid choice for symmetrize")
}

logreg_ising <- function(X, y, means, sds, lambda) {
  dt <- as.data.table(cbind(X,y))
  dt_ <- dt[, list(b = .N, y = sum(y == 1)), by = setdiff(names(dt), "y")]
  X_ <- as.matrix(dt_[,!c("b","y"),with = FALSE])
  y_ <- dt_[["y"]]
  b_ <- dt_[["b"]]
  cpp_out <- logreg_cpp(X_, y_, b_, means, sds, lambda)
  return ( list(wmat = cpp_out$wmat, logliks = cpp_out$logliks) )
}
