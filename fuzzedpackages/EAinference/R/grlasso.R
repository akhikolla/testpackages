grlassoFit <- function(X, Y, lbd, weights = rep(1, max(group)), group = 1:p,
                       Gamma, eps=1e-5, returnGamma = FALSE, initBeta = rep(1, ncol(X)))
{
  n <- nrow(X)
  p <- ncol(X)
  #--------------------
  # Error Handling
  #--------------------
  if (lbd < 0) {
    stop("lbd should be non-negative.")
  }
  if (length(group) != p) {
    stop("group must have a same length with the number of X columns")
  }
  if (length(weights) != length(unique(group))) {
    stop("weights has to have a same length as the number of groups")
  }
  if (any(weights <= 0)) {
    stop("weights should be positive.")
  }
  if (any(!1:max(group) %in% group)) {
    stop("group index has to be a consecutive integer starting from 1.")
  }
  if (!is.null(Y) & length(Y) != n) {
    stop("dimension of X and Y are not conformable.")
  }
  #--------------------
  XY <- X * c(Y)

  if (missing(Gamma)) {
    Gamma <- groupMaxEigen(X, group)
  }

  if (returnGamma) {
    return(list(coef = c(grlasso(X = X, Y = Y, XY = XY, weights = weights, group = group,
                                 lbd = lbd, Gamma = Gamma, initBeta = initBeta,
                                 eps = eps))
    , Gamma = Gamma))
  } else {
    return(list(coef = c(grlasso(X = X, Y = Y, XY = XY, weights = weights, group = group,
                                 lbd = lbd, Gamma = Gamma, initBeta = initBeta,
                                 eps = eps))))
  }
}
