#' Compute the sBIC.
#'
#' Computes the sBIC for a given collection of models.
#'
#' @export
#'
#' @param X the data for which the maximum likelihood estimates will be computed
#'          for the given collection of models. To see how this data should be
#'          formatted check the documentation for \code{setData.YourModelName} (e.g.
#'          \code{setData.LCAs}). If X is NULL then it is assumed that data for the
#'          model has already been set, this can be useful if you want to
#'          compute the sBIC with the same model and data several times (perhaps
#'          after changing some parameter of the model) without having to reset
#'          the data and thus recompute maximum log-likelihoods.
#'
#' @param mod an object representing a poset of models of the same type, e.g.
#'            a collection of binomial mixture models. The currently implemented
#'            models include:
#'            \itemize{
#'              \item Binomial mixtures, see \code{\link{BinomialMixtures}}.
#'              \item Factor analysis, see \code{\link{FactorAnalyses}}.
#'              \item Latent class analysis, see \code{\link{LCAs}}.
#'              \item Latent gaussian forests, see \code{\link{LatentForests}}.
#'              \item Reduced rank regression, see \code{\link{ReducedRankRegressions}}.
#'              \item 1-dimensional gaussian mixtures, see \code{\link{GaussianMixtures}}.
#'            }
#'
#'@return A named list containing the components
#'         \itemize{
#'          \item{logLike - the computed MLE log-likelihoods for each model.}
#'          \item{sBIC - the sBIC score for each model.}
#'          \item{BIC - the usual BIC score for each model.}
#'          \item{modelPoset - the input model poset mod.}
#'         }
sBIC = function(X, mod) {
  numModels = mod$getNumModels()
  topOrder = mod$getTopOrder()
  if (!is.null(X)) {
    mod$setData(X)
  } else {
    X = mod$getData()
  }
  n = mod$getNumSamples()

  # go through the vertices of g in the topological order and
  # find the list of reachable nodes (the set {j:j<i})
  reach = vector(mode = "list", length = numModels)
  for (i in topOrder) {
    parents = mod$parents(i)
    for (j in parents) {
      if (length(j) != 0) {
        reach[[i]] <- union(reach[[i]], union(reach[[j]], j))
      }
    }
  }
  if (is.null(reach[[topOrder[1]]]) == FALSE) {
    print(
      'Warning: There is no minimal element according to the toplogical order, check the construction of graph'
    )
  }

  L <- rep(0, numModels)
  p <- mod$getPrior()

  # compute p and l
  # How to represent Lij (or LogL)? Lij is a list of vectors.
  # Each list element Lij[[i]] corressponds to the larger model i
  # For each i, Lij[[i]] is a vector of Lij's for the submodels given by the "reach" list
  # The log L'_{ii} are stored in logLii
  logL <-  vector(mode = "list", length = numModels)
  logLii <- rep(0, numModels)
  Lij <-  vector(mode = "list", length = numModels)
  logLike <- rep(0, numModels)

  for (i in topOrder) {
    logLike[i] <- mod$logLikeMle(i)
    lf <- mod$learnCoef(i, i)
    logLii[i] <- logLike[i] - lf$lambda * log(n)
    for (j in reach[[i]]) {
      lf <- mod$learnCoef(i, j)
      logL[[i]][j] <-
        logLike[i] - lf$lambda * log(n) + (lf$m - 1) * log(log(n))
    }
  }
  mn = max(c(unlist(logL), logLii), na.rm = TRUE) # note na.rm=TRUE

  # Normalization
  Lii = exp(logLii - mn)
  for (i in topOrder) {
    Lij[[i]] = exp(logL[[i]] - mn)
  }

  # Compute the values of L iteratively for each model using the topological ordering
  L = rep(0, numModels)
  for (i in topOrder) {
    # the first node is the minimal node.
    # There can be more than one minimal nodes!
    if (is.null(reach[[i]]) == TRUE) {
      L[i] <- Lii[i] # Lij[[i]][1]
    } else{
      a <- p[i]
      b <- -Lii[i] * p[i] + sum(L[reach[[i]]] * p[reach[[i]]])
      c <- -sum(Lij[[i]][reach[[i]]] * L[reach[[i]]] * p[reach[[i]]])
      L[i] <- 1 / (2 * a) * (-b + sqrt(b ^ 2 - 4 * a * c))
    }
  }

  if (any(L < 0)) {
    warning(paste("Negative probabilities found, likely due to rounding",
            "errors, rounding these values to 0."))
    L = pmax(L, 0)
  }

  results = list()
  results$logLike = logLike
  results$sBIC = log(L) + mn
  results$BIC = logLike - (mod$getDimension(1:numModels) / 2) * log(n)
  results$modelPoset = mod
  return(results)
}
