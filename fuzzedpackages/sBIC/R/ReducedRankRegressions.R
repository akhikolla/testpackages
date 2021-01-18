#' @include ModelPoset.R
NULL
#' Construct a poset of reduced rank regression models.
#'
#' Creates an object representing a collection of reduced rank regression
#' models. There is one model for each fixed rank from 1 to some specified
#' maximum. In particular each model is identified by a single number
#' specifiying the rank of the regression matrix in the model. Models are
#' naturally ordered by inclusion so that, for example, a model with 2
#' rank 2 comes before before a model with rank 3 or greater.
#'
#' @name ReducedRankRegressions
#' @usage ReducedRankRegressions(numResponses, numCovariates, maxRank)
#' @export ReducedRankRegressions
#'
#' @param numResponses the number of response variables.
#' @param numCovariates the number of covariates.
#' @param maxRank the largest rank model to be included in the collection.
#'
#' @return An object representing the collection.
NULL
setConstructorS3("ReducedRankRegressions",
                 function(numResponses = 1, numCovariates = 1, maxRank = 0) {
                   numModels = maxRank + 1
                   prior = rep(1, numModels)

                   # Generate the partial order of the models
                   if (numModels == 1) {
                     E = matrix(numeric(0), ncol = 2)
                     g = igraph::graph.empty(1)
                   } else {
                     E = cbind(seq(1, numModels - 1), seq(2, numModels))
                     g = igraph::graph.edgelist(E, directed = TRUE)
                   }
                   topOrder = as.numeric(igraph::topological.sort(g))

                   dimension = rep(NA, numModels)

                   extend(
                     ModelPoset(),
                     "ReducedRankRegressions",
                     .numModels = numModels,
                     .prior = prior,
                     .E = E,
                     .posetAsGraph = g,
                     .topOrder = topOrder,
                     .dimension = dimension,
                     .numResponses = numResponses,
                     .numCovariates = numCovariates,
                     .maxRank = maxRank
                   )
                 })

#' @rdname   getTopOrder
#' @name     getTopOrder.ReducedRankRegressions
#' @usage    \method{getTopOrder}{ReducedRankRegressions}(this)
#' @S3method getTopOrder ReducedRankRegressions
#' @export   getTopOrder.ReducedRankRegressions
setMethodS3("getTopOrder", "ReducedRankRegressions", function(this) {
  return(this$.topOrder)
}, appendVarArgs = F)

#' @rdname   getPrior
#' @name     getPrior.ReducedRankRegressions
#' @usage    \method{getPrior}{ReducedRankRegressions}(this)
#' @S3method getPrior ReducedRankRegressions
#' @export   getPrior.ReducedRankRegressions
setMethodS3("getPrior", "ReducedRankRegressions", function(this) {
  return(this$.prior)
}, appendVarArgs = F)

#' @rdname   getNumModels
#' @name     getNumModels.ReducedRankRegressions
#' @usage    \method{getNumModels}{ReducedRankRegressions}(this)
#' @S3method getNumModels ReducedRankRegressions
#' @export   getNumModels.ReducedRankRegressions
setMethodS3("getNumModels", "ReducedRankRegressions", function(this) {
  return(this$.numModels)
}, appendVarArgs = F)

#' Set data for the reduced rank regression models.
#'
#' Sets the data to be used by the reduced rank regression models when computing
#' MLEs.
#'
#' @name     setData.ReducedRankRegressions
#' @usage    \method{setData}{ReducedRankRegressions}(this, data)
#' @S3method setData ReducedRankRegressions
#' @export   setData.ReducedRankRegressions
#'
#' @param this the ReducedRankRegressions object.
#' @param data the data to be set, should be a named list with two components:
#'        \itemize{
#'          \item{X}{A matrix containing the values of covariates for each
#'                  sample. Here each COLUMN represents a single sample from
#'                  all of the covariates.}
#'          \item{Y}{A matrix containing the values of the response variables
#'                   for each sample. Again, each COLUMN is a single sample.}
#'        }
NULL
setMethodS3("setData", "ReducedRankRegressions", function(this, data) {
  X = data$X
  Y = data$Y
  if (nrow(X) != this$.numCovariates || nrow(Y) != this$.numResponses ||
      ncol(X) != ncol(Y)) {
    throw("Input data XY has incorrect dimensions.")
  }
  this$.X = X
  this$.Y = Y
  this$.logLikes = rep(NA, this$getNumModels())
  this$.unconstrainedMLE = NA
}, appendVarArgs = F)

#' @rdname   getData
#' @name     getData.ReducedRankRegressions
#' @usage    \method{getData}{ReducedRankRegressions}(this)
#' @S3method getData ReducedRankRegressions
#' @export   getData.ReducedRankRegressions
setMethodS3("getData", "ReducedRankRegressions", function(this) {
  if (is.null(this$.X)) {
    throw("Data has not yet been set")
  }
  return(list(X = this$.X, Y = this$.Y))
}, appendVarArgs = F)

#' @rdname   getNumSamples
#' @name     getNumSamples.ReducedRankRegressions
#' @usage    \method{getNumSamples}{ReducedRankRegressions}(this)
#' @S3method getNumSamples ReducedRankRegressions
#' @export   getNumSamples.ReducedRankRegressions
setMethodS3("getNumSamples", "ReducedRankRegressions", function(this) {
  return(ncol(this$getData()$X))
}, appendVarArgs = F)

#' @rdname   parents
#' @name     parents.ReducedRankRegressions
#' @usage    \method{parents}{ReducedRankRegressions}(this, model)
#' @S3method parents ReducedRankRegressions
#' @export   parents.ReducedRankRegressions
setMethodS3("parents", "ReducedRankRegressions", function(this, model) {
  if (model > this$getNumModels() ||
      model < 1 || length(model) != 1) {
    throw("Invalid input model.")
  }
  if (model == 1) {
    return(numeric(0))
  } else {
    return(model - 1)
  }
}, appendVarArgs = F)

#' Help compute the MLE.
#'
#' A private method that acts as a helper function for the logLikeMLE method
#' in the ReducedRankRegressions class.
#'
#' @name     logLikeMleHelper
#' @export   logLikeMleHelper
#'
#' @param this the ReducedRankRegressions object.
#' @param model the model number.
logLikeMleHelper <- function(this, model) {
    UseMethod("logLikeMleHelper")
}
#' @rdname   logLikeMleHelper
#' @name     logLikeMleHelper.ReducedRankRegressions
#' @usage    \method{logLikeMleHelper}{ReducedRankRegressions}(this, model)
#' @S3method logLikeMleHelper ReducedRankRegressions
#' @export   logLikeMleHelper.ReducedRankRegressions
setMethodS3("logLikeMleHelper", "ReducedRankRegressions", function(this, model) {
  if (!is.matrix(this$.unconstrainedMLE)) {
    X = this$.X
    Y = this$.Y
    this$.unconstrainedMLE = Y %*% t(X) %*% solve(X %*% t(X))
    if (!is.matrix(this$.unconstrainedMLE)) {
      throw("Unexpected error in logLikeMleHelper.")
    }
    this$.Yhat = this$.unconstrainedMLE %*% X
    this$.S = svd(this$.Yhat)
  }
}, appendVarArgs = F, private = T)

#' @rdname   logLikeMle
#' @name     logLikeMle.ReducedRankRegressions
#' @usage    \method{logLikeMle}{ReducedRankRegressions}(this, model, ...)
#' @S3method logLikeMle ReducedRankRegressions
#' @export   logLikeMle.ReducedRankRegressions
setMethodS3("logLikeMle", "ReducedRankRegressions", function(this, model, ...) {
  if (!is.na(this$.logLikes[model])) {
    return(this$.logLikes[model])
  }
  X = this$.X
  Y = this$.Y
  H = model - 1 # Rank
  this$logLikeMleHelper() # Sets up the variables .C, .Yhat, and .S if they
                          # haven't been computed yet.
  C = this$.C # unconstrained MLE
  Yhat = this$.Yhat # Y predictions under the unconstrainted model
  S = this$.S # SVD of YHat

  if (H == 0) {
    UH = matrix(0, nrow(S$u), 1)
  } else{
    UH = S$u[, 1:H]
  }

  M = this$.numCovariates
  N = this$.numResponses

  if (H < min(M, N)) {
      ell = -1 / 2 * sum((Y - Yhat) ^ 2) - 1 / 2 * sum(S$d[(H + 1):length(S$d)] ^ 2)
  } else {
    ## no singular values
    ell = -1 / 2 * sum((Y - Yhat) ^ 2)
  }
  this$.logLikes[model] = ell
  return(this$.logLikes[model])
}, appendVarArgs = F)

#' @rdname   learnCoef
#' @name     learnCoef.ReducedRankRegressions
#' @usage    \method{learnCoef}{ReducedRankRegressions}(this, superModel, subModel)
#' @S3method learnCoef ReducedRankRegressions
#' @export   learnCoef.ReducedRankRegressions
setMethodS3("learnCoef", "ReducedRankRegressions", function(this, superModel, subModel) {
  ## MxH, NxH matrix sizes
  M = this$.numCovariates
  N = this$.numResponses
  H = superModel
  r = subModel
  if (r > H) {
    return(this$learnCoef(H, H))
  }

  ## case 1
  if ((N + r <= M + H) && (M + r <= N + H) && (H + r <= M + N)) {
    if (((M + H + N + r) %% 2) == 0) {
      m = 1
      lambda = -(H + r) ^ 2 - M ^ 2 - N ^ 2 + 2 * (H + r) * (M + N) + 2 *
        M * N
      lambda = lambda / 8
    }
    else{
      m = 2
      lambda = -(H + r) ^ 2 - M ^ 2 - N ^ 2 + 2 * (H + r) * (M + N) + 2 *
        M * N + 1
      lambda = lambda / 8
    }

  }
  else{
    ## case 2
    if (M + H < N + r) {
      m = 1
      lambda = H * M - H * r + N * r
      lambda = lambda / 2
    }
    else{
      ## case 3
      if (N + H < M + r) {
        m = 1
        lambda = H * N - H * r + M * r
        lambda = lambda / 2
      }
      else{
        ## case 4
        if (M + N < H + r) {
          m = 1
          lambda = M * N / 2
        }
      }
    }
  }
  return(list(lambda = lambda, m = m))
}, appendVarArgs = F)

#' @rdname   getDimension
#' @name     getDimension.ReducedRankRegressions
#' @usage    \method{getDimension}{ReducedRankRegressions}(this, model)
#' @S3method getDimension ReducedRankRegressions
#' @export   getDimension.ReducedRankRegressions
setMethodS3("getDimension", "ReducedRankRegressions", function(this, model) {
  if (!anyNA(this$.dimension[model])) {
   return(this$.dimension[model])
  }
  for (i in model) {
    if (is.na(this$.dimension[i])) {
      this$.dimension[i] = 2 * this$learnCoef(i, i)$lambda
    }
  }
  return(this$.dimension[model])
}, appendVarArgs = F)
