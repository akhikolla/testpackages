#' @include MixtureModels.R
NULL
#' Construct a poset of gaussian mixture models.
#'
#' Creates an object representing a collection of gaussian mixture models. There
#' is one model for each fixed number of components from 1 to some specified
#' maximum. In particular each model is identified by a single number
#' specifiying the number of components in the model. Models are naturally
#' ordered by inclusion so that, for example, a model with 2 components comes
#' before a model with 3 or more components.
#'
#' @name GaussianMixtures
#' @usage GaussianMixtures(maxNumComponents = 1, dim = 1,
#'                  phi = "default", restarts = 50)
#' @export GaussianMixtures
#'
#' @param maxNumComponents the maximum number of gaussian components to
#'                         consider in a mixture.
#' @param dim the ambient dimension in which the gaussian mixtures reside.
#'        Default is 1, corresponding to gaussian mixtures on the real line.
#' @param phi parameter controlling the strength of the sBIC penalty.
#'
#' @param restarts the number of random restarts to perform when computing the
#'        MLE.
#'
#' @return An object representing the collection.
NULL
setConstructorS3("GaussianMixtures",
                 function(maxNumComponents = 1, dim = 1, phi = "default",
                          restarts = 50) {
                   numModels = maxNumComponents
                   prior = rep(1, numModels)

                   # Generate the partial order of the models
                   if (maxNumComponents == 1) {
                     E = matrix(numeric(0), ncol = 2)
                     g = igraph::graph.empty(1)
                   } else {
                     E = cbind(seq(1, numModels - 1), seq(2, numModels))
                     g = igraph::graph.edgelist(E, directed = TRUE)
                   }
                   topOrder = as.numeric(igraph::topological.sort(g))

                   dimension = rep(1, numModels)
                   for (j in topOrder) {
                     dimension[j] = choose(dim + 2, 2) * j - 1
                   }

                   if (phi == "default") {
                     phi = (dimension[1] + 1) / 2
                   }

                   extend(
                     MixtureModels(),
                     "GaussianMixtures",
                     .numModels = numModels,
                     .prior = prior,
                     .E = E,
                     .posetAsGraph = g,
                     .topOrder = topOrder,
                     .ambientDim = dim,
                     .dimension = dimension,
                     .phi = phi,
                     .restarts = restarts
                   )
                 })

#' @rdname   getTopOrder
#' @name     getTopOrder.GaussianMixtures
#' @S3method getTopOrder GaussianMixtures
#' @usage    \method{getTopOrder}{GaussianMixtures}(this)
#' @export   getTopOrder.GaussianMixtures
setMethodS3("getTopOrder", "GaussianMixtures", function(this) {
  return(this$.topOrder)
}, appendVarArgs = F)

#' @rdname   getPrior
#' @name     getPrior.GaussianMixtures
#' @S3method getPrior GaussianMixtures
#' @usage    \method{getPrior}{GaussianMixtures}(this)
#' @export   getPrior.GaussianMixtures
setMethodS3("getPrior", "GaussianMixtures", function(this) {
  return(this$.prior)
}, appendVarArgs = F)

#' @rdname   getNumModels
#' @name     getNumModels.GaussianMixtures
#' @S3method getNumModels GaussianMixtures
#' @usage    \method{getNumModels}{GaussianMixtures}(this)
#' @export   getNumModels.GaussianMixtures
setMethodS3("getNumModels", "GaussianMixtures", function(this) {
  return(this$.numModels)
}, appendVarArgs = F)

#' Set data for the gaussian mixture models.
#'
#' Sets the data to be used by the gaussian mixture models when computing MLEs.
#'
#' @name     setData.GaussianMixtures
#' @S3method setData GaussianMixtures
#' @usage    \method{setData}{GaussianMixtures}(this, data)
#' @export   setData.GaussianMixtures
#'
#' @param this the GaussianMixtures object.
#' @param data the data to be set, a matrix where each row corresponds to a single
#'        multivariate observation. If the corresponding GaussianMixtures object
#'        has ambient dimension 1, then data may be a numeric vector of
#'        observations.
NULL
setMethodS3("setData", "GaussianMixtures", function(this, data) {
  X = data
  if (is.vector(X)) {
    X = matrix(as.numeric(X), ncol = 1)
  }
  if (this$.ambientDim != ncol(X)) {
    stop(paste("Attempting to set data in GaussianMixtures whose ambient",
               "dimesion does not agree with the number of columns in the data."
               ))
  }
  this$.X = X
  this$.logLikes = rep(NA, this$getNumModels())
  this$.mles = rep(list(NA), this$getNumModels())
}, appendVarArgs = F)

#' @rdname   getData
#' @name     getData.GaussianMixtures
#' @S3method getData GaussianMixtures
#' @usage    \method{getData}{GaussianMixtures}(this)
#' @export   getData.GaussianMixtures
setMethodS3("getData", "GaussianMixtures", function(this) {
  if (is.null(this$.X)) {
    throw("Data has not yet been set")
  }
  return(this$.X)
}, appendVarArgs = F)

#' @rdname   getNumSamples
#' @name     getNumSamples.GaussianMixtures
#' @S3method getNumSamples GaussianMixtures
#' @usage    \method{getNumSamples}{GaussianMixtures}(this)
#' @export   getNumSamples.GaussianMixtures
setMethodS3("getNumSamples", "GaussianMixtures", function(this) {
  return(length(this$getData()))
}, appendVarArgs = F)

#' @rdname   logLikeMle
#' @name     logLikeMle.GaussianMixtures
#' @S3method logLikeMle GaussianMixtures
#' @usage    \method{logLikeMle}{GaussianMixtures}(this, model, ...)
#' @export   logLikeMle.GaussianMixtures
setMethodS3("logLikeMle", "GaussianMixtures", function(this, model, ...) {
  if (!is.na(this$.logLikes[model])) {
    return(this$.logLikes[model])
  }
  X = this$getData()
  N = nrow(X)

  if (ncol(X) == 1) {
    mclustModelName = "V"
  } else {
    mclustModelName = "VVV"
  }
  fit = mclust::Mclust(X, G = model, model = mclustModelName)
  logLike = fit$loglik
  mle = fit$parameters
  for (i in 1:this$.restarts) {
    my.z = matrix(rexp(N * model), N, model)
    my.z = my.z / rowSums(my.z)
    temp.fit = mclust::me(modelName = mclustModelName, data = X, z = my.z)
    if (!is.na(temp.fit$loglik)) {
      if (temp.fit$loglik > logLike) {
        logLike = temp.fit$loglik
        mle = temp.fit$parameters
      }
    }
  }
  this$.logLikes[model] = logLike
  this$.mles[[model]] = list(mixWeights = mle$pro, means = mle$mean,
                             vars = mle$vars)
  return(this$.logLikes[model])
}, appendVarArgs = F)

#' @rdname   mle
#' @name     mle.GaussianMixtures
#' @S3method mle GaussianMixtures
#' @usage    \method{mle}{GaussianMixtures}(this, model)
#' @export   mle.GaussianMixtures
setMethodS3("mle", "GaussianMixtures", function(this, model) {
  if (!is.na(this$.mle[[model]])) {
    return(this$.mle[[model]])
  }
  this$logLikeMle(model)
  return(this$.mle[[model]])
}, appendVarArgs = F)

#' @rdname   getDimension
#' @name     getDimension.GaussianMixtures
#' @S3method getDimension GaussianMixtures
#' @usage    \method{getDimension}{GaussianMixtures}(this, model)
#' @export   getDimension.GaussianMixtures
setMethodS3("getDimension", "GaussianMixtures", function(this, model) {
  return(this$.dimension[model])
}, appendVarArgs = F)
