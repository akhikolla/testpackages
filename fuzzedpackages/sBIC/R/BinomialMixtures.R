#' @include MixtureModels.R
NULL
#' Construct a poset of binomial mixture models.
#'
#' Creates an object representing a collection of binomial mixture models. There
#' is one model for each fixed number of components from 1 to some specified
#' maximum. In particular each model is identified by a single number
#' specifiying the number of components in the model. Models are naturally
#' ordered by inclusion so that, for example, a model with 2 components comes
#' before a model with 3 or more components.
#'
#' @name BinomialMixtures
#' @usage BinomialMixtures(maxNumComponents = 1, phi = "default")
#' @export BinomialMixtures
#'
#' @param maxNumComponents the maximum number of components allowed in a model, will
#'                      create a hierarchy of all models with less than or equal
#'                      to this number.
#' @param phi parameter controlling the strength of the sBIC penalty.
#'
#' @return An object representing the collection.
NULL
setConstructorS3("BinomialMixtures",
                 function(maxNumComponents = 1, phi = "default") {
                   numModels =  maxNumComponents
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
                     dimension[j] = 2 * j - 1
                   }

                   if (phi == "default") {
                     phi = (dimension[1] + 1) / 2
                   }

                   extend(
                     MixtureModels(),
                     "BinomialMixtures",
                     .numModels = numModels,
                     .prior = prior,
                     .E = E,
                     .posetAsGraph = g,
                     .topOrder = topOrder,
                     .dimension = dimension,
                     .phi = phi
                   )
                 })

#' @rdname   getTopOrder
#' @name     getTopOrder.BinomialMixtures
#' @usage    \method{getTopOrder}{BinomialMixtures}(this)
#' @S3method getTopOrder BinomialMixtures
#' @export   getTopOrder.BinomialMixtures
setMethodS3("getTopOrder", "BinomialMixtures", function(this) {
  return(this$.topOrder)
}, appendVarArgs = F)

#' @rdname   getPrior
#' @name     getPrior.BinomialMixtures
#' @usage    \method{getPrior}{BinomialMixtures}(this)
#' @S3method getPrior BinomialMixtures
#' @export   getPrior.BinomialMixtures
setMethodS3("getPrior", "BinomialMixtures", function(this) {
  return(this$.prior)
}, appendVarArgs = F)

#' @rdname   getNumModels
#' @name     getNumModels.BinomialMixtures
#' @usage    \method{getNumModels}{BinomialMixtures}(this)
#' @S3method getNumModels BinomialMixtures
#' @export   getNumModels.BinomialMixtures
setMethodS3("getNumModels", "BinomialMixtures", function(this) {
  return(this$.numModels)
}, appendVarArgs = F)

#' Set data for the binomial mixture models.
#'
#' Sets the data to be used by the binomial mixture models when computing MLEs.
#'
#' @name     setData.BinomialMixtures
#' @usage    \method{setData}{BinomialMixtures}(this, data)
#' @S3method setData BinomialMixtures
#' @export   setData.BinomialMixtures
#'
#' @param this the BinomialMixtures object.
#' @param data the data to be set, should be a numeric vector of non-negative
#'        integers.
NULL
setMethodS3("setData", "BinomialMixtures", function(this, data) {

  X = data
  this$.X = X

  flexmixFit = flexmix::initFlexmix(
    X ~ 1,
    k = 1:this$getNumModels(),
    model = flexmix::FLXglm(family = "binomial"),
    control = list(minprior = 0),
    nrep = 10,
    verbose = FALSE
  )

  this$.mles = rep(list(list()), this$getNumModels())
  n = this$getNumSamples()
  for (i in 1:this$getNumModels()) {
    model = flexmix::getModel(flexmixFit, i)
    clusters = as.numeric(flexmix::clusters(model))
    params = as.numeric(flexmix::parameters(model))

    this$.mles[[i]]$binomProbs = exp(params)/(1 + exp(params))
    this$.mles[[i]]$mixWeights = as.numeric(table(factor(clusters, levels = 1:i))) / n
  }

  this$.logLikes = as.numeric(flexmix::logLik(flexmixFit))
}, appendVarArgs = F)

#' @rdname   getData
#' @name     getData.BinomialMixtures
#' @S3method getData BinomialMixtures
#' @usage    \method{getData}{BinomialMixtures}(this)
#' @export   getData.BinomialMixtures
setMethodS3("getData", "BinomialMixtures", function(this) {
  if (is.null(this$.X)) {
    throw("Data has not yet been set")
  }
  return(this$.X)
}, appendVarArgs = F)

#' @rdname   getNumSamples
#' @name     getNumSamples.BinomialMixtures
#' @S3method getNumSamples BinomialMixtures
#' @usage    \method{getNumSamples}{BinomialMixtures}(this)
#' @export   getNumSamples.BinomialMixtures
setMethodS3("getNumSamples", "BinomialMixtures", function(this) {
  return(nrow(this$getData()))
}, appendVarArgs = F)

#' @rdname   logLikeMle
#' @name     logLikeMle.BinomialMixtures
#' @S3method logLikeMle BinomialMixtures
#' @usage    \method{logLikeMle}{BinomialMixtures}(this, model, ...)
#' @export   logLikeMle.BinomialMixtures
setMethodS3("logLikeMle", "BinomialMixtures", function(this, model, ...) {
  return(this$.logLikes[model])
}, appendVarArgs = F)

#' @rdname   mle
#' @name     mle.BinomialMixtures
#' @S3method mle BinomialMixtures
#' @usage    \method{mle}{BinomialMixtures}(this, model)
#' @export   mle.BinomialMixtures
setMethodS3("mle", "BinomialMixtures", function(this, model) {
  return(this$.mle[[model]])
}, appendVarArgs = F)

#' @rdname   getDimension
#' @name     getDimension.BinomialMixtures
#' @S3method getDimension BinomialMixtures
#' @usage    \method{getDimension}{BinomialMixtures}(this, model)
#' @export   getDimension.BinomialMixtures
setMethodS3("getDimension", "BinomialMixtures", function(this, model) {
  return(this$.dimension[model])
}, appendVarArgs = F)
