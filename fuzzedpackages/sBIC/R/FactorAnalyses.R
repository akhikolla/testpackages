#' @include ModelPoset.R
NULL
#' Construct a poset of factor analysis models.
#'
#' Creates an object representing a collection of factor analysis models. There
#' is one model for each fixed number of factors from 0 to some specified
#' maximum number of factors. Each model is identified by a single
#' number, in particular the model with n factors is associated with the model
#' numbered (n + 1), this is because we allow a model with 0 factors and would
#' like all models to be numbered with strictly positive integers for
#' consistency. Models are naturally ordered by inclusion so that, for example,
#' a model with 2 factors comes before a model with 3 or more factors.
#'
#' @name FactorAnalyses
#' @usage FactorAnalyses(numCovariates = 1, maxNumFactors = 0)
#' @export FactorAnalyses
#'
#' @param numCovariates the number of covariates in all of the models.
#' @param maxNumFactors the maximum number of factors allowed in a model, will
#'                      create a hierarchy of all models with less than or equal
#'                      to this number.
#'
#' @return An object representing the collection.
NULL
setConstructorS3("FactorAnalyses",
                 function(numCovariates = 1, maxNumFactors = 0) {
                   numModels = maxNumFactors + 1
                   prior = rep(1, numModels)

                   # Generate the partial order of the models
                   if (maxNumFactors == 0) {
                     E = matrix(numeric(0), ncol = 2)
                     g = igraph::graph.empty(1)
                   } else {
                     E = cbind(seq(1, numModels - 1), seq(2, numModels))
                     g = igraph::graph.edgelist(E, directed = TRUE)
                   }
                   topOrder = as.numeric(igraph::topological.sort(g))

                   dimension = rep(1, numModels)
                   for (j in topOrder) {
                     k = j - 1
                     dimension[j] = (k + 1) * numCovariates - choose(k, 2)
                   }

                   extend(
                     ModelPoset(),
                     "FactorAnalyses",
                     .numModels = numModels,
                     .prior = prior,
                     .E = E,
                     .posetAsGraph = g,
                     .topOrder = topOrder,
                     .dimension = dimension,
                     .numCovariates = numCovariates,
                     .maxNumFactors = maxNumFactors
                   )
                 })

#' @rdname   getTopOrder
#' @name     getTopOrder.FactorAnalyses
#' @usage    \method{getTopOrder}{FactorAnalyses}(this)
#' @S3method getTopOrder FactorAnalyses
#' @export   getTopOrder.FactorAnalyses
setMethodS3("getTopOrder", "FactorAnalyses", function(this) {
  return(this$.topOrder)
}, appendVarArgs = F)

#' @rdname   getPrior
#' @name     getPrior.FactorAnalyses
#' @usage    \method{getPrior}{FactorAnalyses}(this)
#' @S3method getPrior FactorAnalyses
#' @export   getPrior.FactorAnalyses
setMethodS3("getPrior", "FactorAnalyses", function(this) {
  return(this$.prior)
}, appendVarArgs = F)

#' @rdname   getNumModels
#' @name     getNumModels.FactorAnalyses
#' @usage    \method{getNumModels}{FactorAnalyses}(this)
#' @S3method getNumModels FactorAnalyses
#' @export   getNumModels.FactorAnalyses
setMethodS3("getNumModels", "FactorAnalyses", function(this) {
  return(this$.numModels)
}, appendVarArgs = F)

#' Number of factors for a model.
#'
#' Given a model number returns the number of factors in that model
#'
#' @name     getNumFactorsForModel
#' @export   getNumFactorsForModel
#'
#' @param this the FactorAnalyses object.
#' @param model the model number.
getNumFactorsForModel <- function(this, model) {
    UseMethod("getNumFactorsForModel")
}
#' @rdname   getNumFactorsForModel
#' @name     getNumFactorsForModel.FactorAnalyses
#' @S3method getNumFactorsForModel FactorAnalyses
#' @usage    \method{getNumFactorsForModel}{FactorAnalyses}(this, model)
#' @export   getNumFactorsForModel.FactorAnalyses
setMethodS3("getNumFactorsForModel", "FactorAnalyses", function(this, model) {
  if (model < 1 || model > this$getNumModels()) {
    throw("Invalid model number.")
  }
  return(model - 1)
}, appendVarArgs = F)

#' Set data for the factor analysis models.
#'
#' Sets the data to be used by the factor analysis models when computing MLEs.
#'
#' @name     setData.FactorAnalyses
#' @S3method setData FactorAnalyses
#' @usage    \method{setData}{FactorAnalyses}(this, data)
#' @export   setData.FactorAnalyses
#'
#' @param this the FactorAnalyses object.
#' @param data the data to be set, should matrix of observed responses.
#'
NULL
setMethodS3("setData", "FactorAnalyses", function(this, data) {
  if (ncol(data) != this$.numCovariates) {
    throw("Number of covariates in model does not match input matrix.")
  }
  this$.X = data
  this$.logLikes = rep(NA, this$getNumModels())
  this$.mles = rep(list(NA), this$getNumModels())
}, appendVarArgs = F)

#' @rdname   getData
#' @name     getData.FactorAnalyses
#' @S3method getData FactorAnalyses
#' @usage    \method{getData}{FactorAnalyses}(this)
#' @export   getData.FactorAnalyses
setMethodS3("getData", "FactorAnalyses", function(this) {
  if (is.null(this$.X)) {
    throw("Data has not yet been set")
  }
  return(this$.X)
}, appendVarArgs = F)

#' @rdname   getNumSamples
#' @name     getNumSamples.FactorAnalyses
#' @S3method getNumSamples FactorAnalyses
#' @usage    \method{getNumSamples}{FactorAnalyses}(this)
#' @export   getNumSamples.FactorAnalyses
setMethodS3("getNumSamples", "FactorAnalyses", function(this) {
  return(nrow(this$getData()))
}, appendVarArgs = F)

#' @rdname   parents
#' @name     parents.FactorAnalyses
#' @S3method parents FactorAnalyses
#' @usage    \method{parents}{FactorAnalyses}(this, model)
#' @export   parents.FactorAnalyses
setMethodS3("parents", "FactorAnalyses", function(this, model) {
  if (model > this$getNumModels() ||
      model <= 0 || length(model) != 1) {
    throw("Invalid input model.")
  }
  if (model == 1) {
    return(numeric(0))
  } else {
    return(model - 1)
  }
}, appendVarArgs = F)

#' @rdname   logLikeMle
#' @name     logLikeMle.FactorAnalyses
#' @S3method logLikeMle FactorAnalyses
#' @usage    \method{logLikeMle}{FactorAnalyses}(this, model, starts = 1, ...)
#' @export   logLikeMle.FactorAnalyses
#' @param    starts The number of starting values to be tried
setMethodS3("logLikeMle", "FactorAnalyses", function(this, model, starts = 1, ...) {
  if (!is.na(this$.logLikes[model])) {
    return(this$.logLikes[model])
  }
  X = this$getData()
  n = nrow(X)
  m = this$.numCovariates
  S = cov(X)
  k = this$getNumFactorsForModel(model)
  if (k == 0) {
    fit = list()
    R = diag(1, m)
    ell.hat = sum(log(diag(S))) + m
  } else{
    fit = factanal(X, k, control = list(nstart = starts))
    R = tcrossprod(fit$loadings) + diag(fit$uniquenesses)
    ell.hat = sum(log(diag(S))) + log(det(R)) + m
  }
  this$.mles[[model]] = fit
  this$.logLikes[model] = -n / 2 * ell.hat
  return(this$.logLikes[model])
}, appendVarArgs = F)

#' @rdname   mle
#' @name     mle.FactorAnalyses
#' @S3method mle FactorAnalyses
#' @usage    \method{mle}{FactorAnalyses}(this, model)
#' @export   mle.FactorAnalyses
setMethodS3("mle", "FactorAnalyses", function(this, model) {
  if (!is.na(this$.mle[[model]])) {
    return(this$.mle[[model]])
  }
  this$logLikeMle(model)
  return(this$.mle[[model]])
}, appendVarArgs = F)

#' @rdname   learnCoef
#' @name     learnCoef.FactorAnalyses
#' @S3method learnCoef FactorAnalyses
#' @usage    \method{learnCoef}{FactorAnalyses}(this, superModel, subModel)
#' @export   learnCoef.FactorAnalyses
setMethodS3("learnCoef", "FactorAnalyses", function(this, superModel, subModel) {
  m = this$.numCovariates
  k = this$getNumFactorsForModel(superModel)
  l = this$getNumFactorsForModel(subModel)
  return(list(lambda = 1 / 4 * ((k + 2) * m + l * (m - k + 1)), m = 1))
}, appendVarArgs = F)

#' @rdname   getDimension
#' @name     getDimension.FactorAnalyses
#' @S3method getDimension FactorAnalyses
#' @usage    \method{getDimension}{FactorAnalyses}(this, model)
#' @export   getDimension.FactorAnalyses
setMethodS3("getDimension", "FactorAnalyses", function(this, model) {
  return(this$.dimension[model])
}, appendVarArgs = F)
