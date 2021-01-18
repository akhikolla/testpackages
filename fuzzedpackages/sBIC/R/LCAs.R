#' @include ModelPoset.R
NULL
#' Construct a poset of latent class analysis models.
#'
#' Creates an object representing a collection of latent class analysis models.
#' There is one model for each fixed number of latent classes from 1 to some
#' specified maximum. In particular each model is identified by a single number
#' specifiying the number of latent classes in the model. Models are naturally
#' ordered by inclusion so that, for example, a model with 2 latent classes
#' comes before a model with 3 or more latent classes.
#'
#' @name LCAs
#' @usage LCAs(maxNumClasses = 1, numVariables = 2,
#'      numStatesForVariables = 2, phi = "default")
#' @export LCAs
#'
#' @param maxNumClasses the number of classes in the largest LCA model to
#'        considered.
#' @param numVariables the number of observed variables.
#' @param numStatesForVariables the number of states for each observed variable,
#'        at the moment these must all be equal.
#' @param phi parameter controlling the strength of the sBIC penalty.
#'
#' @return An object representing the collection.
NULL
setConstructorS3("LCAs",
                 function(maxNumClasses = 1, numVariables = 2,
                          numStatesForVariables = 2, phi = "default") {
                   numModels = maxNumClasses
                   prior = rep(1, numModels)

                   # Generate the partial order of the models
                   if (maxNumClasses == 1) {
                     E = matrix(numeric(0), ncol = 2)
                     g = igraph::graph.empty(1)
                   } else {
                     E = cbind(seq(1, numModels - 1), seq(2, numModels))
                     g = igraph::graph.edgelist(E, directed = TRUE)
                   }
                   topOrder = as.numeric(igraph::topological.sort(g))

                   dimension = rep(1, numModels)
                   for (j in 1:numModels) {
                     dimension[j] = min(
                       (j - 1) + numVariables * (numStatesForVariables - 1) * j,
                       numStatesForVariables ^ numVariables - 1
                     )
                   }

                   if (phi == "default") {
                     phi = (dimension[1] + 1) / 2
                   }

                   extend(
                     MixtureModels(),
                     "LCAs",
                     .numModels = numModels,
                     .prior = prior,
                     .E = E,
                     .posetAsGraph = g,
                     .topOrder = topOrder,
                     .dimension = dimension,
                     .maxNumClasses = maxNumClasses,
                     .numVariables = numVariables,
                     .numStatesForVariables = numStatesForVariables,
                     .phi = phi
                   )
                 })

#' @rdname   getTopOrder
#' @name     getTopOrder.LCAs
#' @S3method getTopOrder LCAs
#' @usage    \method{getTopOrder}{LCAs}(this)
#' @export   getTopOrder.LCAs
setMethodS3("getTopOrder", "LCAs", function(this) {
  return(this$.topOrder)
}, appendVarArgs = F)

#' @rdname   getPrior
#' @name     getPrior.LCAs
#' @S3method getPrior LCAs
#' @usage    \method{getPrior}{LCAs}(this)
#' @export   getPrior.LCAs
setMethodS3("getPrior", "LCAs", function(this) {
  return(this$.prior)
}, appendVarArgs = F)

#' @rdname   getNumModels
#' @name     getNumModels.LCAs
#' @S3method getNumModels LCAs
#' @usage    \method{getNumModels}{LCAs}(this)
#' @export   getNumModels.LCAs
setMethodS3("getNumModels", "LCAs", function(this) {
  return(this$.numModels)
}, appendVarArgs = F)

#' Set data for the LCA models.
#'
#' Sets the data to be used by the LCA models when computing MLEs.
#'
#' @name     setData.LCAs
#' @export   setData.LCAs
#' @S3method setData LCAs
#' @usage    \method{setData}{LCAs}(this, data)
#'
#' @param this the LCAs object.
#' @param data the data to be set, should be an integer valued matrix where each
#'        row represents a single sample from the observed variables.
NULL
setMethodS3("setData", "LCAs", function(this, data) {
  if (ncol(data) != this$.numVariables) {
    throw("Input data has incorrect number of columns.")
  }
  if (!is.data.frame(data)) {
    data = as.data.frame(data)
  }
  this$.X = data
  this$.logLikes = rep(NA, this$getNumModels())
}, appendVarArgs = F)

#' @rdname   getData
#' @name     getData.LCAs
#' @S3method getData LCAs
#' @usage    \method{getData}{LCAs}(this)
#' @export   getData.LCAs
setMethodS3("getData", "LCAs", function(this) {
  if (is.null(this$.X)) {
    throw("Data has not yet been set")
  }
  return(this$.X)
}, appendVarArgs = F)

#' @rdname   getNumSamples
#' @name     getNumSamples.LCAs
#' @S3method getNumSamples LCAs
#' @usage    \method{getNumSamples}{LCAs}(this)
#' @export   getNumSamples.LCAs
setMethodS3("getNumSamples", "LCAs", function(this) {
  return(nrow(this$getData()))
}, appendVarArgs = F)

#' @rdname   logLikeMle
#' @name     logLikeMle.LCAs
#' @S3method logLikeMle LCAs
#' @usage    \method{logLikeMle}{LCAs}(this, model, ...)
#' @export   logLikeMle.LCAs
setMethodS3("logLikeMle", "LCAs", function(this, model, ...) {
  if (!is.na(this$.logLikes[model])) {
    return(this$.logLikes[model])
  }
  X = this$getData()
  f = as.formula(paste("cbind(", paste(names(X), collapse =  ","), ")~1"))
  fit = poLCA::poLCA(f, X, nclass = model, nrep = 50, maxiter = 8000,
                     verbose = FALSE)
  this$.logLikes[model] = fit$llik
  this$.mles[[model]] = fit$probs
  return(this$.logLikes[model])
}, appendVarArgs = F)

#' @rdname   mle
#' @name     mle.LCAs
#' @S3method mle LCAs
#' @usage    \method{mle}{LCAs}(this, model)
#' @export   mle.LCAs
setMethodS3("mle", "LCAs", function(this, model) {
  if (!is.na(this$.mle[[model]])) {
    return(this$.mle[[model]])
  }
  this$logLikeMle(model)
  return(this$.mle[[model]])
}, appendVarArgs = F)

#' @rdname   getDimension
#' @name     getDimension.LCAs
#' @S3method getDimension LCAs
#' @usage    \method{getDimension}{LCAs}(this, model)
#' @export   getDimension.LCAs
setMethodS3("getDimension", "LCAs", function(this, model) {
  return(this$.dimension[model])
}, appendVarArgs = F)
