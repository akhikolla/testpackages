setConstructorS3("ModelPoset",
                 function() { extend(Object(), "Object") },
                 abstract = T)

#' Topological ordering of models.
#'
#' Returns a topological ordering of models in the collection.
#'
#' @name     getTopOrder
#' @export   getTopOrder
#'
#' @param this the model poset object.
getTopOrder <- function(this) {
    UseMethod("getTopOrder")
}

#' The prior on the models.
#'
#' Returns the unnormalized prior on the collection.
#'
#' @name     getPrior
#' @export   getPrior
#'
#' @param this the model poset object.
getPrior <- function(this) {
    UseMethod("getPrior")
}

#' Number of models.
#'
#' Returns the number of models in the collection.
#'
#' @name     getNumModels
#' @export   getNumModels
#'
#' @param this the model poset object.
getNumModels <- function(this) {
    UseMethod("getNumModels")
}

#' Set data for a model poset.
#'
#' Sets the data to be used by a poset of models when computing MLEs.
#'
#' @name     setData
#' @export   setData
#'
#' @param this the model poset object.
#' @param data the data to be set.
setData <- function(this, data) {
    UseMethod("setData")
}

#' Return the set data.
#'
#' If data has been set for the object using the setData() function
#' then will return that data, otherwise will throw an error.
#'
#' @name     getData
#' @export   getData
#'
#' @param this the object from which to get the data.
getData <- function(this) {
    UseMethod("getData")
}

#' Number of samples in the set data.
#'
#' If data has been set using the setData method then returns the
#' number of samples in the data. Otherwise throws an error.
#'
#' @name     getNumSamples
#' @export   getNumSamples
#'
#' @param this the object from which to get the number of samples.
getNumSamples <- function(this) {
    UseMethod("getNumSamples")
}

#' Parents of a model.
#'
#' Returns the immediate parents of a given model, i.e. those models
#' M that are (in the poset ordering) less than the given model but for
#' which there exists no other model M' such that M < M' < (given model).
#'
#' @name     parents
#' @export   parents
#'
#' @param this the object representing the model poset.
#' @param model the model for which the parents should be found.
parents <- function(this, model) {
    UseMethod("parents")
}

#' Maximum likelihood for data.
#'
#' Computes the maximum likelihood of a model in the model poset for the
#' data set using the setData command.
#'
#' @name     logLikeMle
#' @export   logLikeMle
#'
#' @param this the object representing the model poset.
#' @param model the model for which the maximum likelihood should be computed.
#' @param ... further parameters to be passed to methods
logLikeMle <- function(this, model, ...) {
    UseMethod("logLikeMle")
}

#' Maximum likelihood estimator.
#'
#' Computes the maximum likelihood estimator of the model parameters (for a
#' given model in the collection) given the data set with setData.
#'
#' @name     mle
#' @export   mle
#'
#' @param this the object representing the model poset.
#' @param model the model for which the maximum likelihood should be computed.
mle <- function(this, model) {
    UseMethod("mle")
}


#' Learning coefficient
#'
#' Computes the learning coefficient for a model with respect to one of the
#' model's submodels.
#'
#' @name     learnCoef
#' @export   learnCoef
#'
#' @param this the object representing the model poset.
#' @param superModel the larger model of the two input models.
#' @param subModel the submodel of the larger model.
learnCoef <- function(this, superModel, subModel) {
    UseMethod("learnCoef")
}

#' Model dimension.
#'
#' Computes the dimension of a model in the model poset.
#'
#' @name     getDimension
#' @export   getDimension
#'
#' @param this the object representing the model poset.
#' @param model the model for which the dimension should be computed.
getDimension <- function(this, model) {
    UseMethod("getDimension")
}
