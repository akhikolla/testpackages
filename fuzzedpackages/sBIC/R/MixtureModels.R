#' @include ModelPoset.R
NULL
#' Linear collections of mixture models.
#'
#' An abstract class representing a collection of mixture models that are
#' linearly ordered by the number of components in the mixture. This class
#' should not be instantiated, just extended.
#'
#' @name MixtureModels
#' @usage MixtureModels(phi = "default")
#' @export MixtureModels
#'
#' @param phi parameter controlling the strength of the sBIC penalty.
#'
#' @return An object representing the collection.
#'
#' @seealso \code{\link{GaussianMixtures}}, \code{\link{BinomialMixtures}}, \code{\link{LCAs}}
NULL
setConstructorS3("MixtureModels",
                 function(phi = "default") { extend(ModelPoset(), "MixtureModels", .phi = phi) },
                 abstract = T)

#' @rdname   parents
#' @name     parents.MixtureModels
#' @S3method parents MixtureModels
#' @usage    \method{parents}{MixtureModels}(this, model)
#' @export   parents.MixtureModels
setMethodS3("parents", "MixtureModels", function(this, model) {
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

#' @rdname   learnCoef
#' @name     learnCoef.MixtureModels
#' @usage    \method{learnCoef}{MixtureModels}(this, superModel, subModel)
#' @S3method learnCoef MixtureModels
#' @export   learnCoef.MixtureModels
setMethodS3("learnCoef", "MixtureModels", function(this, superModel, subModel) {
  i = superModel
  j = subModel
  r = this$getDimension(1) # Dimension of a single component
  lambda = 1 / 2 * min(j * r + j - 1 + this$getPhi() * (i - j),
                       r * i + j - 1)
  return(list(lambda = lambda, m = 1))
}, appendVarArgs = F)

#' Get the phi parameter.
#'
#' Gets the phi parameter controlling the strength of the sBIC penalty.
#'
#' @name     getPhi
#' @export   getPhi
#'
#' @param this the MixtureModels object.
#' @param phi the new phi value.
getPhi <- function(this, phi) {
    UseMethod("getPhi")
}
#' @rdname   getPhi
#' @name     getPhi.MixtureModels
#' @usage    \method{getPhi}{MixtureModels}(this, phi)
#' @S3method getPhi MixtureModels
#' @export   getPhi.MixtureModels
setMethodS3("getPhi", "MixtureModels", function(this, phi) {
  return(this$.phi)
}, appendVarArgs = F)

#' Set phi parameter.
#'
#' Set the phi parameter in a mixture model object to a different value.
#'
#' @name     setPhi
#' @export   setPhi
#'
#' @param this the MixtureModels object.
#' @param phi the new phi value.
setPhi <- function(this, phi) {
    UseMethod("setPhi")
}
#' @rdname   setPhi
#' @name     setPhi.MixtureModels
#' @S3method setPhi MixtureModels
#' @usage    \method{setPhi}{MixtureModels}(this, phi)
#' @export   setPhi.MixtureModels
setMethodS3("setPhi", "MixtureModels", function(this, phi) {
  if (!is.numeric(phi) || length(phi) != 1) {
    throw("Invalid phi value.")
  }
  this$.phi = phi
}, appendVarArgs = F)
