#' @include allClasses.R
NULL

#' Prediction on a raw test set of the best logistic regression model on discretized / grouped data.
#'
#' This function discretizes a user-provided test dataset given a discretization scheme provided by an S4 \code{\link{glmdisc}} object.
#' It then applies the learnt logistic regression model and outputs its prediction (see \code{\link{predict.glm}}).
#' @name discretize
#' @rdname discretize
#' @exportMethod discretize
#' @param object glmdisc object
#' @param data the data to discretize according to the provided discretization scheme
#' @author Adrien Ehrhardt.
# #' @description This defines the generic method "discretize" which will discretize a new input dataset given a discretization scheme of S4 class \code{\link{glmdisc}}.
methods::setGeneric("discretize", function(object, data) attributes(object))
