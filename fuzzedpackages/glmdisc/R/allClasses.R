#' Class glmdisc
#'
#' Class \code{glmdisc} represents a discretization scheme associated with its optimal logistic regression model.
#'
#' @slot parameters The parameters associated with the method.
#' @slot best.disc The best discretization scheme found by the method given its parameters.
#' @slot performance The performance obtained with the method given its parameters.
#' @slot disc.data The discretized data: test set if test is TRUE; if test is FALSE and validation is TRUE, then it provides the discretized validation set. Otherwise, it provides the discretized training set.
#' @slot disc.data The continuous data: test set if test is TRUE; if test is FALSE and validation is TRUE, then it provides the discretized validation set. Otherwise, it provides the discretized training set.
#' @name glmdisc-class
#' @rdname glmdisc-class

methods::setClass("glmdisc", representation(parameters = "list", best.disc = "list", performance = "list", disc.data = "data.frame", cont.data = "data.frame"))
