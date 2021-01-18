#' @import methods
#' @name ThetaYList
#' @title ThetaYList-class
#' @description Definiton of ThetaYList parameter sets
#'
#' @slot tao A numeric vector
#' @slot psy A list value
#' @slot M A list value
#' @slot lambda A list value
#' @slot Y A list value
setClass(
  "ThetaYList",
  slots = c(
    tao = "vector",
    psy = "list",
    M = "list",
    lambda = "list",
    Y = "list"
  ),
  prototype = list(
    tao = c(),
    psy = list(),
    M = list(),
    lambda = list(),
    Y = list()
  )
)

setValidity("ThetaYList", function(object) {
  if (any(object@tao > 1) | any(object@tao < 0)) {
    "@tao should be in range 0 to 1"
  }
})


#' #' Getter
#' #' @import methods
#' #' @name tao
#' #' @title tao-getter
#' #' @aliases tao
#' #' @description Definiton of hyper parameter sets
#' setGeneric("tao", function(x) standardGeneric("tao"))
#' setMethod("tao", "ThetaYList", function(x) x@tao)
#' setGeneric("psy", function(x) standardGeneric("psy"))
#' setMethod("psy", "ThetaYList", function(x) x@psy)
#' setGeneric("M", function(x) standardGeneric("M"))
#' setMethod("M", "ThetaYList", function(x) x@M)
#' setGeneric("lambda", function(x) standardGeneric("lambda"))
#' setMethod("lambda", "ThetaYList", function(x) x@lambda)
#' setGeneric("Y", function(x) standardGeneric("Y"))
#' setMethod("Y", "ThetaYList", function(x) x@Y)
#'
#' #' Setter with validation
#' setGeneric("tao<-", function(x, value) standardGeneric("tao<-"))
#' setMethod("tao<-", "ThetaYList", function(x, value) {
#'   x@tao <- value
#'   validObject(x)
#'   x
#' })
#' setGeneric("psy<-", function(x, value) standardGeneric("psy<-"))
#' setMethod("psy<-", "ThetaYList", function(x, value) {
#'   x@psy <- value
#'   validObject(x)
#'   x
#' })
#' setGeneric("M<-", function(x, value) standardGeneric("M<-"))
#' setMethod("M<-", "ThetaYList", function(x, value) {
#'   x@M <- value
#'   validObject(x)
#'   x
#' })
#' setGeneric("lambda<-", function(x, value) standardGeneric("lambda<-"))
#' setMethod("lambda<-", "ThetaYList", function(x, value) {
#'   x@lambda <- value
#'   validObject(x)
#'   x
#' })
#' setGeneric("Y<-", function(x, value) standardGeneric("Y<-"))
#' setMethod("Y<-", "ThetaYList", function(x, value) {
#'   x@Y <- value
#'   validObject(x)
#'   x
#' })
