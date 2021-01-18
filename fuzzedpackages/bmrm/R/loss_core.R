
#' Return or set gradient attribute
#' 
#' gradient attribute is used by loss/risk function to return the gradient of
#' the function at a given point together with the function value
#' 
#' @name gradient
#' @rdname gradient
#' @aliases gradient<-
#' @title Return or set gradient attribute
#' @param x any R object
#' @param value new gradient value to set
#' @param ... additional paramters 
#' @return attr(x,"gradient")
#' @export
gradient <- function(x,...) UseMethod("gradient")

#' @rdname gradient
#' @export
gradient.default <- function(x,...) attr(x, "gradient")

#' @rdname gradient
#' @export
"gradient<-" <- function(x,...,value) UseMethod("gradient<-")

#' @rdname gradient
#' @export
"gradient<-.default" <- function(x,...,value) {attr(x, "gradient") <- value;x}











#' Return or set lvalue attribute
#' 
#' lvalue attribute is used by loss/risk function to return the loss value of
#' the function at a given point together with the function gradient
#' 
#' @name lvalue
#' @rdname lvalue
#' @aliases lvalue<-
#' @title Return or set lvalue attribute
#' @param x any R object
#' @param value new loss value to set
#' @param ... additional paramters 
#' @return attr(x,"lvalue")
#' @export
lvalue <- function(x,...) UseMethod("lvalue")

#' @rdname lvalue
#' @export
lvalue.default <- function(x,...) attr(x, "lvalue")

#' @rdname lvalue
#' @export
"lvalue<-" <- function(x,...,value) UseMethod("lvalue<-")

#' @rdname lvalue
#' @export
"lvalue<-.default" <- function(x,...,value) {attr(x, "lvalue") <- value;x}







#' Return or set is.convex attribute
#' 
#' is.convex attribute is used by loss/risk function to determine if it is convex
#' 
#' @name is.convex
#' @rdname is.convex
#' @aliases is.convex<-
#' @title Return or set is.convex attribute
#' @param x any R object
#' @param value new loss value to set
#' @param ... additional paramters 
#' @return attr(x,"is.convex")
#' @export
is.convex <- function(x,...) UseMethod("is.convex")

#' @rdname is.convex
#' @export
is.convex.default <- function(x,...) {
  if (is.null(attr(x, "is.convex"))) return(FALSE)
  attr(x, "is.convex")
}

#' @rdname is.convex
#' @export
"is.convex<-" <- function(x,...,value) UseMethod("is.convex<-")

#' @rdname is.convex
#' @export
"is.convex<-.default" <- function(x,...,value) {attr(x, "is.convex") <- value;x}






# print.nrbmLoss <- function(x,...) {
#   cat(sprintf("w[%d]: [%s ...]\n",length(x),paste(pretty(head(x)),collapse=" "),"\n"))
#   cat(sprintf("attributes: %s",paste(names(attributes(x)),collapse=", "),"\n"))
#   invisible(x)
# }


