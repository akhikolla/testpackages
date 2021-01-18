


#' Creates a probability model for a latent ordered network model
#'
#'
#' @param formula A LOLOG formula. See \code{link{lolog}}
#' @param theta Parameter values.
#'
#'
#' @return
#' An Rcpp object representing the likelihood model
#'
#'
#' @examples
#' # See the methods of the objects returned by this function
#' UndirectedLatentOrderLikelihood
#' DirectedLatentOrderLikelihood
#'
#' # A Barabasi-Albert type graph model with 1000 vertices
#' el <- matrix(0, nrow=0, ncol=2)
#' net <- new(UndirectedNet, el, 1000L)
#' lolik <- createLatentOrderLikelihood(net ~ preferentialAttachment(), theta=1)
#' banet <- lolik$generateNetwork()$network # generate a random network from the model
#' degrees <- banet$degree(1:1000)
#' hist(degrees, breaks=100) # plot the degree distribution
#' order <- banet[["__order__"]] # The vertex inclusion order
#' 
#' # Earlier nodes have higher degrees
#' library(ggplot2)
#' qplot(order, degrees, alpha=I(.25)) + geom_smooth(method="loess")
#' 
createLatentOrderLikelihood <- function(formula, theta = NULL) {
  env <- environment(formula)
  net <- as.BinaryNet(eval(formula[[2]], envir = env))
  model <- createCppModel(formula)
  clss <- class(net)
  networkEngine <- substring(clss, 6, nchar(clss) - 3)
  LikType <-
    eval(parse(
      text = paste("lolog::", networkEngine, "LatentOrderLikelihood", sep = "")
    ))
  lik <- new(LikType, model)
  if (!is.null(theta)) {
    lik$setThetas(theta)
  }
  lik
}


.createLatentOrderLikelihoodFromTerms <- function(terms, net, theta = NULL) {
  net <- as.BinaryNet(net)
  model <- .makeCppModelFromTerms(terms, net, theta)
  clss <- class(net)
  networkEngine <- substring(clss, 6, nchar(clss) - 3)
  LikType <-
    eval(parse(
      text = paste("lolog::", networkEngine, "LatentOrderLikelihood", sep = "")
    ))
  lik <- new(LikType, model)
  if (!is.null(theta)) {
    lik$setThetas(theta)
  }
  lik
}
