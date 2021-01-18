


#' Network conversion
#' @param x The object
#' @param ... Additional parameters
as.network <- function(x, ...) {
  UseMethod("as.network")
}

#' Convert a UndirectedNet to a network object
#' @param x the object
#' @param ... unused
#' @return A network object
#' @examples
#' el <- matrix(c(1,2),ncol=2)
#'
#' #make an UndirectedNet with one edge and 5 nodes
#' net <- new(UndirectedNet, el, 5L)
#' net[1:5,1:5]
#'
#' nw <- as.network(net)
#' nw
#' @seealso \code{\link{UndirectedNet}}
#' @method as.network Rcpp_UndirectedNet
as.network.Rcpp_UndirectedNet <- function(x, ...) {
  el <- x$edges()
  attr(el, "n") <- n <- x$size()
  
  if (nrow(el) > 0)
    nw <- network(el, directed = FALSE)
  else
    nw <- network.initialize(n, directed = FALSE)
  
  for (i in which(x$nMissing(1:n) > 0)) {
    nas <- which(is.na(x[i, 1:n]))
    nw[i, nas] <- NA
  }
  
  vn <- x$variableNames(TRUE)
  if (length(vn) > 0) {
    for (i in 1:length(vn)) {
      vals <- x[[vn[i]]]
      if (vn[i] == "vertex.names") {
        network.vertex.names(nw) <- as.character(vals)
      } else if (vn[i] == "na") {
        
      } else{
        nw %v% vn[i] <-
          if (is.factor(vals))
            as.character(vals)
        else
          as.vector(vals)
      }
    }
  }
  nw
}

#' Convert a DirectedNet to a network object
#' @param x the object
#' @param ... unused
#' @return A network object
#' @examples
#' el <- matrix(c(1,2),ncol=2)
#'
#' #make an UndirectedNet with one edge and 5 nodes
#' net <- new(UndirectedNet, el, 5L)
#'
#' nw <- as.network(net)
#' nw
#' @seealso \code{\link{DirectedNet}}
#' @method as.network Rcpp_DirectedNet
as.network.Rcpp_DirectedNet <- function(x, ...) {
  el <- x$edges()
  attr(el, "n") <- n <- x$size()
  
  if (nrow(el) > 0)
    nw <- network(el, directed = TRUE)
  else
    nw <- network.initialize(n, directed = TRUE)
  
  for (i in which(x$nMissing(1:n) > 0)) {
    nas <- which(is.na(x[i, 1:n]))
    nw[i, nas] <- NA
  }
  
  vn <- x$variableNames(TRUE)
  if (length(vn) > 0) {
    for (i in 1:length(vn)) {
      vals <- x[[vn[i]]]
      if (vn[i] == "vertex.names") {
        network.vertex.names(nw) <- as.character(vals)
      } else if (vn[i] == "na") {
        
      } else{
        nw %v% vn[i] <-
          if (is.factor(vals))
            as.character(vals)
        else
          as.vector(vals)
      }
    }
  }
  nw
}

#' plot an DirectedNet object
#' @param x the Rcpp_DirectedNet object
#' @param ... additional parameters for plot.network
#' @details
#' This is a thin wrapper around \code{\link{plot.network}}.
#' @examples
#' data(ukFaculty)
#' net <- as.BinaryNet(ukFaculty)
#' plot(net, vertex.col=net[["Group"]]+1)
#' @method plot Rcpp_DirectedNet
plot.Rcpp_DirectedNet <- function(x, ...) {
  x <- as.network(x)
  plot(x, ...)
}

#' Plot an UndirectedNet object
#' @param x the object
#' @param ... additional parameters for plot.network
#' @details
#' This is a thin wrapper around \code{\link{plot.network}}.
#' @examples
#' el <- matrix(c(1,2),ncol=2)
#' net <- new(UndirectedNet, el, 5L)
#' net[1,5] <- 1
#' net[2,5] <- 1
#' plot(net)
#' @method plot Rcpp_UndirectedNet
plot.Rcpp_UndirectedNet <- function(x, ...) {
  x <- as.network(x)
  plot(x, ...)
}

#' Convert to either an UndirectedNet or DirectedNet object
#' 
#' @param x the object
#' @param ... unused
#' @return either an Rcpp_UndirectedNet or Rcpp_DirectedNet object
#' @details 
#' Converts network objects to BinaryNets. This function also converts
#' other graph formats, such as igraph and tidygraph, utilizing
#' intergraph::asNetwork.
#' @examples
#' data(ukFaculty)
#' net <- as.BinaryNet(ukFaculty)
#' net
as.BinaryNet <- function(x, ...) {
  UseMethod("as.BinaryNet")
}

#' Convert to either an UndirectedNet or DirectedNet object
#' 
#' @param x the object
#' @param ... unused
#' @return either an Rcpp_UndirectedNet or Rcpp_DirectedNet object
#' @details 
#' Converts network objects to BinaryNets. This function also converts
#' other graph formats, such as igraph and tidygraph, utilizing
#' intergraph::asNetwork.
#' @examples
#' data(ukFaculty)
#' net <- as.BinaryNet(ukFaculty)
#' net
#' @method as.BinaryNet default
as.BinaryNet.default <- function(x, ...) {
  if (inherits(x, "Rcpp_UndirectedNet"))
    return(x)
  if (inherits(x, "Rcpp_DirectedNet"))
    return(x)
  if (!inherits(x, "network")){
    x <- intergraph::asNetwork(x, ...)
  }
  
  if(is.bipartite(x))
    stop("Bipartite BinaryNets not supported")
  if(has.loops(x))
    stop("network object contains loops")
  if(!is.null(x$gal$multiple) && x$gal$multiple)
    stop("network object contains multi-edges")
  
  directed <- is.directed(x)
  el <- as.matrix(x, matrix.type = "edgelist")
  n <- attr(el, "n")
  if (directed)
    net <- new(DirectedNet, el, n)
  else
    net <- new(UndirectedNet, el, n)
  vn <- list.vertex.attributes(x)
  for (v in vn)
    net[[v]] <- x %v% v
  
  net
}

#' indexing
#' @name [
#' @aliases [,Rcpp_DirectedNet-method [,Rcpp_DirectedNet,ANY,ANY,ANY-method \S4method{[}{Rcpp_DirectedNet,ANY,ANY,ANY}
#' @param x object
#' @param i indices
#' @param j indices
#' @param ... unused
#' @param maskMissing should missing values be masked by NA
#' @param drop unused
#' @docType methods
#' @examples
#' data(ukFaculty)
#' net <- as.BinaryNet(ukFaculty)
#'
#'
#' #dyad Extraction
#' net[1:2,1:5]
#' net$outNeighbors(c(1,2,3))
#'
#' #dyad assignment
#' net[1,1:5] <- rep(NA,5)
#' net[1:2,1:5]
#' net[1:2,1:5,maskMissing=FALSE] #remove the mask over missing values and see
#' #nothing was really changed
#'
#' #node variables
#' net$variableNames()
#' net[["Group"]]
#' net[["rnorm"]] <- rnorm(net$size())
#' net[["rnorm"]]
#' @rdname extract-methods
setMethod("[", c("Rcpp_DirectedNet"),
          function(x,
                   i,
                   j,
                   ...,
                   maskMissing = TRUE,
                   drop = TRUE)
          {
            x$`[`(i, j, maskMissing)
          })

#' indexing
#' @name [
#' @aliases [,Rcpp_UndirectedNet-method [,Rcpp_UndirectedNet,ANY,ANY,ANY-method \S4method{[}{Rcpp_UndirectedNet,ANY,ANY,ANY}
#' @docType methods
#' @rdname extract-methods
setMethod("[", c("Rcpp_UndirectedNet"),
          function(x,
                   i,
                   j,
                   ...,
                   maskMissing = TRUE,
                   drop = TRUE)
          {
            x$`[`(i, j, maskMissing)
          })

#' indexing
#' @name [<-
#' @aliases [<-,Rcpp_DirectedNet-method [<-,Rcpp_DirectedNet,ANY,ANY,ANY-method \S4method{[<-}{Rcpp_DirectedNet,ANY,ANY,ANY}
#' @param value values to assign
#' @docType methods
#' @rdname extract-methods
setMethod("[<-", c("Rcpp_DirectedNet"),
          function(x, i, j, ..., value)
          {
            if (is.vector(value)) {
              if (length(value) == length(i) && length(j) == 1)
                value <- as.matrix(as.logical(value))
              else if (length(value) == length(j) && length(i) == 1)
                value <- t(as.matrix(as.logical(value)))
              else
                stop("invalid assignment")
            }
            x$`[<-`(i, j, value)
            x
          })

#' indexing
#' @name [<-
#' @aliases [<-,Rcpp_UndirectedNet-method [<-,Rcpp_UndirectedNet,ANY,ANY,ANY-method \S4method{[<-}{Rcpp_UndirectedNet,ANY,ANY,ANY}
#' @docType methods
#' @rdname extract-methods
setMethod("[<-", c("Rcpp_UndirectedNet"),
          function(x, i, j, ..., value)
          {
            if (is.vector(value)) {
              if (length(value) == length(i) && length(j) == 1)
                value <- as.matrix(as.logical(value))
              else if (length(value) == length(j) && length(i) == 1)
                value <- t(as.matrix(as.logical(value)))
              else
                stop("invalid assignment")
            }
            x$`[<-`(i, j, value)
            x
          })
