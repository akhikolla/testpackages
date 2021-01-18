#' @title Hierarchy Structure
#' 
#' @description Stores variable indexes of clustering hierarchies in a fast accessible 
#' manner.
#' 
#' @param x A S3 object e.g. from hclust or dendrogram.
#' @param max.height Is the maximal height below the height of the global node which is 
#' considered.
#' @param height A vector of heights at which nodes are grouped.
#' @param names Variable names in the order in which the indexes should be assigned to the 
#' variables.
#' @param ... Further arguments.
#' 
#' @details For the HIT algorithm it is important to have the hierarchical clustering 
#' structure in a fast accessible format. This is provided by the hierarchy object 
#' generated with this function.
#' 
#' @examples 
#' ##
#' set.seed(123)
#' n <- 80
#' p <- 90
#' # x with correlated columns
#' corMat <- toeplitz((p:1/p)^3)
#' corMatQ <- chol(corMat)
#' x <- matrix(rnorm(n * p), nrow = n) %*% corMatQ
#' colnames(x) <- paste0("x", 1:p)
#' # hierarchy
#' hc <- hclust(dist(t(x)))
#' hier <- as.hierarchy(hc)
#' 
#' @export
as.hierarchy <- function(x, max.height, height, names, ...)
  UseMethod("as.hierarchy")


# #' @export
# as.hierarchy.hierarchy <- function(x, max.height, height, names, ...) x


#' @importFrom stats as.dendrogram
#' @export
as.hierarchy.hclust <- function(x, max.height, height, names, ...) {
  as.hierarchy(as.dendrogram(x), max.height, height, names, ...)
}


#' @export
as.hierarchy.dendrogram <- function(x, max.height, height, names, ...) {
  if (missing(height)) {
    height <- heightDendrogram(x)
    if (missing(max.height))
      max.height <- attr(x, "height")
  } else {
    height <- sort(height, decreasing = TRUE)
    max.height <- min(height[1L], attr(x, "height"))
  }
  height <- height[height <= max.height]
  if (attr(x, "height") > max.height)
    height <- c(attr(x, "height"), height)
  if (missing(names)) {
    names <- labels(x)
  } else if (length(setdiff(labels(x), names))) {
    stop("'x' includs variabels not in 'names'")
  }
  out <- unname(dend2hier(x, as.numeric(height), as.character(names)))
  out[[1L]][] <- sort(out[[1L]])
  names(out[[1L]]) <- names[out[[1L]]]
  class(out) <- "hierarchy"
  out
}


#' @title Heights of Dendrogram
#' 
#' @description All heights from a dendrogram. 
#' 
#' @param x A \code{\link[stats]{dendrogram}}.
#' 
#' @keywords internal
heightDendrogram <- function(x) {
  node.height <- function(d) {
    if (is.list(d)) {
      r <- attributes(d)$height
      return(c(r, node.height(d[[1L]]), node.height(d[[2L]])))
    }
    attributes(d)$height
  }
  if (!inherits(x, "dendrogram")) 
    stop("'x' is not a dendrogram")
  sort(unique(round(node.height(x), 8)), decreasing = TRUE)
}


#' @title Names of Hierarchy
#' 
#' @description Names of variables of an hierarchy.
#' 
#' @param x A \code{\link{as.hierarchy}}.
#' 
#' @export
names.hierarchy <- function(x) {
  names(x[[1L]])
}


#' @title Reorder Hierarchy
#' 
#' @description Reorder indexes according to a vector of names.
#' 
#' @param x A \code{\link{as.hierarchy}}.
#' @param names Variable names in the order in which the indexes should be assigned to the 
#' variables.
#' @param ... Further arguments passed to or from other methods (not used).
#' 
#' @importFrom stats reorder
#' @export
reorder.hierarchy <- function(x, names, ...) {
  if (!inherits(x, "hierarchy")) 
    stop("'x' is not a hierarchy")
  if (!all(names(x[[1L]]) %in% names))
    stop("'x' includs variabels not in 'names'")
  newOrder <- match(names(x[[1L]]), names)
  out <- lapply(x, function(xi, x1, newOrder) {
    xi[] <- sort(newOrder[match(xi, x1)])
    xi
  }, x1 = x[[1L]], newOrder = newOrder)
  names(out[[1L]]) <- names[out[[1L]]]
  class(out) <- "hierarchy"
  out
}
