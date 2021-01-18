## ---------------------------------------------------------
##                EXPORTED HELPERS
## ---------------------------------------------------------
#' Subgraph
#'
#' Construct a subgraph with a given set of nodes removed
#'
#' @param x Character vector of nodes
#' @param g Adjacency list (named) or a adjacency matrix with dimnames given as the nodes
#' @return An adjacency list or adjacency matrix. 
#' @examples
#' adj1 <- list(a = c("b", "d"), b = c("a", "c", "d"), c = c("b", "d"), d = c("a", "c", "b"))
#' d <- data.frame(a = "", b = "", c ="", d = "") # Toy data so we can plot the graph
#' g <- gengraph(d, type = "gen", adj = adj1)
#' plot(g)
#' subgraph(c("c", "b"), adj1)
#' subgraph(c("b", "d"), as_adj_mat(adj1))
#' @export
subgraph <- function(x, g) {
  # x: vector of nodes to delete
  if (inherits(g, "matrix")) {
    keepers <- setdiff(dimnames(g)[[1]], x)
    g <- g[keepers, keepers]
    return(g)
  }
  else if (inherits(g, "list")) {
    l <- list(a = "a", b = "b")
    g <- g[-match(x, names(g))]
    g <- lapply(g, function(e) {
      rm_idx <- as.vector(stats::na.omit(match(x, e)))
      if (neq_empt_int(rm_idx)) return(e[-rm_idx])
      return(e)
    })
    return(g)
  }
  else {
    stop("g must either be a matrix of an adjacency list.", call. = FALSE)
  }
}

#' A test for decomposability in undirected graphs
#'
#' This function returns \code{TRUE} if the graph is decomposable and \code{FALSE} otherwise
#'
#' @param adj Adjacency list of an undirected graph
#' @return Logial describing whether or not \code{adj} is decomposable
#' @examples
#' # 4-cycle:
#' adj1 <- list(a = c("b", "d"), b = c("a", "c"), c = c("b", "d"), d = c("a", "c"))
#' is_decomposable(adj1) # FALSE
#' # Two triangles:
#' adj2 <- list(a = c("b", "d"), b = c("a", "c", "d"), c = c("b", "d"), d = c("a", "c", "b"))
#' is_decomposable(adj2) # TRUE
#' @export
is_decomposable <- function(adj) {
  m <- try(mcs(adj), silent = TRUE)
  if( inherits(m, "list") ) return(TRUE)
    else return(FALSE)
}

#' Finds the components of a graph
#'
#' @param adj Adjacency list or \code{gengraph} object
#' @return A list where the elements are the components of the graph
#' @export
components <- function(adj) {
  if (inherits(adj, "gengraph")) adj <- adj_lst(adj)
  nodes <- names(adj)
  comps <- list()
  comps[[1]] <- dfs(adj, nodes[1])
  while (TRUE) {
    new_comp  <- setdiff(nodes, unlist(comps))
    if (identical(new_comp, character(0))) return(comps)
    comps <- c(comps, list(dfs(adj[new_comp], new_comp[1])))
  }
  return(comps)
}


#' Make a complete graph
#'
#' A helper function to make an adjacency list corresponding to a complete graph
#'
#' @param nodes A character vector containing the nodes to be used in the graph
#' @return An adjacency list of a complete graph
#' @examples
#' d  <- derma[, 5:8]
#' cg <- make_complete_graph(colnames(d))
#' @export
make_complete_graph <- function(nodes) {
  structure(lapply(seq_along(nodes), function(k) {
    nodes[-which(nodes == nodes[k])]
  }), names = nodes)
}

#' Make a null graph
#'
#' A helper function to make an adjacency list corresponding to a null graph (no edges)
#'
#' @param nodes A character vector containing the nodes to be used in the graph
#' @return An adjacency list the null graph with no edges
#' @examples
#' d  <- derma[, 5:8]
#' ng <- make_null_graph(colnames(d))
#' @export
make_null_graph <- function(nodes) {
  structure(lapply(seq_along(nodes), function(x) {
    character(0)
  }), names = nodes)
}
