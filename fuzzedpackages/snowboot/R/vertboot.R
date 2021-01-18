#' Bootstrapping a Network with Vertex Bootstrap
#'
#' This function bootstraps the original network using a vertex bootstrap technique.
#' See \insertCite{snijders_borgatti_1999;textual}{snowboot}
#' and \insertCite{chen_etal_2018_snowboot;textual}{snowboot}.
#'
#' @param m1 An adjacency matrix of the original network.
#' @param boot_rep A positive integer number, the number of bootstrap replications.
#' @references
#' \insertAllCited{}
#' @return A list of bootstrapped networks as adjacency matrices.
#' @export
#' @examples
#' graph_ex <- igraph::graph_from_edgelist(artificial_networks[[1]]$edges)
#' m1 <- igraph::as_adjacency_matrix(graph_ex)
#' m1 <- as.matrix(m1)
#' vertboot_out <- vertboot(m1, 20)

vertboot <- function(m1, boot_rep){
  res <- list()
  for (i in 1:boot_rep) {
    blist <- sample(0:(dim(m1)[1]-1), replace = TRUE)
    res <- c(res, list(vertboot_matrix_rcpp(m1,blist)))
  }
  res
}
