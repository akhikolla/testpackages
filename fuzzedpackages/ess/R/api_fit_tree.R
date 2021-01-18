tree_weights <- function(x, df) {
  # x: gengraph object
  nodes       <- colnames(df) 
  n           <- length(nodes) 
  pairs       <- utils::combn(nodes, 2,  simplify = FALSE) 
  weights     <- structure(vector(mode = "numeric", length = n * (n - 1) / 2), names = "")
  for (j in 1:n) x$MEM[[nodes[j]]] <- entropy(df[nodes[j]])
  for (p in seq_along(pairs)) {
    edge   <- sort_(pairs[[p]])
    ed     <- entropy_difference(edge, character(0), df, x$MEM)
    weights[p]        <- ed$ent
    names(weights)[p] <- edge
  }
  return(sort(weights, decreasing = TRUE))
}


kruskal      <- function(x) UseMethod("kruskal")

kruskal.tree <- function(x) {
  n          <- length(x$G_adj)
  nodes      <- names(x$G_adj)
  node_pairs <- es_to_vs(names(x$WGT))
  number_of_nodes_total <- n
  number_of_nodes_added <- 0L
  for (e in seq_along(x$WGT)) {
    if (number_of_nodes_added == number_of_nodes_total - 1) return(x)
    node1 <- node_pairs[[e]][1]
    node2 <- node_pairs[[e]][2]
    component1 <- dfs(x$G_adj, node1)
    component2 <- dfs(x$G_adj, node2)
    if (!neq_empt_chr(intersect(component1, component2))) {
      x$G_adj[[node1]] <- c(x$G_adj[[node1]], node2)
      x$G_adj[[node2]] <- c(x$G_adj[[node2]], node1)
      x$G_A[node1, node2] <- 1L                     
      x$G_A[node2, node1] <- 1L                     
      number_of_nodes_added <- number_of_nodes_added + 1L
    }
  }
  # FIX: Update CG - the print method is not correct
  return(x)
}


## as_fwd <- function(adj, mat, lst, ...) UseMethod("as_fwd")

tree_as_fwd <- function(x, df) {
  x$CG   <- rip(x$G_adj, check = FALSE)$C
  x$e    <- new_edge()
  nC     <- length(x$CG)
  x$CG_A <- Matrix::Matrix(0L, nC, nC)
  msi    <- vector("list", 0L)
  k         <- 1L
  if (nC > 1) {
    for (i in 2:nC) {
      for (j in 1:(i-1)) {
        Ci   <- x$CG[[i]]
        Cj   <- x$CG[[j]]
        Sij  <- intersect(Ci, Cj)
        if (neq_empt_chr(Sij)) { # Note: This ONLY work for trees
          x$CG_A[i,j] = 1L
          x$CG_A[j,i] = 1L
          Ci_minus_Sij <- setdiff(Ci, Sij)
          Cj_minus_Sij <- setdiff(Cj, Sij)
          edge_ij      <- sort_(c(Ci_minus_Sij, Cj_minus_Sij))
          ed           <- entropy_difference(edge_ij, Sij, df, x$MEM)
          ent_ij       <- ed$ent
          x$MEM        <- ed$mem
          if (ent_ij >= attr(x$e, "d_qic")) {
            x$e <- new_edge(edge_ij, ent_ij, k, c(i, j))
          }
          msi[[k]] <- list(C1 = Ci, C2 = Cj, S = Sij, e = structure(ent_ij, names = edge_ij))
          k <- k + 1L
        }
      }
    }
  }
  x$MSI    <- msi
  class(x) <- setdiff(c("fwd", class(x)), "tree")
  return(x)
}

fit_tree <- function(x, df, wrap = TRUE) {
  if (wrap) return(tree_as_fwd(kruskal(x), df))
  else return(kruskal(x))
}
