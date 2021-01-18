elimination <- function(x, nlvls, opt = "min_sp") {
  # x: dim-named adjacency matrix
  # opt: optimization strategy; either min_sp or min_fill

  # out: set of edges that, when added to the graph, makes it triangulated

  if (opt %ni% c("min_fill", "min_sp")) {
    stop("opt must be one of 'min_fill' org 'min_sp'", call. = FALSE)
  }
  
  # triangulation set
  triang_edges <- list()
  
  # Sort accoring to x
  nlvls <- nlvls[dimnames(x)[[1]]]
  
  # non-eliminated nodes
  nodes_left = 1:ncol(x)

  while (length(nodes_left) > 1L) {
    
    # Find the optimal node to eliminate
    optimal_index <- integer(0)
    optimal_crit  <- Inf
    optimal_fam   <- integer(0)
    optimal_nei   <- integer(0)

    
    # Optimize by only consider the leaves (if any)!
    leaves <- unname(which(apply(x, 2L, sum) == 1L))
    nodes <- if (neq_empt_int(leaves)) leaves else 1:ncol(x)
    
    for (k in nodes) {
      nei_k <- which(x[, k] == 1L)
      family <- c(k, nei_k)
      
      crit <- if (opt == "min_sp") {
        prod(nlvls[family])
        # prod(.map_int(dim_names[family], length))
      } else {
        x_nei_k <- x[nei_k, nei_k]
        all_edges <- length(nei_k) * (length(nei_k) - 1L) / 2
        existing_edges <- sum(x_nei_k) / 2L
        all_edges - existing_edges
      }

      if (crit < optimal_crit) {
        optimal_crit  <- crit
        optimal_index <- k
        optimal_fam   <- family
        optimal_nei   <- nei_k
      }
      
    }

    # Complete the optimal family
    x_fam <- x[optimal_fam, optimal_fam, drop = FALSE]
    # Suffices to use x_nei!
    # x_nei <- x[optimal_nei, optimal_nei]
    max_edges <- ncol(x_fam) * (ncol(x_fam) - 1) / 2
      
    if (sum(x_fam) / 2 < max_edges) {
      x_fam_nodes <- colnames(x_fam)
      for (k in 1:ncol(x_fam)) {

        x_nei_k <- x_fam[-c(1:k), k]
        non_nei_k <- which(x_nei_k == 0L) + k # +k to compensate for -1:k!
        pairs <- lapply(seq_along(non_nei_k), function(i) {
          c(x_fam_nodes[k], x_fam_nodes[non_nei_k[i]])
        })
        
        for (p in pairs) {
          triang_edges <- push(triang_edges, p)
        }

        for (i in seq_along(pairs)) {
          p <- pairs[[i]]
          x[p[1], p[2]] <- 1L
          x[p[2], p[1]] <- 1L
        }
      }
    }
    
    x <- x[-optimal_index, -optimal_index, drop = FALSE]
    nodes_left <- nodes_left[-optimal_index]
  }

  triang_edges
}

triangulate_adjacency_matrix <- function(x, nlvls, opt = "min_sp") {
  # TODO: Test if x is already decomposable?
  edges <- elimination(x, nlvls, opt)
  for (e in edges) {
    x[e[1], e[2]] <- 1L
    x[e[2], e[1]] <- 1L
  }
  x
}
