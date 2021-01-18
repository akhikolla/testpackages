as_undirected_igraph <- function(g) igraph::as.undirected(g)

parents_igraph <- function(g) {
  # g: dag
  # out: named list where a name is the child and elements are parents
  A  <- igraph::as_adjacency_matrix(g)
  cn <- colnames(A)
  if (is.null(cn)) {
    stop("The vertices in the igraph object must have names")
  }
  apply(A, 2, function(x) {
    names(which(x == 1L))
  })
}

parents_cpt_list <- function(x) {
  parents <- structure(lapply(seq_along(x), function(i) {
    child   <- names(x)[i]
    setdiff(names(attr(x[[i]], "dim_names")), child)
  }), names  = names(x))
}


graph_from_cpt_list <- function(x) {
  pairs <- lapply(seq_along(x), function(i) {
    child <- names(x)[i]
    parents <- setdiff(names(attr(x[[i]], "dim_names")), child)
    as.matrix(expand.grid(parents, child, stringsAsFactors = FALSE))
  })
  el <- do.call(rbind, pairs)
  igraph::graph_from_edgelist(el)
}

moralize_igraph <- function(g, parents) {
  g <- igraph::as.undirected(g)
  for (p in parents) {
    if (length(p) > 1) {
      pairs <- utils::combn(p, 2,  simplify = FALSE)
      for (ps in pairs) {
        if (!igraph::are_adjacent(g, ps[1], ps[2])) {
          g <- g + igraph::edge(ps[1], ps[2]) 
        }
      }
    }
  }
  g
}

add_joint_vars_igraph <- function(g, joint_vars) {
  if (length(joint_vars) < 1) return(g)
  pairs <- utils::combn(joint_vars, 2,  simplify = FALSE)
  for (ps in pairs) {
    if (!igraph::are_adjacent(g, ps[1], ps[2])) {
      g <- g + igraph::edge(ps[1], ps[2]) 
    }
  }
  g
}

triangulate_igraph <- function(g) {
  igraph::is.chordal(g, fillin = FALSE, newgraph = TRUE)$newgraph
}


## construct_cliques_and_parents <- function(adj, root_node = "") {
##   rip_ <- rip(adj, start_node = root_node, check = FALSE)
##   cliques <- rip_$C
##   names(cliques) <- paste("C", 1:length(cliques), sep = "")
##   return(list(cliques = cliques, parents = rip_$P))
## }

construct_cliques <- function(adj, root_node = "") {
  rip_ <- rip(adj, start_node = root_node, check = FALSE)
  structure(rip_$C, names = paste("C", 1:length(rip_$C), sep = ""))
}

construct_cliques_int <- function(adj_mat) {
  # cliques_int is needed to construct the junction tree in new_jt -> new_schedule
  dimnames(adj_mat) <- lapply(dimnames(adj_mat), function(x) 1:nrow(adj_mat))
  adj_lst_int       <- as_adj_lst(adj_mat)
  cliques_int       <- rip(adj_lst_int)$C
  lapply(cliques_int, as.integer)
}
