#' Junction Tree
#'
#' Construction of a junction tree and message passing
#' 
#' @param x An object return from \code{compile}
#' @param evidence A named vector. The names are the variabes and the elements
#' are the evidence
#' @param flow Either "sum" or "max"
#' @param propagate Either "no", "collect" or "full".
#' @return A \code{jt} object
#' @seealso \code{\link{query_belief}}, \code{\link{mpe}},
#' \code{\link{get_cliques}}, \code{\link{get_clique_root}},
#' \code{\link{propagate}}
#' @examples
#'
#' # Setting up the network
#' # ----------------------
#'
#' library(igraph)
#' el <- matrix(c(
#' "A", "T",
#' "T", "E",
#' "S", "L",
#' "S", "B",
#' "L", "E",
#' "E", "X",
#' "E", "D",
#' "B", "D"),
#'  nc = 2,
#'  byrow = TRUE
#' )
#' 
#' g <- igraph::graph_from_edgelist(el)
#' plot(g)
#' # -----------------------
#'
#' # Data
#' # ----
#' # We use the asia data; see the man page (?asia)
#'
#' # Compilation
#' # -----------
#' cl <- cpt_list(asia, g) # Checking and conversion
#' cp <- compile(cl)
#'
#' # After the network has been compiled, the graph has been triangulated and
#' # moralized. Furthermore, all conditional probability tables (CPTs) has been
#' # designated one of the cliques (in the triangulated and moralized graph).
#'
#' # Example 1: sum-flow without evidence
#' # ------------------------------------
#' jt1 <- jt(cp)
#' plot(jt1)
#' print(jt1)
#' query_belief(jt1, c("E", "L", "T"))
#' query_belief(jt1, c("B", "D", "E"), type = "joint")
#'
#'
#' # Notice, that jt1 is equivalent to:
#' # jt1 <- jt(cp, propagate = "no")
#' # jt1 <- propagate(jt1, prop = "full")
#' # That is; it is possible to postpone the actual propagation
#' 
#' # Example 2: sum-flow with evidence
#' # ---------------------------------
#' e2  <- c(A = "y", X = "n")
#' jt2 <- jt(cp, e2) 
#' query_belief(jt2, c("B", "D", "E"), type = "joint")
#'
#' # Notice that, the configuration (D,E,B) = (y,y,n) has changed
#' # dramatically as a consequence of the evidence
#'
#' # We can get the probability of the evidence:
#' query_evidence(jt2)
#' 
#' # Example 3: max-flow without evidence
#' # ------------------------------------
#' jt3 <- jt(cp, flow = "max")
#' mpe(jt3)
#'
#' 
#' # Example 4: max-flow with evidence
#' # ---------------------------------
#' e4  <- c(T = "y", X = "y", D = "y")
#' jt4 <- jt(cp, e4, flow = "max")
#' mpe(jt4)
#' 
#' # Notice, that T, E, S, B, X and D has changed from "n" to "y"
#' # as a consequence of the new evidence e4
#'
#' 
#' # Example 5: specifying a root node and only collect to save run time
#' # -------------------------------------------------------------------------
#'
#' \donttest{
#'   cp5 <- compile(cpt_list(asia, g), "X")
#'   jt5 <- jt(cp5, propagate = "collect")
#'   query_belief(jt5, get_cliques(jt5)$C1, "joint")
#' }
#'
#' # We can only query from the root clique now (clique 1)
#' # but we have ensured that the node of interest, "X", does indeed live in
#' # this clique
#' 
#' # Example 6: Compiling from a list of conditional probabilities
#' # -------------------------------------------------------------------------
#'
#' # * We need a list with CPTs which we extract from the asia2 object
#' #    - the list must be named with child nodes
#' #    - The elements need to be array-like objects
#'
#' cl  <- cpt_list(asia2)
#' cp6 <- compile(cl, save_graph = TRUE)
#'
#' # Inspection; see if the graph correspond to the cpts
#' # g <- dag(cp6)
#' # plot(g) 
#'
#' jt6 <- jt(cp6)
#' query_belief(jt6, c("either", "smoke"))
#' 
#' @export
jt <- function(x, evidence = NULL, flow = "sum", propagate = "full") UseMethod("jt")

#' @rdname jt
#' @export
jt.charge <- function(x, evidence = NULL, flow = "sum", propagate = "full") {

  if (!is.null(evidence)) {
    if (!valid_evidence(attr(x, "dim_names"), evidence)) {
      stop("evidence is not on correct form", call. = FALSE)
    }
  }

  j <- new_jt(x, evidence, flow)
  attr(j, "propagated") <- "no"

  if (propagate == "no") {
    return(j)
  } else if (propagate == "collect") {
    m <- send_messages(j, flow)    
    while (attr(m, "direction") != "distribute") m <- send_messages(m, flow)
    attr(m, "propagated") <- "collect"
    return(m)
  } else {
    m <- send_messages(j, flow)
    while (attr(m, "direction") != "full") m <- send_messages(m, flow)
    attr(m, "propagated") <- "full"
    return(m)
  }
  stop("propagate must be either 'no', 'collect' or full", call. = TRUE)
}

#' Propagation of junction trees
#'
#' Given a junction tree object, propagation is conducted
#' 
#' @param x A junction tree object \code{jt}
#' @param prop Either "collect" or "full".
#' @seealso \code{\link{jt}}
#' @examples
#' # See Example 1 in the 'jt' function
#' @export
propagate <- function(x, prop = "full") UseMethod("propagate")

#' @rdname propagate
#' @export
propagate.jt <- function(x, prop = "full") {

  if (prop == "collect") {

    if (attr(x, "propagated") == "collect") return(x)

    if (attr(x, "propagated") == "full") {
      stop("the junction tree is already propageted fully", call. = FALSE)
    }
    
    m <- send_messages(x, attr(x, "flow"))    
    while (attr(m, "direction") != "distribute") m <- send_messages(m, attr(x, "flow"))

    attr(m, "propagated") <- "collect"
    return(m)
    
  } else if (prop == "full"){

    if (attr(x, "propagated") == "full") return(x)
    m <- send_messages(x, attr(x, "flow"))

    while (attr(m, "direction") != "full") m <- send_messages(m, attr(x, "flow"))
    attr(m, "propagated") <- "full"
    return(m)
  } else {
    
    stop("propagate must be either 'collect' or full", call. = TRUE)  
  }
}



#' Most Probable Explanation
#'
#' Returns the most probable explanation given the evidence entered in the
#' junction tree
#' 
#' @param x A junction tree object, \code{jt}, with max-flow.
#' @seealso \code{\link{jt}}
#' @examples
#' # See the 'jt' function
#' @export
mpe <- function(x) UseMethod("mpe")

#' @rdname mpe
#' @export
mpe.jt <- function(x) {
  if (attr(x, "flow") != "max") stop("The flow of the junction tree is not 'max'.")
  attr(x, "mpe")
}


#' Return the cliques of a junction tree
#'
#' @param x A junction tree object, \code{jt}.
#' @seealso \code{\link{jt}}
#' @examples
#' # See Example 5 of the 'jt' function 
#'
#' @rdname get_cliques
#' @export
get_cliques <- function(x) UseMethod("get_cliques")

#' @rdname get_cliques
#' @export
get_cliques.jt <- function(x) x$cliques

#' @rdname get_cliques
#' @export
get_clique_root <- function(x) UseMethod("get_clique_root")

#' @rdname get_cliques
#' @export
get_clique_root.jt <- function(x) attr(x, "clique_root")


#' Query Evidence 
#'
#' Get the probability of the evidence entered in the junction tree object
#'
#' @param x A junction tree object, \code{jt}.
#' @examples
#' # See the 'jt' function
#' @seealso \code{\link{jt}}, \code{\link{mpe}}
#' @export
query_evidence <- function(x) UseMethod("query_evidence")

#' @rdname query_evidence
#' @export
query_evidence.jt <- function(x) {
  if(attr(x, "flow") != "sum") {
    stop("The flow of the junction tree must be 'sum'.")
  }
  if (attr(x, "propagated") == "no") {
    stop("In order to query the probabilty of evidence, ",
      "the junction tree must at least be propagted to ",
      "the root node.")
  }
  return(attr(x, "probability_of_evidence"))
}


#' Query probabilities
#'
#' Get probabilities from a junction tree object
#'
#' @param x A junction tree object, \code{jt}.
#' @param nodes The nodes for which the probability is desired
#' @param type Either 'marginal' or 'joint'
#' @examples
#' # See the 'jt' function
#' @seealso \code{\link{jt}}, \code{\link{mpe}}
#' @export
query_belief <- function(x, nodes, type = "marginal") UseMethod("query_belief")


#' @rdname query_belief
#' @export
query_belief.jt <- function(x, nodes, type = "marginal") {
  
  if (type %ni% c("marginal", "joint")) {
    stop("Type must be 'marginal' or 'joint'.", call. = FALSE)
  }
  
  if (attr(x, "flow") == "max") {
    stop("It does not make sense to query probablities from a junction tree with max-flow. ",
      "Use 'mpe' to obtain the max configuration.", call. = FALSE)
  }

  has_rn <- has_root_node(x)

  if (has_rn) if (!all(nodes %in% names(attr(x$charge$C$C1, "dim_names")))) {
    stop(
      "All nodes must be in the root node (clique 1) ",
      "since the junction tree has only collected! ",
      "See get_cliques(x) to find the nodes in the root node.",
      call. = FALSE
    )
  }
  
  node_lst <- if (type == "joint") {
    list(nodes)
  } else {
    as.list(nodes)
  }
  
  .query <- lapply(node_lst, function(z) {
    
    if (has_rn) {
      sd <- setdiff(x$cliques$C1, z)
      if (!neq_empt_chr(sd)) return(x$charge$C$C1)
      return(sparta::marg(x$charge$C$C1, sd))
    }
    
    # TODO: Also check the separators! They may be much smaller!!!
    # - especially for type = "marginal"!

    # TODO: Use mapply over cliques and separators?

    in_which_cliques <- .map_lgl(x$cliques, function(clq) all(z %in% clq))
    
    if (!any(in_which_cliques) && type == "joint") {
      stop("The function does not, at the moment, support queries of ",
        "nodes that belong to different cliques. ",
        "Use plot(x) or get_cliques(x) to see ",
        "the cliques of the junction tree."
      )
    }

    index_in_which_cliques <- which(in_which_cliques)
    length_of_possible_cliques <- .map_int(x$cliques[in_which_cliques], length)
    idx <- index_in_which_cliques[which.min(length_of_possible_cliques)]
    pot <- x$charge$C[[idx]]
    rm_var <- setdiff(names(attr(pot, "dim_names")), z)
    
    return(sparta::marg(pot, rm_var))
  })

  if (type == "joint") {
    return(sparta::as_array(.query[[1]]))
  } else {
    return(structure(lapply(.query, function(z) sparta::as_array(z)), names = nodes))
  }
}

#' A print method for junction trees
#'
#' @param x A junction tree object, \code{jt}.
#' @param ... For S3 compatability. Not used.
#' @seealso \code{\link{jt}}
#' @export
print.jt <- function(x, ...) {
  cls <- paste0("<", paste0(class(x), collapse = ", "), ">")
  direction <- attr(x, "direction")
  flow <- attr(x, "flow")
  nv  <- ncol(x$clique_graph)
  ne  <- sum(x$clique_graph)/2  
  clique_sizes <- .map_int(x$cliques, length)
  max_C <- max(clique_sizes)
  min_C <- min(clique_sizes)
  avg_C <- mean(clique_sizes)

  cat(" Junction Tree",
    "\n -------------------------",
    "\n  Propagated:", attr(x, "propagated"),
    "\n  Flow:", flow,
    "\n  Nodes:", nv,
    "\n  Edges:", ne, "/", nv*(nv-1)/2,
    "\n  Cliques:", length(x$cliques),
    "\n   - max:", max_C,
    "\n   - min:", min_C,
    "\n   - avg:", round(avg_C, 2), 
    paste0("\n  ", cls),
    "\n -------------------------\n"
  )
  
}

#' A plot method for junction trees
#'
#' @param x A junction tree object, \code{jt}.
#' @param ... For S3 compatability. Not used.
#' @seealso \code{\link{jt}}
#' @export
plot.jt <- function(x, ...) {
  direction <- attr(x, "direction")
  y <- if (direction == "collect") {
    list(
      cliques   = x$schedule$collect$cliques,
      tree      = x$schedule$collect$tree,
      type      = "directed"
    )
  } else if (direction == "distribute") {
    list(
      cliques = x$schedule$distribute$cliques,
      tree    = x$schedule$distribute$tree,
      type    = "directed"
    )
  } else {
    list(
      cliques = x$cliques,
      tree    = x$clique_graph,
      type    = "undirected"
    )
  }

  .names <- unlist(lapply(y$cliques, function(z) paste(z, collapse = "\n")))
  dimnames(y$tree) <- list(.names, .names)
  g <- igraph::graph_from_adjacency_matrix(y$tree, y$type)
  graphics::plot(g, ...)
}
