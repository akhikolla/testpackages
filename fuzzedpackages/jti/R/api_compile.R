#' Conditional probability list
#'
#' A check and conversion of cpts to be used in the junction tree algorithm
#'
#' @param x Either a named list with cpts in form of array-like object(s)
#' where names must be the child node or a \code{data.frame}
#' @param g Either a directed acyclic graph (DAG) as an igraph object or a
#' decomposable graph as an igraph object. If \code{x} is a list,
#' \code{g} must be \code{NULL}. The procdure then deduce the graph
#' from the conditional probability tables.
#' @examples
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
#' cpt_list(asia, g)
#' @export
cpt_list <- function(x, g = NULL) UseMethod("cpt_list")


#' @rdname cpt_list
#' @export
cpt_list.list <- function(x, g = NULL) { 

  if (!is_named_list(x)) {
    stop(
      "x must be a named list of cpts. ",
      "A name should be the name of the corresponding child node."
    )
  }

  if (!is.null(g)) stop("g must be 'NULL'")

  dim_names <- list()
  
  y <- lapply(seq_along(x), function(i) {
    l <- x[[i]]
    class_allowed <- any(.map_lgl(sparta::allowed_class_to_sparta(), function(x) inherits(l, x)))
    if (!class_allowed) stop("one ore more elements in x is not an array-like object")
    spar <- sparta::as_sparta(l)
    # This ensures, that the CPTs and dim_names have the same ordering of the lvls!
    dim_names <<- push(dim_names, attr(spar, "dim_names"))
    spar
  })

  names(y)  <- names(x)
  dim_names <- unlist(dim_names, FALSE)
  dim_names <- dim_names[unique(names(dim_names))]

  g <- graph_from_cpt_list(y)
  if (!igraph::is_dag(g)) stop("The cpts does not induce an acyclic graph.")

  structure(y,
    names = names(x),
    dim_names = dim_names,
    parents = parents_cpt_list(y),
    graph = g,
    class = c("cpt_list", "list")
  )
}


#' @rdname cpt_list
#' @export
cpt_list.data.frame <- function(x, g) {

  if (!igraph::is.igraph(g)) stop("g must be an igraph object")

  is_dag <- igraph::is_dag(g)
  if (!is_dag) {
    if (!igraph::is_chordal(g)$chordal) {
      stop("undirected graphs must be be decomposable")
    }
  }

  parents <- if (is_dag) {
    parents_igraph(g)
  } else {
    rip(as_adj_lst(igraph::as_adjacency_matrix(g)), check = FALSE)$P 
  }
  
  dns <- list()

  y <- lapply(seq_along(parents), function(i) {
    child <- names(parents)[i]
    pars  <- parents[[i]]
    spar  <- sparta::as_sparta(x[, c(child, pars), drop = FALSE])
    spar  <- sparta::as_cpt(spar, pars)
    # This ensures, that the CPTs and dim_names have the same ordering of the lvls!
    dns   <<- push(dns, sparta::dim_names(spar))
    spar
  })

  dns <- unlist(dns, FALSE)
  dns <- dns[unique(names(dns))]

  structure(
    y,
    names = names(parents),
    dim_names = dns,
    parents = parents,
    graph = g,
    class = c("cpt_list", "list")
  )

}

#' Compile information
#'
#' Compiled objects are used as building blocks for junction tree inference
#'
#' @param x An object returned from \code{cpt_list}
#' @param root_node A node for which we require it to live in the root
#' clique (the first clique)
#' @param joint_vars A vector of variables for which we require to be in
#' the same clique
#' @param save_graph Logical.
#' @param opt The optimization strategy used for triangulation. Either
#' 'min_fill' or 'min_sp'
#' @details The Junction Tree Algorithm performs both a forward and inward
#' message passsing (collect and distribute). However, when the forward
#' phase is finish, the root clique potential is guaranteed to be the
#' joint pmf over the variables involved in the root clique. Thus, if
#' it is known in advance that a specific variable is of interest, the
#' algortihm can be terminated after the forward phase. Use the \code{root_node}
#' to specify such a variable.
#'
#' Moreover, if interest is in some joint pmf for variables that end up
#' being in different cliques these variables must be specified in advance
#' using the \code{joint_vars} argument. The compilation step then
#' adds edges between all of these variables to ensure that at least one
#' clique contains all of them.
#' @examples
#' cptl <- cpt_list(asia2)
#' compile(cptl, joint_vars = c("bronc", "tub"))
#' @export
compile <- function(x, root_node = "", joint_vars = NULL, save_graph = FALSE, opt = "min_fill") {
  UseMethod("compile")
}

#' @rdname compile
#' @export
compile.cpt_list <- function(x, root_node = "", joint_vars = NULL, save_graph = FALSE, opt = "min_fill") {

  g       <- attr(x, "graph")
  parents <- attr(x, "parents")

  gm      <- moralize_igraph(g, parents)
  if (!is.null(joint_vars)) gm <- add_joint_vars_igraph(gm, joint_vars)

  gmt     <- triangulate_adjacency_matrix(
   igraph::as_adjacency_matrix(gm),
   .map_int(attr(x, "dim_names"), length),
   opt
  )

  adj_lst <- as_adj_lst(gmt)

  # cliques_int is needed to construct the junction tree in new_jt -> new_schedule
  dimnames(gmt) <- lapply(dimnames(gmt), function(x) 1:nrow(gmt))
  adj_lst_int       <- as_adj_lst(gmt)

  root_node_int <- ifelse(root_node != "", as.character(match(root_node, names(adj_lst))), "")
  cliques_int <- lapply(rip(adj_lst_int, root_node_int)$C, as.integer)

  cliques <- construct_cliques(adj_lst, root_node)
  charge  <- new_charge(x, cliques, parents)
  out     <- structure(
    list(charge = charge, cliques = cliques),
    root_node   = root_node,
    dim_names   = attr(x, "dim_names"),
    cliques_int = cliques_int,
    class       = c("charge", "list")
  )

  if (save_graph) attr(out, "graph") <- g
  out
}



#' DAG
#'
#' Retrieve the DAG from a compiled object
#'
#' @param x A compiled object
#' @return A dag as an \code{igraph} object

#' @rdname dag
#' @export
dag <- function(x) UseMethod("dag")

#' @rdname dag
#' @export
dag.charge <- function(x) attr(x, "graph")

#' @rdname dag
#' @export
dag.cpt_list <- function(x) attr(x, "graph")

## #' @export
## triangulate.igraph <- function(g) triangulate_igraph(g)

## #' @export
## moralize.igraph <- function(g) moralize_igraph(g, parents_igraph(g))


#' A print method for cpt lists
#'
#' @param x A \code{cpt_list} object
#' @param ... For S3 compatability. Not used.
#' @seealso \code{\link{compile}}
#' @export
print.cpt_list <- function(x, ...) {
  cls <- paste0("<", paste0(class(x), collapse = ", "), ">")
  nn  <- length(names(x))
  cat(" List of CPTs",
    "\n -------------------------\n")

  for (child in names(x)) {
    parents <- setdiff(names(sparta::dim_names(x[[child]])), child)
    if (neq_empt_chr(parents)) {
      cat("  P(", child, "|" , paste(parents, collapse = ", "), ")\n")      
    } else {
      cat("  P(", child, ")\n")
    }
  }
  
  cat(paste0("\n  ", cls),
    "\n -------------------------\n"
  )
  
}


#' A print method for compiled objects
#'
#' @param x A compiled object
#' @param ... For S3 compatability. Not used.
#' @seealso \code{\link{jt}}
#' @export
print.charge <- function(x, ...) {
  cls <- paste0("<", paste0(class(x), collapse = ", "), ">")
  nn  <- length(names(attr(x, "dim_names")))
  clique_sizes <- .map_int(x$cliques, length)
  max_C <- max(clique_sizes)
  min_C <- min(clique_sizes)
  avg_C <- mean(clique_sizes)
  cat(" Compiled network",
    "\n -------------------------",
    "\n  Nodes:", nn,
    "\n  Cliques:", length(x$cliques),
    "\n   - max:", max_C,
    "\n   - min:", min_C,
    "\n   - avg:", round(avg_C, 2), 
    paste0("\n  ", cls),
    "\n -------------------------\n"
  )
  
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#               gRain TRIANGULATION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## triangulate_grbase <- function(g, nlvls) {
##   gu    <- igraph::as_adjacency_matrix(igraph::as.undirected(g))
##   nlvls <- nlvls[dimnames(gu)[[2]]]
##   tg    <- gRbase::triang(gu, control = list(method="mcwh", nLevels = nlvls))
##   igraph::graph_from_adjacency_matrix(tg, "undirected")
## }


## compile_grbase <- function(x, root_node = "", save_graph = FALSE) {
##   g       <- attr(x, "graph")
##   parents <- attr(x, "parents")

##   nlvls   <- jti:::.map_int(attr(x, "dim_names"), length)
##   gmt     <- triangulate_grbase(jti:::moralize_igraph(g, parents), nlvls)
##   adj_mat <- igraph::as_adjacency_matrix(gmt)
##   adj_lst <- jti:::as_adj_lst(adj_mat)

##   # cliques_int is needed to construct the junction tree in new_jt -> new_schedule
##   dimnames(adj_mat) <- lapply(dimnames(adj_mat), function(x) 1:nrow(adj_mat))
##   adj_lst_int       <- jti:::as_adj_lst(adj_mat)
##   cliques_int       <- jti:::rip(adj_lst_int)$C
##   cliques_int       <- lapply(cliques_int, as.integer)

##   cliques <- jti:::construct_cliques(adj_lst, root_node)
##   charge  <- jti:::new_charge(x, cliques, parents)
##   out     <- structure(
##     list(charge = charge, cliques = cliques),
##     root_node   = root_node,
##     dim_names   = attr(x, "dim_names"),
##     cliques_int = cliques_int,
##     class       = c("charge", "list")
##   )
##   if (save_graph) attr(out, "graph") <- g
##   out
## }
