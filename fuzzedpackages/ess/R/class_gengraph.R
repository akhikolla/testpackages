# Give more arguments so we can call these constructors in the internal functions
new_gengraph <- function(df, adj, cg = NULL, ...) {
  if (!setequal(colnames(df), names(adj))) stop("column names of df does not correspond to adj")
  structure(list(
    G_adj = adj,                                            # Graph as adjacency list
    G_A   = as_adj_mat(adj),                                # Graph as adjacency matrix
    CG    = cg,                                             # Clique list
    LV    = .map_int(df, function(x) length(unique(x))),  # Level vector (for stopping criteria)
    MEM   = new.env(hash = TRUE)),                          # Memoiser - saving entropies to reuse
    class = c("gengraph", "list")
  )
}

new_bwd <- function(df, adj = NULL, q = 0.5) {
  if (is.null(adj)) adj <- make_complete_graph(colnames(df))
  g    <- new_gengraph(df, adj)
  g$CG <- rip(adj, check = FALSE)$C
  g$e  <- NULL # The newly deleted edge
  structure(g, class = c("bwd", class(g)))
}

new_fwd <- function(df, adj = NULL, q = 0.5) {
  is_graph_null <- is.null(adj)
  if (is_graph_null) adj <- make_null_graph(colnames(df))
  g <- new_gengraph(df, adj)
  if (is_graph_null) {
    g$CG_A <- as_adj_mat(make_complete_graph(colnames(df))) # Can be more efficient!
    g$CG   <- as.list(names(adj))
  }
  else {
    stop("Not implementet for general graphs yet.") # Fix this to handle the conversion
  } 
  g$MSI <- list(S = NULL, max = list(e = character(0), idx = numeric(0), ins = vector("numeric", 2L)))
  g$e   <- new_edge()
  g     <- fwd_init(g, df , q)
  structure(g, class = c("fwd", class(g)))
}

new_tree <- function(df) {
  adj    <- make_null_graph(colnames(df))
  g      <- new_gengraph(df, adj)
  g$CG   <- as_adj_lst(adj)
  g$WGT  <- tree_weights(g, df) # Weights to use in kruskal procedure
  structure(g, class = c("tree", class(g)))
}

new_tfwd <- function(df) {
  g <- fit_tree(new_tree(df), df, wrap = TRUE)
  structure(g, class = setdiff(c("tfwd", class(g)), "tree"))
}

new_edge <- function(e = character(0), d_qic = 0, idx = integer(0), ins = vector("integer", 2L)) {
  # e     : edge to be deletede or added
  # d_aic : entropy difference in the two competing models
  # idx   : in fwd procedure this is the index in MSI where e lives
  # ins   : in fwd procedure this is the indicies in CG where a new clique must be inserted
  structure(e, d_qic = d_qic, idx = idx, ins = ins)
}

#' A generic and extendable structure for decomposable graphical models
#' @description A generic structure for decomposable graphical models
#' @param df data.frame
#' @param type character ("fwd", "bwd", "tree", "tfwd", "gen")
#' @param adj A user-specified adjacency list
#' @param q Penalty term in the stopping criterion (\code{0} = AIC and \code{1} = BIC)
#' @param ... Not used (for extendibility)
#' @return A \code{gengraph} object with child class \code{type} used for model selection.
#' @examples
#'
#' gengraph(derma, type = "fwd")
#' gengraph(derma, type = "bwd")
#' 
#' @seealso \code{\link{adj_lst.gengraph}}, \code{\link{adj_mat.gengraph}}, \code{\link{fit_graph}}, \code{\link{walk.fwd}}, \code{\link{walk.bwd}}
#' @export
gengraph <- function(df, type = "gen", adj = NULL, q = 0.5, ...) {
  switch(type,
    "fwd"  = new_fwd(df, adj, q),
    "bwd"  = new_bwd(df, adj, q),
    "tree" = new_tree(df),
    "tfwd" = new_tfwd(df),
    "gen"  = new_gengraph(df, adj, cg = rip(adj)$C)
  ) 
}

.types <- function() return(c("fwd", "bwd", "tree", "tfwd"))
.types_msg <- function() {
  paste0("Types must be in one of ",
    paste0(paste0(.types()[-length(.types())], collapse = ", "), " or ", .types()[length(.types())])
  )  
}
