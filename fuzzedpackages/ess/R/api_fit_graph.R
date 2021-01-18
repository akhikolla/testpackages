#' Fit a decomposable graphical model
#' @description A generic method for structure learning in decomposable graphical models
#' @param df Character data.frame
#' @param type Character ("fwd", "bwd", "tree" or "tfwd")
#' @param adj Adjacency list of a decomposable graph
#' @param q Penalty term in the stopping criterion (\code{0} = AIC and \code{1} = BIC)
#' @param trace Logical indicating whether or not to trace the procedure
#' @param thres A threshold mechanism for choosing between two different ways of calculating the entropy.
#' @param wrap logical specifying if the result of a run with type = "tree" should be converted to a "fwd" object
#' @return A \code{gengraph} object representing a decomposable graph.
#' @examples
#'
#' g <- fit_graph(derma, trace = FALSE, q = 0)
#' print(g)
#' plot(g)
#'
#' # Adjacency matrix and adjacency list
#' adjm <- adj_mat(g)
#' adjl <- adj_lst(g)
#'
#' # Cliques in the graph
#' rip(adjl)$C
#'
#' # Components of the graph
#' components(adjl) # only one here
#' 
#' @details
#' The types are
#' \itemize{
#' \item "fwd": forward selection
#' \item "bwd": backward selection
#' \item "tree": Chow-Liu tree (first order interactions only)
#' \item "tfwd": A combination of "tree" and "fwd". This can speed up runtime considerably in high dimensions.
#' }
#' Using \code{adj_lst} on an object returned by \code{fit_graph} gives the adjacency list corresponding to the graph. Similarly one can use \code{adj_mat} to obtain an adjacency matrix. Applying the \code{rip} function on an adjacency list returns the cliques and separators of the graph.
#' @references \url{https://arxiv.org/abs/1301.2267}, \url{https://doi.org/10.1109/ictai.2004.100} 
#' @seealso \code{\link{fit_components}}, \code{\link{adj_lst.gengraph}}, \code{\link{adj_mat.gengraph}},
#' \code{\link{walk.fwd}}, \code{\link{walk.bwd}}, \code{\link{gengraph}}
#' @export
fit_graph <- function(df,
                      type  = "fwd",
                      adj   = NULL,
                      q     = 0.5,
                      trace = TRUE,
                      thres = 5,
                      wrap  = TRUE)
{

  n        <- ncol(df)
  complete <- n * (n-1L) / 2L
  null     <- 0L

  if (!is.data.frame(df)) stop("df must be a data.frame.")
  
  if (!(type %in% .types())) stop(.types_msg())
  
  if (q < 0 || q > 1) stop("q must be between 0 and 1")
  
  if (n == 1L) {
    adj <- structure(list(character(0)), names = colnames(df))
    x   <- gengraph(df, type = "gen", adj)
    return(x)
  } 

  x <- gengraph(df, type, adj)

  if (inherits(x, "fwd")) {
    if (!neq_empt_chr(as.vector(x$e))) {
      # If no edges are added in fwd_init, x$e = character(0)
      if (trace) msg(0L, complete, 0L, "delta-qic")
      return(x)
    }
  }  
  
  if (inherits(x, "tree")) return(fit_tree(x, df, wrap))
    
  triv     <- trivial(x, null, complete)
  update_k <- update_iteration(x)
  k        <- sum(x$G_A)/2

  x <- walk(x = x, df = df, q = q, thres = thres)
  k <- update_k(k)

  if (k == triv) return(x)

  stp      <- stop_condition(x)
  stop_val <- attr(x$e, "d_qic")
  if (stp(stop_val)) return(x)

  while (!stp(stop_val)) {
    if (trace) msg(k, complete, stop_val, "delta-qic")
    x <- walk(x = x, df = df, q = q, thres = thres)
    k  <- update_k(k)
    if (k == triv) {
      if (trace) msg(k, complete, stop_val, "delta-qic")
      return(x)
    }
    stop_val <- attr(x$e, "d_qic")
  }
  if (trace) msg(k, complete, stop_val, "delta-qic")
  return(x)
}


#' Fit a decomposable graphical model on each component
#' @description Structure learning in decomposable graphical models on several components
#' @param df data.frame
#' @param comp A list with character vectors. Each element in the list is a component in the graph (using expert knowledge)
#' @param type Character ("fwd", "bwd", "tree" or "tfwd")
#' @param q Penalty term in the stopping criterion (\code{0} = AIC and \code{1} = BIC)
#' @param as_gen Logical. Convert to gengraph or not. If true, the graph can be plotted.
#' @param trace Logical indicating whether or not to trace the procedure
#' @param thres A threshold mechanism for choosing between two different ways of calculating the entropy.
#' @param wrap logical specifying if the result of a run with type = "tree" should be converted to a "fwd" object
#' @return A \code{gengraph} object
#' @seealso \code{\link{fit_graph}}, \code{\link{adj_lst.gengraph}}, \code{\link{adj_mat.gengraph}}, \code{\link{walk.fwd}},
#' \code{\link{walk.bwd}}, \code{\link{gengraph}}
#' @export
fit_components <- function(df,
                      comp,
                      type   = "fwd",
                      q      = 0.5,
                      as_gen = TRUE,
                      trace  = TRUE,
                      thres  = 5,
                      wrap   = TRUE)
{
  adj <- lapply(unname(comp), function(x) {
    fit_graph(df[, x, drop = FALSE], type = type,  q = q, trace = trace, thres = thres, wrap = wrap)
  })
  CG <- NULL
  if (as_gen) CG <- unlist(lapply(adj, function(x) x$CG), recursive = FALSE)
  adj <- unlist(lapply(adj, adj_lst), recursive = FALSE)
  if (as_gen) adj <- new_gengraph(df, adj, CG)
  return(adj)
}
