## ---------------------------------------------------------
##                 NON-EXPORTED HELPERS
## ---------------------------------------------------------
parents <- function(g) {
  # out: the parents in a perfect ordering
  # g  : a pure discrete graph
  ps <- rip(g)$P
  structure(
    lapply(seq_along(ps), function(i) {
      setdiff(ps[[i]], names(ps)[i])
    }),
    names = names(ps)
  )
}

.parr <- function(pas, lvls, cell_rate = 0.5) {
  # pas: parents of a perfect ordering
  p_arr <- lapply(seq_along(pas), function(i) {
    p         <- names(pas)[i]
    pasi      <- pas[[i]]
    dim_names <- if (neq_empt_chr(pasi)) c(pasi, p) else p
    dim_vec   <- .map_int(lvls[dim_names], length)
    pi_arr    <- cbind(expand.grid(lvls[c(pasi, p)], stringsAsFactors = FALSE), prob = 0L)

    for (x in pasi) { # x in character(0) convieniently skips the loop
      pi_arr <- pi_arr[order(pi_arr[, x]), ]
    }

    reps_of_p_lvls <- 1:(nrow(pi_arr) / dim_vec[p])
    pi_arr[, "prob"] <- as.vector(
      sapply(reps_of_p_lvls, function(x) {
        px <- stats::rbeta(unname(dim_vec[p]), 2, cell_rate) # 2 is magic!
        px / sum(px)
      })
    )
    
    return(pi_arr)
  })
  structure(p_arr, names = names(pas))
}


#' Simulate observations from a decomposable graphical model
#'
#' @param g An adjacency list
#' @param lvls Named list with levels of the discrete variables
#' @param nsim Number of simulations
#' @param cell_rate Control discrete cell probabilities
#' @return This function returns a matrix of dimension where each row correspond
#' to a simulated observation from a DGM represented by \code{g}.
#' @examples
#'
#' g = list(
#'   A = c("B", "X", "Y"),
#'   B = c("A", "Y"),
#'   X = c("A", "Y"),
#'   Y = c("A", "X", "B")
#' )
#'
#' lvls <- list(
#'   A = c("0", "1"),
#'   B = c("0", "1"),
#'   X = c("a", "b", "c"),
#'   Y = c("0", "1", "2")
#' )
#'
#' dgm_sim_from_graph(g, lvls, nsim = 10)
#'
#' # Converting the g to a gengraph object and plot
#'
#' d <- data.frame(A = "", B = "", X = "", Y = "") # auxillary data
#' g <- gengraph(d, adj = g)
#' plot(g)
#'
#' @export
dgm_sim_from_graph <- function(g, lvls, nsim = 1000, cell_rate = 0.5) {

  pas  <- parents(g)
  parr <- .parr(pas, lvls, cell_rate)
  y    <- matrix("", nrow = nsim, ncol = length(names(g)), dimnames = list(NULL, names(g)))

  for (i in 1:nsim) {
    yi <- y[i, ]
    for (p in parr) {
      curr_var <- colnames(p)[(ncol(p)-1)]
      if (ncol(p) == 2L) {
        yi[curr_var] <- sample(p[, curr_var], 1L, prob = p[, "prob"])
      } else {
        curr_pas <- pas[[curr_var]]
        curr_val <- yi[curr_pas]
        for (cp in curr_pas) {
          p <- p[p[, cp, drop = TRUE] == curr_val[cp], ]
        }
        yi[curr_var] <- sample(p[, curr_var], 1L, prob = p[, "prob"]) 
      }
    }
    y[i, ] <- yi
  }
  return(y)
}

## #' Simulate observations from a decomposable graphical model given data
## #'
## #' @param A Character matrix or data frame (data)
## #' @param adj Adjacency list of a decomposable graph
## #' @param nsim Number of simulations
## #' @param ncores Number of cores to use in parallelization
## #' @return This function returns a matrix of dimension \code{nsim x ncol(A)} where each row correspond to a simulated observation from a DGM represented by \code{adj}.
## #' @examples
## #'
## #' g = list(
## #'   A = c("B", "X", "Y"),
## #'   B = c("A", "Y"),
## #'   X = c("A", "Y"),
## #'   Y = c("A", "X", "B")
## #' )
## #'
## #' d <- derma[, 1:4]
## #' colnames(d) <- names(g)
## #'
## #' dgm_sim_from_data(d, g)
## #'
## #' @export
## dgm_sim_from_data <- function(A, adj, nsim = 1000, ncores = 1) {
##   A   <- if (inherits(A, "data.frame")) as.matrix(A) else {stopifnot(is.matrix(A)); A}
##   adj <- if (inherits(adj, "gengraph")) adj_lst(adj) else adj
##   RIP <- ess::rip(adj) # the rip (or actually mcs) will check for decomposability here
##   cms <- a_marginals(A, RIP$C)
##   sms <- a_marginals(A, RIP$S)
##   out <- .sim_internal(A, cms, sms, nsim = nsim, type = "raw", ncores = ncores)
##   row.names(out) <- NULL
##   return(out)
## }
