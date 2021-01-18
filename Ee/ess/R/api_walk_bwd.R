#' Stepwise backward selection
#' @description Stepwise backward selection in decomposable graphical models
#' @param x \code{gengraph}
#' @param df data.frame
#' @param q Penalty term in the stopping criterion  (\code{0} = AIC and \code{1} = BIC)
#' @param thres A threshold mechanism for choosing between two different ways of calculating
#' the entropy. Can Speed up the procedure with the "correct" value.
#' @details A \code{bwd} object can be created using the \code{gengraph} constructor with \code{type = "bwd"}
#' @return A \code{bwd} object; a subclass of \code{gengraph}) used for backward selection.
#' @examples
#'
#' d <- derma[, 10:25]
#'
#' g <- gengraph(d, type = "bwd")
#' s <- walk(g, d)
#' print(s)
#' plot(s)
#' adj_lst(s)
#' adj_mat(s)
#' 
#' @seealso \code{\link{fit_graph}}, \code{\link{walk.fwd}}, \code{\link{gengraph}}
#' @export
walk.bwd <- function(x, df, q = 0.5, thres = 5) {
  nodes   <- names(x$G_adj)
  M       <- nrow(df)
  penalty <- log(M)*q + (1 - q)*2 
  e_min   <- Inf
  for (i in 1:ncol(x$G_A)) {
    for (j in 1:i) {
      if (x$G_A[i, j] == 1L) {
        pair <- c(nodes[i], nodes[j])
        pair_clique_idx <- 1L
        pair_in_cliques <- sapply(seq_along(x$CG), function(z) {
          pair_in_z <- all(pair %in% x$CG[[z]])
          if (pair_in_z) pair_clique_idx <<- z
          pair_in_z
        })
        if (sum(pair_in_cliques) < 2) { # Eligeble for deletion
          C     <- x$CG[[pair_clique_idx]]
          va    <- pair[1]
          vb    <- pair[2]
          Ca    <- setdiff(C, vb)
          Cb    <- setdiff(C, va)
          S     <- intersect(Ca, Cb)
          sp    <- sort_(pair)
          ed    <- entropy_difference(sp, S, df, x$MEM, thres)
          x$MEM       <- ed$mem
          HM_HM_prime <- ed$ent
          dev         <- 2*M*HM_HM_prime
          d_parms     <- -prod(x$LV[pair] - 1) * prod(x$LV[S])
          d_qic       <- dev + penalty * d_parms
          if (d_qic <= e_min) {
            e_min <- d_qic
            x$e   <- structure(c(va, vb), "d_qic" = d_qic)
            x$S   <- S
          }
        }
      }
    }  
  }
  x$G_adj[[x$e[1]]] <- setdiff(x$G_adj[[x$e[1]]], x$e[2])
  x$G_adj[[x$e[2]]] <- setdiff(x$G_adj[[x$e[2]]], x$e[1])
  del_idx <- match(x$e, colnames(x$G_A))
  x$G_A[del_idx[1], del_idx[2]] <- 0L
  x$G_A[del_idx[2], del_idx[1]] <- 0L
  x$CG <- rip(x$G_adj, check = FALSE)$C
  return(x)
}
