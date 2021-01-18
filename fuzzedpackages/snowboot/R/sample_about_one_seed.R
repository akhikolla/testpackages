#' Snowball Sampling with Multiple Inclusions around One Network Node
#'
#' This function obtains a labeled snowball with multiple inclusions (LSMI) sample,
#' starting from a single network node called seed. See Figure 1 by
#' \insertCite{thompson_etal_2016;textual}{snowboot} illustrating the algorithm
#' of sampling around one seed.
#'
#' @inheritParams lsmi
#' @param seed numeric ID of a seed to start the LSMI.
#'
#' @return \code{sample_about_one_seed} returns a list of length \code{n.wave + 1}
#' containing ID of the seed (1st element of the output list), IDs of nodes in the
#' 1st wave (2nd element of the list), \ldots, IDs of nodes in the wave \code{n.wave}
#' ((\code{n.wave + 1})th element of the list). If a wave has no nodes in it, the
#' corresponding element of the output contains \code{NA}.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{lsmi}}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' a <- sample_about_one_seed(net, seed = 1, n.wave = 2)
#'
sample_about_one_seed <- function(net, seed, n.wave = 1) {
  if (n.wave < 1) {stop("Number of waves, n.wave, should be >= 1.")}
  effEdges <- net$edges
  nodes.waves <- as.list(c(seed, rep(NA, n.wave)))
  wave <- 1
  while (wave <= n.wave & nrow(effEdges) >= 0) {
    tmp <- is.element(effEdges, nodes.waves[[wave]])
    if (any(tmp)) {
      tmp <- which(matrix(tmp, dim(effEdges)[1], 2), arr.ind = TRUE)
      nodes.waves[[wave + 1]] <- sort(effEdges[cbind(tmp[, 1], sapply(tmp[, 2], FUN = switch, 2, 1))])
      effEdges <- effEdges[-tmp[, 1], ]
      if (is.vector(effEdges)) {
        effEdges <- t(effEdges)
      }
    }
    wave <- wave + 1
  }
  nodes.waves
}
