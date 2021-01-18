#' Labeled Snowball with Multiple Inclusions (LSMI)
#'
#' Obtain LSMI samples around several seeds, which can be selected randomly or
#' pre-specified. See Figure 1 by
#' \insertCite{gel_etal_2017;textual}{snowboot} or
#' Figure 2 by \insertCite{chen_etal_2018_snowboot;textual}{snowboot}
#' illustrating the algorithm of sampling around multiple seeds.
#'
#' If \code{seeds} specified, \code{n.seed} is not used.
#'
#' @param net a network object that is a list containing:
#' \describe{
#'   \item{\code{degree}}{the degree sequence of the network, which is
#'      an \code{integer} vector of length \eqn{n};}
#'   \item{\code{edges}}{the edgelist, which is a two-column
#'      matrix, where each row is an edge of the network;}
#'   \item{\code{n}}{the network order (i.e., number of nodes in the network).}
#' }
#' The network object can be simulated by \code{\link{random_network}},
#' selected from the networks available in \code{\link{artificial_networks}},
#' converged from an \code{igraph} object with \code{\link{igraph_to_network}},
#' etc.
#' @param n.seed an integer defining the number of nodes to randomly sample
#' from the network to start an LSMI sample around each of them.
#' @param n.wave an integer defining the number of waves (order of the neighborhood)
#' to be recorded around the seed in the LSMI. For example, \code{n.wave = 1} corresponds to
#' an LSMI with the seed and its first neighbors. Note that the algorithm allows for
#' multiple inclusions.
#' @param seeds a vector of numeric IDs of pre-specified seeds. If specified,
#' LSMIs are constructed around each such seed.
#'
#'
#' @return A list of length \code{n.seed} (or, if \code{seeds} are specified,
#' of length \code{length(unique(seeds))}), where each element is a list
#' of length \code{n.wave + 1} representing an LSMI produced by
#' \code{\link{sample_about_one_seed}}.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{sample_about_one_seed}}, \code{\link{lsmi_union}}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' a <- lsmi(net, n.seed = 20, n.wave = 2)
#'
lsmi <- function(net, n.seed = 10, n.wave = 1, seeds = NULL) {
  if (is.null(seeds)) {
    seeds <- sort(sample(1:net$n, n.seed, replace = FALSE))
  } else {
    seeds <- sort(unique(seeds))
  }
  lapply(seeds, function(x) sample_about_one_seed(net, x, n.wave))
}
