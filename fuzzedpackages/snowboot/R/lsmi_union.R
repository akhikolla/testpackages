#' Snowball Sampling with Multiple Inclusions around Several Subsets of Seeds
#'
#' Obtain one big LSMI -- with \code{max(n.seeds)} seeds and \code{n.wave}
#' waves around each -- and subsample seeds to create smaller LSMIs (with less
#' seeds and/or waves). The function is primarily used in cross-validation.
#'
#' Note that the produced LSMIs are slightly different from those described by
#' \insertCite{gel_etal_2017;textual}{snowboot}. The current R implementation
#' produces smaller LSMIs by subsetting the seeds, not by new sampling of
#' seeds from the network and growing completely new LSMIs, as it was done by
#' \insertCite{gel_etal_2017;textual}{snowboot}. See the details in Figure 3 by
#' \insertCite{chen_etal_2018_snowboot;textual}{snowboot}
#'
#' @inheritParams lsmi
#' @param n.seeds an integer vector of numbers of seeds for snowball sampling
#' (cf. a single integer \code{n.seed} in \code{\link{lsmi}}). Only
#' \code{n.seeds <= n} are retained. If \code{seeds} is
#' specified, only values \code{n.seeds < length(unique(seeds))} are retained
#' and automatically supplemented by \code{length(unique(seeds))}.
#'
#' @return A list with two elements:
#'   \item{lsmi_big}{LSMI with \code{max(n.seeds)} seeds (see the argument definition
#'   above) and \code{n.wave} waves produced by the \code{\link{lsmi}} function.}
#'   \item{sequence_seeds}{A list of length equal to \code{length(n.seeds)};
#'   each element of the list is a random subset of the seeds' IDs, starting from
#'   the largest (a set of size \code{max(n.seeds)}) to the smallest
#'   (a set of size \code{min(n.seeds)}).}
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{sample_about_one_seed}}, \code{\link{lsmi}}, \code{\link{lsmi_cv}}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' a <- lsmi_union(net, n.seeds = c(5, 10, 15), n.wave = 2)
#'
lsmi_union <- function(net, n.seeds, n.wave, seeds = NULL){
  #Make this function sample just the biggest LSMI, and remember indices for subsamples.
  n.seeds <- n.seeds[n.seeds <= net$n]
  if (is.null(seeds)) {
    max_seeds <- max(n.seeds)
    seeds <- sample(1:net$n, max_seeds, replace = FALSE)
  } else {
    seeds <- unique(seeds)
    if (length(seeds) > net$n) {
      seeds <- seeds[1:net$n]
    }
    max_seeds <- length(seeds)
    n.seeds <- c(n.seeds[n.seeds < max_seeds], max_seeds)
  }
  seeds <- sort(seeds)
  n.wave <- round(max(n.wave))
  n.seeds <- sort(unique(n.seeds), decreasing = TRUE)
  #Get an LSMI for the combination of max_seeds and largest n.wave
  lsmi_big <- lsmi(net, n.wave = n.wave, seeds = seeds)
  sequence_seeds <- list(seeds)
  if (length(n.seeds) > 1){
    for (i in 2:length(n.seeds)) {
      sequence_seeds[[i]] <- sort(sample(sequence_seeds[[i - 1]], n.seeds[i], replace = FALSE))
    }
  }
  list(lsmi_big = lsmi_big, sequence_seeds = sequence_seeds)
}
