#' Cross-validation to Select an Optimal Combination of n.seed and n.wave
#'
#' From the vector of specified \code{n.seeds} and possible waves \code{1:n.wave} around each
#' seed, the function selects a single number \code{n.seed} and an \code{n.wave}
#' (optimal seed-wave combination) that produce
#' a labeled snowball with multiple inclusions (LSMI) sample with desired
#' bootstrap confidence intervals for a parameter of interest. Here by `desired'
#' we mean that the interval (and corresponding seed-wave combination) are selected
#' as having the best coverage (closest to the specified level \code{prob}), based on
#' a cross-validation procedure with proxy estimates of the parameter.
#' See Algorithm 2 by \insertCite{gel_etal_2017;textual}{snowboot} and Details
#' below.
#'
#' @inherit boot_ci details
#'
#' @inheritParams lsmi
#' @inheritParams lsmi_union
#' @inheritParams boot_dd
#' @inheritParams boot_ci
#' @param param The parameter of interest for which to run a cross-validation
#' and select optimal \code{n.seed} and \code{n.wave}. Currently, only one
#' selection is possible: \code{"mu"} (the network mean degree).
#' @param proxyRep The number of times to repeat proxy sampling. Default is 19.
#' @param proxySize The size of the proxy sample. Default is 30.
#'
#' @return A list consisting of:
#' \item{bci}{A numeric vector of length 2 with the bootstrap confidence interval
#' (lower bound, upper bound) for the parameter of interest. This interval is
#' obtained by bootstrapping node degrees in an LSMI with the optimal combination
#' of \code{n.seed} and \code{n.wave}
#' (the combination is reported in \code{best_combination}).}
#' \item{estimate}{Point estimate of the parameter of interest
#' (based on the LSMI with \code{n.seed} seeds and \code{n.wave} waves
#' reported in the \code{best_combination}).}
#' \item{best_combination}{An integer vector of lenght 2 containing the optimal
#' \code{n.seed} and \code{n.wave} selected via cross-validation.}
#' \item{seeds}{A vector of numeric IDs of the seeds that were used
#' in the LSMI with the optimal combination of \code{n.seed} and \code{n.wave}.}
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{lsmi}}, \code{\link{lsmi_union}}, \code{\link{boot_dd}}, \code{\link{boot_ci}}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' a <- lsmi_cv(net, n.seeds = c(10, 20, 30), n.wave = 5, B = 100)
#'
lsmi_cv <- function(net, n.seeds, n.wave, seeds = NULL,
                    B = 100, prob = 0.95, cl = 1,
                    param = c("mu"),
                    method = c("percentile", "basic"),
                    proxyRep = 19, proxySize = 30)
{
  method <- match.arg(method)
  param <- match.arg(param)
  param <- ifelse(param == "mu", 2, 1) #mu is the 2nd output in the functions, 1st is fk
  clStop <- FALSE #cl is supplied externally, we will not stop it
  if (!is.list(cl)) {
    if (is.null(cl)) {
      cores <- parallel::detectCores()
    } else {
      cores <- cl
    }
    if (cores > 1) {
      cl <- parallel::makeCluster(cores)
      clStop <- TRUE #cl is created internally, we will stop it later
    }
  }
  patch <- lsmi_union(net, n.seeds, n.wave, seeds)
  subpatch <- subsubpatch <- patch$lsmi_big
  ss <- patch$sequence_seeds
  ns <- length(ss)
  if (proxySize < 0.7*length(ss[[1]])) {
    #use only seeds for proxies
    proxyReplace <- FALSE
    used <- ss[[1]]
  } else {
    #use all unique LSMI nodes for proxies
    proxyReplace <- TRUE
    used <- unique(na.omit(unlist(patch$lsmi_big)))
  }
  used <- net$degree[used]
  CIs <- estimate <- as.list(rep(NA, ns * n.wave))
  counter <- 1
  for (i in 1:ns) {
    if (i > 1) {
      subpatch <- subsubpatch <- subpatch[is.element(ss[[i - 1]], ss[[i]])]
    }
    for (j in 1:(n.wave)) {
      if (j > 1) {
        subsubpatch <- lapply(subsubpatch, function(x) x[1:(n.wave - j + 1)])
      }
      dd <- lsmi_dd(subsubpatch, net)
      bdd <- boot_dd(dd, B, cl = cl)
      CIs[[counter]] <- boot_ci(bdd, prob, method)[[param]]
      estimate[[counter]] <- dd[[param]]
      counter <- counter + 1
    }
  }
  if (param == 2) { #calculate proxies for mu
    proxies <- sapply(1:proxyRep, function(x) mean(sample(used, proxySize, replace = proxyReplace)) )
  }
  proxiCoverage <- sapply(CIs, function(x) mean(x[1] <= proxies & proxies <= x[2]))
  SWcombinations <- cbind(rep(1:ns, each = n.wave), rep(n.wave:1, ns))
  #if several options are chosen above, by default, prefer the one with more seeds and more waves:
  best_comb_index <- which.min(abs(proxiCoverage - prob))[1]
  SeedWave <- SWcombinations[best_comb_index,]
  seeds <- ss[[SeedWave[1]]]
  best_combination <- SeedWave
  best_combination[1] <- length(seeds)
  names(best_combination) <- c("n.seed", "n.wave")
  if(clStop) { parallel::stopCluster(cl) }
  list(bci = CIs[[best_comb_index]],
       estimate = estimate[[best_comb_index]],
       best_combination = best_combination,
       seeds = seeds)
}
