#' Bootstrapping Empirical Degree Distribution
#'
#' This function delivers bootstrap estimates of network degree distribution
#' based on an LSMI sample. The bootstrap scheme is non-weighted for seeds
#' (resampling with replacement) and weighted for non-seeds (resampling with
#' replacement, with weights proportional to inverse of the degrees),
#' as described in Section 3.3 by \insertCite{thompson_etal_2016;textual}{snowboot}
#' and in Algorithm 1 by \insertCite{gel_etal_2017;textual}{snowboot}.
#'
#' @param x a list that is the output of \code{\link{lsmi_dd}}, i.e., an estimate
#' of the degree distribution together with all degrees of seeds and non-seeds
#' from an LSMI.
#' @param B a positive integer, the number of bootstrap replications to perform.
#' Default is 100.
#' @param cl parameter to specify computer cluster for bootstrapping, passed to
#' the package \code{parallel} (default is \code{1}, meaning no cluster is used).
#' Possible values are:
#' \itemize{
#'   \item cluster object (list) produced by \link[parallel]{makeCluster}.
#'   In this case, new cluster is not started nor stopped;
#'   \item \code{NULL}. In this case, the function will attempt to detect
#'   available cores (see \link[parallel]{detectCores}) and, if there are
#'   multiple cores (\eqn{>1}), a cluster will be started with
#'   \link[parallel]{makeCluster}. If started, the cluster will be stopped
#'   after computations are finished;
#'   \item positive integer defining the number of cores to start a cluster.
#'   If \code{cl = 1}, no attempt to create a cluster will be made.
#'   If \code{cl > 1}, cluster will be started (using \link[parallel]{makeCluster})
#'   and stopped afterwards (using \link[parallel]{stopCluster}).
#' }
#'
#' @return A list object of class "\code{snowboot}" consisting of:
#'    \item{fkb}{A matrix of dimensions \code{length(x$fk)}\eqn{\times}\code{B}
#'    with \code{B} bootstrap estimates of the degree distribution.
#'    The bootstrap estimates are computed according to
#'    Equation 1 by \insertCite{gel_etal_2017;textual}{snowboot}, also
#'    see \insertCite{chen_etal_2018_snowboot;textual}{snowboot}.}
#'    \item{mub}{A vector of length \code{B} with bootstrapped estimates
#'    of the network mean degree.
#'    The bootstrap estimates are computed according to
#'    Equation 2 by \insertCite{gel_etal_2017;textual}{snowboot}.}
#'    \item{fk}{A vector with an estimate of the degree distribution, copied
#'    from the input \code{x$fk}.}
#'    \item{mu}{An estimate of the mean degree, copied from the input \code{x$mu}.}
#'    \item{B}{The number of bootstrap replications performed.}
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{lsmi}}, \code{\link{lsmi_dd}}, \code{\link{boot_ci}}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' lsmiEstimate <- lsmi_dd(net = net, n.seed = 5, n.wave = 3)
#' bootEstimates <- boot_dd(lsmiEstimate, B = 10)
#'
boot_dd <- function(x, B = 100, cl = 1) {
  bootparallel <- FALSE
  if (is.list(cl)) {
    bootparallel <- TRUE
    clStop <- FALSE #cl is supplied externally, we will not stop it
  } else {
    if (is.null(cl)) {
      cores <- parallel::detectCores()
    } else {
      cores <- cl
    }
    if (cores > 1) {
      bootparallel <- TRUE
      cl <- parallel::makeCluster(cores)
      clStop <- TRUE #cl is created internally, we will stop it later
    }
  }
  #x is the output of lsmi_dd
  ds <- x$ds
  dns <- x$dns
  ns <- length(ds)
  nns <- length(dns)
  fk <- x$fk
  dmax <- length(fk) - 1
  k <- c(0:dmax)
  if (bootparallel) {
    #degree distribution of bootstrapped seeds
    ddsb <- parallel::parSapply(cl, X = 1:B, FUN = function(b)
      table(c(sample(ds, replace = TRUE), k)) - 1
    )
    #degree distribution of bootstrapped non-seeds
    ddnsb <- boot_dd_ns_cl(dns, k[-1], B, cl)
    if(clStop) { parallel::stopCluster(cl) }
  } else {
    #degree distribution of bootstrapped seeds
    ddsb <- sapply(1:B, FUN = function(b)
      table(c(sample(ds, replace = TRUE), k)) - 1
    )
    #degree distribution of bootstrapped non-seeds
    ddnsb <- boot_dd_ns(dns, k[-1], B)
  }
  fkb <- (ddsb[-1,] + sweep(ddnsb, MARGIN = 2, (1 - ddsb[1,]), `*`)) / (ns + nns)
  fkb <- rbind(ddsb[1,], fkb)
  rownames(fkb) <- k
  mub <- apply(sweep(fkb, MARGIN = 1, k, `*`), 2, sum)
  res <- list(fkb = fkb, mub = mub, fk = fk, mu = x$mu, B = B)
  class(res) <- "snowboot"
  return(res)
}


