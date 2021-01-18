#' Threshold estimate
#' @description Returns the estimate of the threshold.
#'
#' @param gfi an output of \code{\link{gfigpd2}}
#'
#' @return The estimated threshold.
#' @export
thresholdEstimate <- function(gfi){
  attr(gfi, "threshold")
}

#' Join MCMC chains
#' @description Joins multiple MCMC chains into a single chain.
#'
#' @param gfi an output of \code{\link{gfigpd1}} or \code{\link{gfigpd2}} 
#'   containing more than one chain
#'
#' @return A \code{mcmc} object.
#' @export
#' @importFrom stats start
joinMCMCchains <- function(gfi){
  if(inherits(gfi, "mcmc.list")){
    coda::mcmc(as.matrix(gfi), start = start(gfi), thin = coda::thin(gfi))
  }else{
    if(inherits(gfi, "mcmc")){
      gfi
    }else{
      stop(
        "Not a `mcmc.list` object."
      )
    }
  }
}
