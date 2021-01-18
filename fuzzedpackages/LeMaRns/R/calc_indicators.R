#' Calculate community indicators
#'
#' @description Calculates the Large Fish Indicator (LFI), Typical Length (TyL) or Length Quantile (LQ) for the community or a subset of the species.
#' @param wgt A matrix with dimensions \code{nsc} and \code{nfish} representing the weight of each species in each length class.
#' @param N A matrix with dimensions \code{nsc} and \code{nfish} representing the number of individuals in each length class for the current time step.
#' @param bry A numeric vector representing the length classes that are larger than \code{length_LFI}.
#' @param prop A numeric value between 0 and 1 representing how far along, as a proportion, the value of the large indicator threshold is in the length class that contains it.
#' @param u_bound A numeric vector of length \code{nsc} representing the upper bounds of the length classes.
#' @param prob A numeric value or vector between 0 and 1 denoting the LQ to be calculated.
#' @param mid A numeric vector of length \code{nfish} representing the mid-point of the length classes in the model.
#' @seealso \code{\link{get_indicators}}
#' @return \code{calc_LFI} returns a numeric value or vector representing the proportion of the biomass in the length classes above \code{bry} and in \code{prop} of the biomass in the length class \code{bry}.
#' @return \code{calc_TYL} returns a numeric value or vector representing the biomass-weighted geometric mean length of the community.
#' @return \code{calc_LQ} returns a numeric value or vector representing the length at which biomass exceeds a given proportion \code{prob} of the total biomass.
#' @examples
#' # Set up and run the model
#' NS_params <- LeMansParam(NS_par, tau=NS_tau, eta=rep(0.25, 21), L50=NS_par$Lmat, other=1e11)
#' effort <- matrix(0.5, 10, dim(NS_params@Qs)[3])
#' model_run <- run_LeMans(NS_params, years=10, effort=effort)
#'
#' # Calculate the LFI for 40cm
#' bry <- which(NS_params@l_bound<=40 & NS_params@u_bound>40)
#' length_LFI <- 40
#' prop <- (length_LFI-NS_params@l_bound[bry])/(NS_params@u_bound[bry]-NS_params@l_bound[bry])
#' LFI <- calc_LFI(model_run@N[,,101], NS_params@wgt, bry, prop)
#'
#' # Calculate TyL for the final time step
#' calc_TyL(wgt=NS_params@wgt, mid=NS_params@mid, N=model_run@N[,,101])
#'
#' # Calculate the LQ for the final time step
#' calc_LQ(wgt=NS_params@wgt, u_bound=NS_params@u_bound, N=model_run@N[,,101], prob=0.5)
#' @export
calc_LFI <- function(wgt, N, bry, prop){
  biomass <- N*wgt
  if (is.matrix(biomass)){
    biomass <- rowSums(biomass)
  }
  return((sum(biomass[(bry+1):length(biomass)])+prop*biomass[bry])/sum(biomass))
}

#' @rdname calc_LFI
#' @export
calc_LQ <- function(wgt, u_bound, N, prob){
  biomass <- N*wgt
  if (is.matrix(biomass)){
    biomass <- rowSums(biomass)
  }
  u_bound <- c(0, u_bound)
  probs <- c(0, cumsum(biomass)/sum(biomass))
  ov <- min(which(probs>prob))
  return(u_bound[ov-1]+(prob-probs[ov-1])/(probs[ov]-probs[ov-1])*(u_bound[ov]-u_bound[ov-1]))
}

#' @rdname calc_LFI
#' @export
calc_TyL <- function(wgt, mid, N){
  biomass <- N*wgt
  if (is.matrix(biomass)){
    biomass <- rowSums(biomass)
  }
  return(exp(sum(biomass*log(mid))/sum(biomass)))
}
