#' Combines prey preference and prey suitability
#'
#' @description Calculates a combined value for prey preference and prey suitability standardised to a value between 0 and 1.
#' @param nsc A numeric value representing the number of length classes in the model.
#' @param nfish A numeric value representing the number of species in the model.
#' @param sc_Linf A numeric vector of length \code{nsc} representing the length class at which each species reaches its asymptotic length.
#' @param prefs An array of dimensions \code{nsc}, \code{nfish}, \code{nsc} and \code{nfish}. The first and second dimensions represent the prey whereas the the third and fourth dimensions represent the predator.
#' @param tau A matrix of dimensions \code{nfish} and \code{nfish}. Row indices represent predators and column indices represent prey. A value of 1 at location \code{i}, \code{j} indicates prey \code{j} is eaten by predator \code{i}.
#' @details \code{tau} values are assigned to an array of dimensions \code{nsc}, \code{nfish}, \code{nsc} and \code{nfish} and multiplied by the array \code{prefs}. This creates an array of dimensions \code{nsc}, \code{nfish}, \code{nsc} and \code{nfish} indicating prey suitability. Prey suitability is then standardised to sum to 1 for each predator species in each length class.
#' @return A list object of length \code{nfish}. Each element in the list is an array of dimensions \code{nsc}, \code{nsc} and \code{nfish} containing a value between 0 and 1 that represents prey preference and prey suitability for each species and length class.
#' @seealso \code{\link{calc_M2}}.
#' @references Hall, S. J., Collie, J. S., Duplisea, D. E., Jennings, S., Bravington, M., & Link, J. (2006). A length-based multispecies model for evaluating community responses to fishing. \emph{Canadian Journal of Fisheries and Aquatic Sciences}, 63(6):1344-1359.
#' @examples
#' # Set up the inputs to the function - species-independent parameters
#' nfish <- nrow(NS_par)
#' nsc <- 32
#' maxsize <- max(NS_par$Linf)*1.01 # the biggest size is 1% bigger than the largest Linf
#' l_bound <- seq(0, maxsize, maxsize/nsc); l_bound <- l_bound[-length(l_bound)]
#' u_bound <- seq(maxsize/nsc, maxsize, maxsize/nsc)
#' mid <- l_bound+(u_bound-l_bound)/2
#'
#' # Set up the inputs to the function - species-specific parameters
#' Linf <- NS_par$Linf # the von-Bertalanffy asymptotic length of each species (cm).
#' W_a <- NS_par$W_a # length-weight conversion parameter.
#' W_b <- NS_par$W_b # length-weight conversion parameter.
#' k <- NS_par$k # the von-Bertalnaffy growth parameter.
#' Lmat <- NS_par$Lmat # the length at which 50\% of individuals are mature (cm).
#'
#' # Get phi_min
#' tmp <- calc_phi(k, Linf, nsc, nfish, u_bound, l_bound, calc_phi_min=FALSE,
#'                   phi_min=0.1) # fixed phi_min
#' phi_min <- tmp$phi_min
#'
#' # Calculate growth increments
#' tmp <- calc_ration_growthfac(k, Linf, nsc, nfish, l_bound, u_bound, mid, W_a, W_b, phi_min)
#' sc_Linf <- tmp$sc_Linf
#' wgt <- tmp$wgt
#'
#' # Calculate predator-prey size preferences
#' prefs <- calc_prefs(pred_mu=-2.25, pred_sigma=0.5, wgt, sc_Linf)
#'
#' # Calculate prey preference and prey suitability
#' suit_M2 <- calc_suit_vect(nsc, nfish, sc_Linf, prefs, NS_tau)
#' @export
calc_suit_vect <- function(nsc, nfish, sc_Linf, prefs, tau) {
  # sc_Linf - calculated in calc_ration
  suit <- tau_temp <- array(0, dim = c(nsc, nfish, nsc, nfish))
  denom <- matrix(0, nsc, nfish)
  for (l in 1:nsc) {
    for (j in 1:nsc) {
      tau_temp[l,,j,] <- tau
    }
  }
  suit <- prefs*tau_temp

  # Standardize such that suitabilities sum to one
  # If the species is not a predator of anything, denom is zero
  for(sp1 in 1:nfish) {
    sc <- sc_Linf[sp1]
    for (sc1 in 1:sc) {
      #denom[sc1,sp1] <- sum(suit[sc1,sp1,,2:nfish])
      denom[sc1,sp1] <- sum(suit[sc1,sp1,,]) ## check this with Rob
      if (denom[sc1,sp1] > 0) {
        suit[sc1,sp1,,] <- suit[sc1,sp1,,]/denom[sc1,sp1]
      }
    }
  }

  suit[1,,,] <- 0
  suit_M2 <- lapply(seq_len(nfish), function(x) {return(suit[,x,,])})
  return(suit_M2)
}
