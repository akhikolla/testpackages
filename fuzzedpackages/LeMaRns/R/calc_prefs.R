#' Calculate predator-prey size preferences
#'
#' @description Calculates the size preference of each predator species in each length class for each prey species in each length class.
#' @param pred_mu A numeric value representing the preferred predator-prey mass ratio.
#' @param pred_sigma A numeric value representing the width of the weight preference function.
#' @param wgt A matrix with dimensions \code{nsc} and \code{nfish} representing the weight of each species in each length class.
#' @param sc_Linf A numeric vector of length \code{nsc} representing the length class at which each species reaches its asymptotic length.
#' @details A predator of species \code{i} in length class \code{j} has a size preference for species \code{k} in length class \code{l} equal to
#' @details \code{exp(-(log10(wgt[l, k]/wgt[j, i])-pred_mu)^2/(2*pred_sigma))}.
#' @return An array of dimensions \code{nsc}, \code{nfish}, \code{nsc} and \code{nfish}. The first and second dimensions represent the prey species whereas the third and fourth dimensions represent the predator species.
#' @references Hall, S. J., Collie, J. S., Duplisea, D. E., Jennings, S., Bravington, M., & Link, J. (2006). A length-based multispecies model for evaluating community responses to fishing. \emph{Canadian Journal of Fisheries and Aquatic Sciences}, 63(6):1344-1359.
#' @references Ursin, E. (1973). On the prey size preferences of cod and dab. \emph{Meddelelser fra Danmarks Fiskeri-og Havundersgelser}, 7:85-98.
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
#' @export
calc_prefs <- function(pred_mu, pred_sigma, wgt, sc_Linf) {
  prefs <- outer(wgt, wgt, function(x, y) {exp(-(log10(x/y)-pred_mu)^2/(2*pred_sigma))})
  prefs <- aperm(prefs, c(3, 4, 1, 2))
  prefs[is.na(prefs)] <- 0
return(prefs)
}
