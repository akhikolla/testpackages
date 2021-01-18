#' Calculate growth
#'
#' @description Calculates the number of individuals of each species in each length class for the next time step.
#' @param N A matrix with dimensions \code{nsc} and \code{nfish} representing the number of individuals in each length class for the current time step.
#' @param phi A matrix with dimensions \code{nsc} and \code{nfish} representing the proportion of individuals that leave each length class.
#' @param nfish A numeric value representing the number of species in the model.
#' @param nsc A numeric value representing the number of length classes in the model.
#' @return A matrix with dimensions \code{nsc} and \code{nfish} representing the number of individuals of each species in each length class for the next time step.
#' @examples
#' # Set up the inputs to the function - species-independent parameters
#' nfish <- nrow(NS_par)
#' nsc <- 32
#' maxsize <- max(NS_par$Linf) * 1.01 # the biggest size is 1% bigger than the largest Linf
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
#' phi <- tmp$phi
#' phi_min <- tmp$phi_min
#'
#' # Calculate growth increments
#' tmp <- calc_ration_growthfac(k, Linf, nsc, nfish, l_bound, u_bound, mid, W_a, W_b, phi_min)
#' sc_Linf <- tmp$sc_Linf
#' wgt <- tmp$wgt
#'
#' # Get an initial population
#' N0 <- get_N0(nsc, nfish, mid, wgt, sc_Linf, intercept=1e10, slope=-5)
#'
#' # Calculate growth
#' growth <- calc_growth(N0, phi, nfish, nsc)
#' @export
calc_growth <- function(N, phi, nfish, nsc) {
  stay <- (1-phi)*N
  leave <- phi*N
  res <- stay
  res[2:nsc, ] <- res[2:nsc, ] + leave[1:nsc-1, ]
  return(res)
}
