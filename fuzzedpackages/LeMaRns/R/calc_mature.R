#' Calculate the proportion of mature individuals
#'
#' @description Calculates the proportion of individuals that are mature for each species and length class.
#' @param Lmat A numeric vector of length \code{nfish} representing the length at which 50\% of individuals are mature for each species.
#' @param nfish A numeric value representing the number of species in the model.
#' @param mid A numeric vector of length \code{nfish} representing the mid-point of the length classes in the model.
#' @param kappa A numeric vector of length \code{nfish} representing the rate of change from immaturity to maturity for each species.
#' @param sc_Linf A numeric vector of length \code{nsc} representing the length class at which each species reaches its asymptotic length.
#' @param eps A numeric value specifying a numerical offset. The default is \code{1e-5}.
#' @param force_mature A logical statement indicating whether to force all fish in the largest length class to be mature. The default is \code{TRUE}.
#' @details The proportion of individuals in the \code{j}th length class of the \code{i}th species that are mature is described by a logistic model
#' @details \code{1/(1+exp(-kappa[i]*(mid[j]-Lmat[i])))}
#' @return A matrix with dimensions \code{nsc} and \code{nfish} and elements in the range 0-1 representing the proportion of individuals that are mature for each species and length class.
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
#' phi <- tmp$phi
#' phi_min <- tmp$phi_min
#'
#' # Calculate growth increments
#' tmp <- calc_ration_growthfac(k, Linf, nsc, nfish, l_bound, u_bound, mid, W_a, W_b, phi_min)
#' sc_Linf <- tmp$sc_Linf
#'
#' # Calculate the proportion of mature individuals
#' mature <- calc_mature(Lmat, nfish, mid, kappa=rep(10, nfish), sc_Linf)
#' @export
calc_mature <- function(Lmat, nfish, mid, kappa, sc_Linf, eps=1e-5, force_mature=TRUE) {
  mature <- t(outer(1:nfish, mid, function(x, y, kappa, Lmat){
    return(1/(1+exp(-kappa[x]*(y-Lmat[x]))))}, kappa=kappa, Lmat=Lmat))

  if (force_mature){
    mature[cbind(sc_Linf, 1:nfish)] <- 1
  }

  mature[mature<eps] <- 0
  return(mature)
}

