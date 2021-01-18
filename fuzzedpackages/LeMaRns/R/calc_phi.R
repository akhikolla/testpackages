#' Calculate the proportion of individuals that leave each length class for each species each length class each time step.
#'
#' @description Calculate the proportion of individuals of each species that leave each length class at each time step.
#' @param k A numeric vector of length \code{nfish} representing the von Bertalanffy growth parameter \code{(1/yr)} for each species.
#' @param Linf A numeric vector of length \code{nfish} representing the asymptotic length of each species.
#' @param nsc A numeric value representing the number of length classes in the model.
#' @param nfish A numeric value representing the number of species in the model.
#' @param u_bound A numeric vector of length \code{nsc} representing the upper bounds of the length classes.
#' @param l_bound A numeric vector of length \code{nsc} representing the lower bounds of the length classes.
#' @param calc_phi_min A logical statement indicating whether \code{phi_min} should be calculated within the function. The default is \code{FALSE}.
#' @param phi_min A fixed numeric value of \code{phi_min}, which represents the time step of the model. This is only required if \code{calc_phi_min=FALSE}. The default is \code{0.1}.
#' @details Calculates the time (yrs) for an average fish to grow from the lower to the upper bound of a length class assuming von Bertalanffy growth. The values are scaled to the fastest growing fish and length class combination in order to calculate the proportion of individuals leaving each length class in a time step.
#' @return A list object containing \code{phi} and \code{phi_min}. \code{phi} is a matrix of dimensions \code{nsc} and \code{nfish} representing the proportion of individuals of each species that leave each length class. \code{phi_min} is a numeric value representing the time step of the model.
#' @references Hilborn, R. & Walters, C.J. (1992). Quantitative Fisheries Stock Assessment. Springer.
#' @references von Bertalanffy, L. (1957). Quantitative Laws in Metabolism and Growth. \emph{The Quarterly Review of Biology}, 32:217-231.
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
#' # Calculate the proportion of individuals that leave each length class
#' # with and without a fixed value for phi_min
#' tmp <- calc_phi(k, Linf, nsc, nfish, u_bound, l_bound, calc_phi_min=FALSE,
#'                   phi_min=0.1) # fixed phi_min
#' phi <- tmp$phi
#' phi_min <- tmp$phi_min
#'
#' tmp <- calc_phi(k, Linf, nsc, nfish, u_bound, l_bound, calc_phi_min=TRUE) # without fixed phi_min
#' phi <- tmp$phi
#' phi_min <- tmp$phi_min
#' @export
calc_phi <- function(k, Linf, nsc, nfish, u_bound, l_bound, calc_phi_min=FALSE, phi_min=0.1) {
  L_low <- matrix(l_bound, nsc, nfish)
  L_up <- matrix(u_bound, nsc, nfish)
  ktemp <- t(matrix(k, nfish, nsc))
  Ltemp <- t(matrix(Linf, nfish, nsc))

  phi <- suppressWarnings((1/ktemp)*log((Ltemp-L_low)/(Ltemp-L_up)))
  phi[phi<0|is.na(phi)] <- 0

  if (calc_phi_min) {
    phi_min <- min(phi[phi>0])
  }

  if(any(phi[phi>0] < phi_min)){
    warning(paste("phi_min, ",phi_min,", is larger than the smallest phi,",min(phi[phi>0]),".",  sep=""))
  }

  phi[phi>0] <- phi_min/phi[phi>0]
  return(list(phi=phi, phi_min=phi_min))
}
