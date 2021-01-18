#' Calculate growth increments
#'
#' @description Calculates the amount of food required for fish of a given species and length class to grow according to the von Bertalanffy growth curve in a time step.
#' @param k A numeric vector of length \code{nfish} representing the von Bertalanffy growth parameter \code{(1/yr)} for each species.
#' @param Linf A numeric vector of length \code{nfish} representing the asymptotic length of each species.
#' @param nsc A numeric value representing the number of length classes in the model.
#' @param nfish A numeric value representing the number of species in the model.
#' @param l_bound A numeric vector of length \code{nsc} representing the lower bounds of the length classes.
#' @param u_bound A numeric vector of length \code{nsc} representing the upper bounds of the length classes.
#' @param mid A numeric vector of length \code{nfish} representing the mid-point of the length classes.
#' @param W_a A numeric vector of length \code{nfish} representing the parameter \code{a} in the length-weight conversion. See 'Details' for more information.
#' @param W_b A numeric vector of length \code{nfish} representing the parameter \code{b} in the length-weight conversion. See 'Details' for more information.
#' @param phi_min A numeric value representing the time step of the model.
#' @param vary_growth A logical statement indicating whether growth efficiency should vary for each species (\code{vary_growth=TRUE}) or be fixed at the value given by \code{fixed_growth} (\code{vary_growth=FALSE}). The default is \code{FALSE}.
#' @param growth_eff If \code{vary_growth==TRUE}, \code{growth_eff} is a numeric representing the growth efficiencies of a fish of length 0. If \code{vary_growth==FALSE}, \code{growth_eff} is a numeric value of length \code{1} representing a fixed growth efficiency for all fish. The default is 0.5.
#' @param growth_eff_decay A numeric value specifying the rate at which growth efficiency decreases as length approaches \code{Linf}. The default is 0.11.
#' @details The weight increments of the \code{i}th species in the \code{j}th length class is calculated by determining the amount an individual will grow in one time step, \code{phi_min}, if it were to follow the von Bertalanffy growth curve
#' @details \code{L22=(Linf[i]-mid[j])*(1-exp(-k[i]*phi_min))}.
#' @details The weight of a fish at the mid-point of the size class is calculated using the length-weight relationship
#' @details \code{wgt[j,i] = a[i]*mid[j]^b[i]},
#' @details and similarly the expected change in weight of the the fish is calculated as
#' @details \code{growth_inc = (W_a[i]*L22^W_b[i])}.
#' @details It also has a growth efficiency
#' @details  \code{g_eff[j, i]=(1-(wgt[j,i]/(W_a[i]*Linf[i]^W_b[i]))^growth_eff_decay)*growth_eff}
#' @details if \code{vary_growth==TRUE} or \code{g_eff[j, i]=growth_eff} otherwise.
#' @details \code{ration} is then calculated by
#' @details \code{growth_inc*(1/g_eff[j, i])}.
#' @return A list object containing \code{ration}, \code{sc_Linf}, \code{wgt} and \code{g_eff}. \code{ration} is a matrix with dimensions \code{nsc} and \code{nfish} representing the amount of food required for fish of a given species and length class to grow according to the von Bertalanffy growth curve in a time step. \code{sc_Linf} is a numeric vector of length \code{nfish} representing the length class at which each species reaches \code{Linf}. \code{wgt} is a matrix with dimensions \code{nsc} and \code{nfish} representing the weight of each species in each length class. \code{g_eff} is a matrix with dimensions \code{nsc} and \code{nfish} representing the growth efficiency of each species in each length class.
#' @references von Bertalanffy, L. (1957). Quantitative Laws in Metabolism and Growth. \emph{The Quarterly Review of Biology}, 32:217-231
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
#' ration <- tmp$ration
#' sc_Linf <- tmp$sc_Linf
#' wgt <- tmp$wgt
#' g_eff <- tmp$g_eff
#' @export
calc_ration_growthfac <- function(k, Linf, nsc, nfish, l_bound, u_bound, mid, W_a, W_b, phi_min, vary_growth=TRUE, growth_eff=0.5, growth_eff_decay=0.11) {
  # Check growth_eff_decay
  if (growth_eff_decay<=0) {
    stop("growth_eff_decay cannot be less than or equal to 0. Please select a value greater than 0 but less than or equal to 1.")
  } else if (growth_eff_decay>1) {
    warning("A growth_eff_decay of greater than 1 is not meaningful. Please select a value greater than 0 but less than or equal to 1.")
  }

  ration <- wgt <- g_eff <- matrix(0, nsc, nfish)
  sc_Linf <- rep(0, nfish)

  sc_Linf <- sapply(1:nfish, function(x, l_bound, Linf) {
    max(which(l_bound<Linf[x]))}, l_bound=l_bound, Linf=Linf)

  for (i in 1:nfish) {
    for (j in 1:sc_Linf[i]) {
      # Calculate the length after 1 time step when the initial length is the mid-point of a length class
        L22 <- mid[j]+(Linf[i]-mid[j])*(1-exp(-k[i]*phi_min))
        wgt[j, i] <- (W_a[i]*mid[j]^W_b[i])
        growth_inc <- (W_a[i]*L22^W_b[i])-wgt[j, i]

        if (vary_growth) {
          g_eff[j, i] <- (1-(wgt[j,i]/(W_a[i]*Linf[i]^W_b[i]))^growth_eff_decay)*growth_eff
          } else {
          g_eff[j, i] <- growth_eff
        }
        ration[j, i] <- growth_inc*(1/g_eff[j, i])
      }

    fmid <- l_bound[sc_Linf[i]]+((Linf[i]-l_bound[sc_Linf[i]])/2)
    L2 <- fmid+(Linf[i]-fmid)*(1-exp(-k[i]*phi_min))
    wgt[sc_Linf[i], i] <- (W_a[i]*fmid**W_b[i])
    growth_inc <- (W_a[i]*L2**W_b[i])-wgt[sc_Linf[i], i]
    g_eff[sc_Linf[i], i] <- 1-(wgt[sc_Linf[i], i]/(W_a[i]*Linf[i]**W_b[i]))**growth_eff_decay
    g_eff[sc_Linf[i], i] <- g_eff[sc_Linf[i], i]*0.5
    ration[sc_Linf[i], i] <- growth_inc*(1/g_eff[sc_Linf[i], i])
  }

  return(list(ration=ration, sc_Linf=sc_Linf, wgt=wgt, g_eff=g_eff))
}
