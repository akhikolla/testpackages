#' A function to calculate metabolic rates
#'
#' All model parameters below were estimated by Barneche & Allen 2018 Ecology
#' Letters doi: 10.1111/ele.12947. These parameters are for the best model
#' (Model 2 in the paper online supplementary material) of fish resting
#' metabolic rates reported in the paper, which also includes trophic level as
#' a covariate.
#' 
#' @param temp  Temperature in degrees Celsius
#' @param m_max Maximum biomass fish (in g)
#' @param m     Wet weight fish (in g)
#' @param asp   The caudal fin aspect ratio , a proxy for activity level
#' @param troph Trophic level (from 1 to 5)
#' @param f     Activity scope (from 1 to 4)
#' @param a     Resting metabolic rate mass-scaling exponent
#' @param B0    Constant for resting metabolic rate. If NA, function will calculate an average.
#' @param growth_g_day   Daily growth in grams of wet weight
#' 
#' @details     All model parameters below were estimated by 
#' Barneche & Allen 2018 Ecology Letters doi: 10.1111/ele.12947. 
#' These parameters are for the best model 
#' (Model 2 in the paper online supplementary material) 
#' of fish resting metabolic rates reported in the paper, 
#' which also includes trophic level as a covariate.
#' 
#' @keywords fish metabolism
#' 
#' @examples
#' library(fishflux)
#' fishflux::metabolic_rate(temp = 27, m_max = 600, m = 300, asp = 3,
#' troph = 2, f = 2, growth_g_day = 0.05, B0 = 0.2, a = 0.6 )
#' 
#' @export
metabolic_rate <- function(temp, troph, asp, B0, m_max, m, a, growth_g_day, f) {

  # the cost of growth (Em, in J / g)
  Em <- exp(4.38 + 0.1032 * log(temp) + 0.73 * log(troph) + 0.41 * log(asp + 1))

  # metabolic rates, B, are given in g C / day, while growth rates, G, are given in g / day,
  # so B / G results g C / g.
  # Em values, however, are expressed in Joules / g, so we need to convert g C to Joules.
  # This conversion should be calculated based on the standard enthalpy combustion of glucose,
  # i.e. a biochemical argument which corresponds to metabolism --
  # it includes both the energy to be sequestrated by ATP and the energy lost during this process.
  # (2805 kJ / 1 mole Glucose) * (1 mole of Glucose / 6 moles C) * (1 mole C / 12 g C) = 2805 / 6 / 12 approx. 39 kJ per g C

  gC_to_J  <- 39e3

  Ec       <- 24e3 # combustion energy of biomass (Joules / g)
  Bm       <- B0 * gC_to_J * m_max ^ (a - 1)
  B_main   <- Bm * m  # maintenance metabolic rate
  B_syn    <- Em * growth_g_day # cost of growth
  B_rest   <- B_main + B_syn # resting metabolic rate (Joules / day)
  B_tot    <- B_rest * f # assimilation rate in joule per day
  Cm       <- B_tot / gC_to_J # amount of mass C needed for metabolism

  data.frame(Total_metabolic_rate_C_g_d   = Cm,
             Total_metabolic_rate_j_d     = B_tot,
             Resting_metabolic_rate_j_d   = B_rest,
             Cost_growth_j_g              = B_syn)
}
