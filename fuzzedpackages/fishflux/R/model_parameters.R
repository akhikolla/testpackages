#' A function to find a set of parameters
#'
#' @param sp Species name
#' @param family family
#' @param otolith TRUE or FALSE, if TRUE, function will only search fishbase for growth parameters that are based upon otolith analysis
#' @param temp temperature
#' @param ... Additional arguments to \code{\link{find_lw}}.
#' 
#' @details Returns a dataframe with all parameters that can be estimated
#' 
#' 
#' @examples
#' \donttest{
#' library(fishflux)
#' model_parameters(sp = "Scarus psittacus", family = "Scaridae", temp = 27)}
#' 
#' @export
model_parameters <- function(sp, family, otolith = TRUE, temp, ...) {
  #check species
  check_name_fishbase(sp)

  #dry_weight / wet_weight
  wprop <- wprop(family = family)

  #length weight
  length_weight <- find_lw(sp, ...)

  #growth parameters
  growth <- growth_params(sp = sp, otolith = otolith)

  #trophic level
  troph <- trophic_level(sp)

  #aspect ratio
  asp <- aspect_ratio(sp)

  #metabolism
  met <- metabolism(family = family, troph_m = troph$trophic_level,
                    temp = temp)

  data.frame(species  = sp,
             t0       = mean(growth$t0, na.rm = TRUE),
             Linf     = mean(growth$Linf, na.rm = TRUE),
             k        = mean(growth$k, na.rm = TRUE),
             asp      = asp$aspect_ratio,
             troph    = troph$trophic_level,
             lwa_m    = length_weight$lwa_m,
             lwa_sd   = length_weight$lwa_sd,
             lwb_m    = length_weight$lwb_m,
             lwb_sd   = length_weight$lwb_sd,
             mdw_m    = wprop$mdw,
             f0_m     = met$f0_m,
             f0_sd    = met$f0_sd,
             alpha_m  = met$alpha_m,
             alpha_sd = met$alpha_sd)
}
