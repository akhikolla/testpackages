#' @title Build a BLRM formula with linear interaction term in logit-space
#'
#' @description \code{blrm_formula_linear} is a convenience function for generating
#' a formula for \code{blrm_trial} and \code{blrm_exnex} with an interaction of the form:
#'
#' \deqn{\eta \, \prod_{i=1}^N (d_i \big / d_i^*))}
#'
#' @param ref_doses Numeric vector of reference doses with names corresponding to drug names
#' @param max_interaction_level Highest interaction order to consider [1 - Inf]. Default: 2
#'
#' @return The function returns an object of class \code{blrm_formula}.
#'
#' @examples
#'
#' ref_doses <- c(drug_A=10, drug_B=20)
#'
#' # can be used with blrm_trial
#' blrm_formula_linear(ref_doses)
#'
#' @export
blrm_formula_linear <- function(
  ref_doses,
  max_interaction_level = 2
)
{
    assert_int(max_interaction_level)

    assert_numeric(ref_doses, lower=0, finite=TRUE, any.missing=FALSE, names="strict")

    num_components <- length(ref_doses)
    component_names <- names(ref_doses)

  # Generate blrm_exnex formula
  blrm_formula <- "cbind(num_toxicities, num_patients - num_toxicities) ~ "

  # Add individual component terms
  for (component_index in seq(1, num_components))
  {
    blrm_formula <- paste0(
      blrm_formula,
      "1 + I(log(", component_names[component_index], "/",
      ref_doses[[component_names[component_index]]], ")) | "
    )
  }

  # Assemble interaction term
  blrm_formula <- paste0(blrm_formula, "0 ")

  num_interaction_terms <- 0

  if (num_components >= 2 && max_interaction_level >=2) {
    interaction_orders <- seq(2, min(num_components, max_interaction_level))
    for (num_interacting_components in interaction_orders)
    {
      for (interaction_compund_names in combn(component_names, num_interacting_components, simplify=FALSE))
      {
        blrm_formula <- paste0(blrm_formula, "+ I(")

        # Generate terms that are multiplied in the interaction
        i <- 0
        for (interaction_component_name in interaction_compund_names)
        {
          if (i > 0)
          {
            blrm_formula <- paste0(blrm_formula, ' * ')
          }
          blrm_formula <- paste0(blrm_formula,
            interaction_component_name ,"/", ref_doses[[interaction_component_name]]
          )
          i <- i + 1
        }

        blrm_formula <- paste0(blrm_formula, ") ")
      }
      num_interaction_terms <- num_interaction_terms + choose(num_components, num_interacting_components)
    }
  }

  blrm_formula <- paste0(blrm_formula, "| stratum_id / group_id")

  result <- list()
  result$blrm_formula <- blrm_formula

  result$num_interaction_terms <- num_interaction_terms
  result$num_components <- num_components

  result$component_names <- component_names

  result$max_interaction_level <- max_interaction_level
  result$ref_doses <- ref_doses

  structure(result, class="blrm_formula")
}
