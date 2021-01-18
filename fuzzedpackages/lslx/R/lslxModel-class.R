## define R6 class \code{lslxModel} to store model information. ##
lslxModel <-
  R6::R6Class(
    classname = "lslxModel",
    public = list(
      model = "character",
      numeric_variable = "character",
      ordered_variable = "character",
      weight_variable = "character",
      auxiliary_variable = "character",
      group_variable = "character",
      reference_group = "character",
      level_group = "character",
      nlevel_ordered = "numeric",
      name_response = "character",
      name_factor = "character",
      name_eta = "character",
      name_exogenous = "character",
      name_endogenous = "character",
      name_threshold = "character",
      specification = "data.frame"
    )
  )
