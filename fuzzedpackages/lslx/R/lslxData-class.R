## define R6 class \code{lslxData} to store data. ##
lslxData <-
  R6::R6Class(
    classname = "lslxData",
    public = list(
      index = "list",
      response = "list",
      pattern = "list",
      weight = "list",
      auxiliary = "list",
      sample_cov = "list",
      sample_mean = "list",
      sample_size = "list",
      sample_moment_acov = "list",
      numeric_variable = "character",
      ordered_variable = "character",
      weight_variable = "character",
      auxiliary_variable = "character",
      group_variable = "character",
      reference_group = "character",
      level_group = "character",
      nlevel_ordered = "numeric",
      name_response = "character"
    )
  )
