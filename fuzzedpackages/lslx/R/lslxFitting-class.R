## define R6 class \code{lslxFitting} to store fitting result. ##
lslxFitting <-
  R6::R6Class(
    classname = "lslxFitting",
    public = list(
      control = "list",
      reduced_model = "list",
      reduced_data = "list",
      supplied_result = "list",
      fitted_result = "list"
    )
  )

