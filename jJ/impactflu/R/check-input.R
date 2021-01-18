# Input checking functions
# Arseniy Khvorov
# Created 2019/01/03
# Last edit 2019/01/03

#' Checks input for simulation and methods
#'
#' @param vaccinations Integer vector counts of vaccinations
#' @param cases Integer vector counts of cases
#' @param ve Numeric vector vaccine effectiveness (0 to 1).
#'
#' @importFrom rlang abort
#' @importFrom glue glue
#'
#' @noRd
check_counts <- function(vaccinations, cases, ve) {
  if (length(cases) == 1)
    abort("length of cases should be greater than 1")
  if (length(vaccinations) != length(cases))
    abort(glue(
      "length of cases ({length(cases)}) should match length of ",
      "vaccinations ({length(vaccinations)})"
    ))
  if (length(ve) != length(vaccinations))
    abort(glue(
      "length of ve ({length(ve)}) should be either 1 or the same as ",
      "the length of vaccinations and cases ({length(cases)})"
    ))
}
