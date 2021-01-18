

#' Calculate network statistics from a formula
#' @param formula A lolog formula (See \code{\link{lolog}}).
#' @examples
#' data(ukFaculty)
#' calculateStatistics(ukFaculty ~ edges + mutual + triangles)
calculateStatistics <- function(formula) {
  createCppModel(formula, cloneNet = FALSE)$statistics()
}