#' Example AIL pedigree
#'
#' Example matrix describing the pedigree for advanced intercross lines
#'
#' @docType data
#'
#' @usage data(AILped)
#'
#' @format A data frame with five columns: individual id, mom, dad, sex (0
#' for females and 1 for males) and generation.
#'
#' @source Derived from the pedF8 dataset in the QTLRel package,
#' <https://cran.r-project.org/package=QTLRel>
#'
#' @keywords datasets
#'
#' @examples
#' data(AILped)
#' x <- sim_from_pedigree(AILped)
"AILped"
