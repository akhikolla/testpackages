#' ciders data
#'
#' Case study pertaining to Quantitative Descriptive Analysis (QDA) applied to ten varieties of cider. The sensory panel consists of seven trained assessors who were asked to rate ten varieties of cider using a list of ten sensory attributes.
#'
#' @docType data
#'
#' @usage data(ciders)
#'
#' @format An object of class \code{"array"} with 10 ciders (mode 1), 10 sensory attributes (mode 2) and 7 assessors (mode 3):
#' \describe{
#' \item{ciders}{1 to 10}
#' \item{sensory attributes}{sweet, acid, bitter, astringency, odor strength, pungent, alcohol, perfume, intensity, and fruity }
#' \item{Panel}{Judge.1 to Judge.7 }
#' }
#' @references Ledauphin, S., Hanafi, M., & Qannari, E. M. (2006). Assessment of the agreement among the subjects in fixed vocabulary profiling. Food quality and preference, 17(3-4), 277-280.
#' (\href{https://doi.org/10.1016/j.foodqual.2005.03.017}{ScienceDirect})
#' @keywords datasets
#'
#'
#' @examples
#' data(ciders)
#' str(ciders)
"ciders"
