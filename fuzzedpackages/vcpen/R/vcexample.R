#' Example data for Penalized Variance Component method
#'
#' Datasets for an example run of vcpen with 4 variance components calculated as kernel matrices from genotype dosage (dose) on 100 subjects with two covariates (covmat), and a continuous response.
#'
#' @format The example contains three data.frames and a response vector for 100 subjects at 70 SNPs accross 4 variance components:
#'   \describe{
#'     \item{\code{covmat}}{two arbitrary covariates (columns) for 100 subjects (rows)}
#'     \item{\code{dose}}{genotype dosage at 70 SNPs (columns) and 100 subjects (rows)}
#'     \item{\code{doseinfo}}{2-column matrix with indices for grouping SNPs into variance components (for Kernel Matrix)}
#'     \item{\code{response}}{continuous response vector for 100 subjects}
#'   }
#' @examples
#' data(vcexample)
#' dim(dose)
#' dim(doseinfo)
#' dim(covmat)
#' length(response)
#' @name vcexample
NULL
#> NULL

#' @rdname vcexample
#' @name covmat
NULL
#> NULL
#' @rdname vcexample
#' @name dose
NULL
#> NULL
#' @rdname vcexample
#' @name doseinfo
NULL
#> NULL
#' @rdname vcexample
#' @name response
NULL
#> NULL
