#' @useDynLib HDPenReg, .registration = TRUE
#' @import rtkore
#' @import methods
#' @import Matrix
#' @importFrom graphics abline axis lines points matplot
#' @importFrom stats rbeta rbinom rpois
#'
#' @title Algorithms for lasso and fused-lasso problems.
#' @docType package
#' @aliases HDPenReg-package, HDPenReg
#' @name HDPenReg-package
#' @description This package contains algorithms for lasso and fused-lasso problems.
#' It contains an implementation of the lars algorithm [1],
#' for the lasso and fusion penalization and EM-based algorithms for (logistic) lasso and fused-lasso.
#'
#' @details
#'
#' \tabular{ll}{
#' Package: \tab HDPenReg\cr
#' Type: \tab Package\cr
#' Version: \tab 0.94.5\cr
#' Date: \tab 2019-03-29\cr
#' License: \tab GPL (>=2) \cr
#' }
#'
#' The main function is \link{HDlars}.
#'
#'
#' @author Maintainer: Quentin Grimonprez <quentin.grimonprez@@inria.fr>
#'
#' @examples
#' \dontrun{
#' # see vignette
#' vignette("HDPenReg")
#' }
#'
#' @seealso \code{\link{HDlars}} \code{\link{HDcvlars}}
#'
#' @keywords package
NULL
