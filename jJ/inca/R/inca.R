#' @importFrom stats model.matrix
#' @importFrom Matrix Matrix
#' @importFrom Rcpp evalCpp
#' @useDynLib inca
#'
#' @name inca-package
#' @aliases inca
#' 
#' @docType package
#' 
#' @encoding UTF-8
#' 
#' @title
#' Integer Calibration
#' 
#' @description
#' Specific functions are provided for rounding real weights to integers and performing integer programming algorithms for calibration problems.
#' 
#' @details
#' \tabular{ll}{
#' Package: \tab inca\cr
#' Type: \tab Package\cr
#' Version: \tab 0.0.4\cr
#' Date: \tab 2019-09-18\cr
#' License: \tab GPL (>= 2)\cr
#' }
#' 
#' Calibration forces the weighted estimates of calibration variables to match known totals. 
#' This improves the quality of the design-weighted estimates. It is used to adjust for 
#' non-response and/or under-coverage. The commonly used methods of calibration produce 
#' non-integer weights. In cases where weighted estimates must be integers, one must 
#' "integerize" the calibrated weights. However, this procedure often produces final weights 
#' that are very different for the "sample" weights. To counter this problem, the \pkg{inca} 
#' package provides specific functions for rounding real weights to integers, and performing 
#' an integer programming algorithm for calibration problems with integer weights.
#'
#' For a complete list of exported functions, use \code{library(help = "inca")}.
#' 
#' This research was supported in part by the U.S. Department of Agriculture, 
#' National Agriculture Statistics Service. The findings and conclusions in this
#' publication are those of the authors and should not be construed to represent
#' any official USDA or U.S. Government determination or policy.
#' 
#' @author
#' Luca Sartore \email{luca.sartore@@usda.gov} and Kelly Toppin \email{kelly.toppin@@nass.usda.gov}
#'
#' Maintainer: Luca Sartore \email{drwolf85@@gmail.com}
#' 
#' @references
#' Theberge, A. (1999). Extensions of calibration estimators in survey sampling. \emph{Journal of the American Statistical Association}, \bold{94}(446), 635-644.
#'
#' Little, R. J., & Vartivarian, S. (2003). On weighting the rates in non-response weights.
#'
#' Kish, L. (1992). Weighting for unequal Pi. \emph{Journal of Official Statistics}, \bold{8}(2), 183.
#'
#' Rao, J. N. K., & Singh, A. C. (1997). A ridge-shrinkage method for range-restricted weight calibration in survey sampling. \emph{In Proceedings of the section on survey research methods} (pp. 57-65). American Statistical Association Washington, DC.
#'
#' Horvitz, D. G., & Thompson, D. J. (1952). A generalization of sampling without replacement from a finite universe. \emph{Journal of the American Statistical Association}, \bold{47}(260), 663-685.
#'
#' Kalton, G., & Flores-Cervantes, I. (2003). Weighting methods. \emph{Journal of Official Statistics}, \bold{19}(2), 81-98.
#'
#' Sartore, L., Toppin, K., Young, L., Spiegelman, C. (2019). Developing integer calibration weights for the Census of Agriculture. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, \bold{24}(1), 26-48.
#' 
#' @keywords
#' integer calibration rounding
#' 
#' @examples
#' library(inca)
#' 
NULL
