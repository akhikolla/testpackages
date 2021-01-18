#' cnum: Working with Chinese Numerals
#'
#' This R package provides useful functions to work with Chinese numerals in R,
#' such as conversion between Chinese numerals and Arabic numerals as well as
#' detection and extraction of Chinese numerals in character objects and string.
#'
#' @section Warnings: This package supports conversion of numbers with absolute
#'   value not greater than 1e+18. Note that numbers in R are in double
#'   precision that carries approximately 16 significant digits. The conversion
#'   accuracy for numbers beyond this limit is therefore not guaranteed.
#'
#' @note Due to technical limitation, R package documentation cannot contain any
#'   non-ASCII characters. Therefore, Chinese characters are represented in
#'   romanized Chinese \emph{pinyin} in the documentation. Visit the GitHub page
#'   for examples in Chinese.
#'
#' @author Elgar Teo (\email{elgarteo@@connect.hku.hk})
#'
#' @seealso GitHub page: \url{https://github.com/elgarteo/cnum}
#'
#' @docType package
#' @keywords internal
#' @importFrom stringr str_detect str_extract_all
#' @importFrom Rcpp sourceCpp
#' @useDynLib cnum, .registration = TRUE
#'
"_PACKAGE"
