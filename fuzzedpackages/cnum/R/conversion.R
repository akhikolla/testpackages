#' Chinese Numerals Conversion
#'
#' Functions to convert between Chinese and Arabic numerals.
#'
#' @param x the Arabic/Chinese numerals to be converted, or a vector of them.
#'   The absolute value must not be greater than 1e+18.
#'
#' @param lang the language of the Chinese numerals. \code{"tc"} for Traditional
#'   Chinese. \code{"sc"} for Simplified Chinese. The default is \code{"tc"},
#'   but this can be changed by setting \code{\link[base]{options}(cnum.lang =
#'   "sc")}.
#'
#' @param mode the scale naming system to be enforced. See the ‘Details’ section
#'   for the list of supported modes.
#'
#' @param financial logical: should the financial numerals be used (\emph{daxie
#'   shuzi})?
#'
#' @param literal logical: should the numerals be converted literally? (e.g. 721
#'   to be converted to "\emph{qi er yi}" instead of "\emph{qibai ershiyi}" and
#'   vice versa)
#'
#' @param single logical: should the return result with one scale character
#'   only? (e.g. 1.5e+08 as "\emph{yi dian wuyi}" instead of "\emph{yiyi
#'   wuqianwan}")
#'
#' @section Details: The following scale naming systems are supported: \itemize{
#'   \item {\code{"casual"}: }{the casual naming system used outside of mainland
#'   China, i.e. 1e+09 is referred to as "\emph{yi zhao}".} \item
#'   {\code{"casualPRC"}: }{the casual naming system used in mainland China,
#'   i.e. 1e+9 is referred to as "\emph{yi wanyi}".} \item {\code{"SIprefix"}:
#'   }{the SI prefix system used in Taiwan as stipulated in the document
#'   \emph{Names, Definitions and Symbols of the Legal Units of Measurement and
#'   the Decimal Multiples and Submultiples}.} \item{\code{"SIprefixPRC"}: }{the
#'   SI prefix system used in mainland China as stipulated in the document
#'   \emph{China Statutory Measurement Units}.} \item{\code{"SIprefixPRClong"}:
#'   }{ a variant of \code{"SIprefixPRC"} with long prefixes, e.g. 1e+09 is
#'   referred to as "\emph{yi jika}" instead of "\emph{yi ji}".} }
#'
#' @inheritSection cnum-package Warnings
#' @section Warnings: The modes \code{"casual"} and \code{"casualPRC"}
#'   implements a “myriad scale” with an interval of 1e+04 for large numbers,
#'   i.e. "\emph{yi}" is 10,000 times of "\emph{wan}", which is different from
#'   some of the interval systems used in ancient Chinese writings.
#'
#' @references The standard for mode \code{"SIprefix"} \emph{Names, Definitions
#'   and Symbols of the Legal Units of Measurement and the Decimal Multiples and
#'   Submultiples} is available from
#'   \url{https://gazette.nat.gov.tw/egFront/detail.do?metaid=108965} (in
#'   Traditional Chinese).
#'
#'   The standard for mode \code{"SIprefixPRC"} \emph{China Statutory
#'   Measurement Units} is available from
#'   \url{http://gkml.samr.gov.cn/nsjg/jls/201902/t20190225_291134.html} (in
#'   Simplified Chinese).
#'
#' @seealso \link[=tools]{Functions for detetction and extraction}
#'
#' @name conversion
#' @order 1
#'
NULL
