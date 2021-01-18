#' Chinese Numerals Detection and Extraction
#'
#' Functions to detect and extract Chinese numerals in character object and
#' string.
#'
#' @param x the character object or string to be tested or to extract from.
#'
#' @param strict logical: Should the Chinese numerals format be strictly
#'   enforced? A casual test only checks if \code{x} contains Chinese numerals
#'   characters. A strict test checks if \code{x} is valid Chinese numerals.
#'   (e.g. "\emph{yi bai yi}" will pass the casual test and fail the strict
#'   test)
#'
#' @param ... optional arguments to be passed to \code{\link[base]{grepl}} (for
#'   \code{is_cnum} and \code{has_cnum}) or
#'   \code{\link[stringr]{str_extract_all}} (for \code{extract_cnum}).
#'   Disregarded when \code{strict = TRUE}.
#'
#' @inheritParams conversion
#'
#' @inheritSection conversion Details
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
#' @seealso \link[=conversion]{Functions for conversion}
#'
#' @name tools
#' @order 1
#'
NULL

#' @describeIn tools Test if character object is Chinese numerals. A wrapper
#'   around \code{\link[base]{grepl}}.
#'
#' @return \code{is_cnum} returns a logical vector indicating is Chinese
#'   numerals or not for each element of \code{x}).
#'
#' @examples
#' is_cnum("yibai ershiyi")
#'
#' @export
#'
is_cnum <- function(x, lang = default_cnum_lang(), mode = "casual", financial = FALSE,
                    literal = FALSE, strict = FALSE, ...) {
  if (strict) {
    if (length(x) > 1)
      return(sapply(x, function(y) is_cnum(y, lang, mode, financial, literal, strict, ...)))
    tryCatch(as.logical(c2num(x, lang, mode, financial, literal)),
             error = function(cnd) FALSE)
  } else
    grepl(return_regex(lang, mode, financial, TRUE), x, ...)
}

#' @describeIn tools Test if string contains Chinese numerals. A wrapper around
#'   \code{\link[base]{grepl}}.
#'
#' @return \code{has_cnum} returns a logical vector indicating contains Chinese
#'   numerals or not for each element of \code{x}.
#'
#' @examples
#' has_cnum("yibai bashi yuan")
#'
#' @export
#'
has_cnum <- function(x, lang = default_cnum_lang(), mode = "casual", financial = FALSE, ...) {
  if (length(x) > 1)
    return(sapply(x, function(y) has_cnum(y, lang, mode, financial, ...)))
  grepl(return_regex(lang, mode, financial, FALSE), x, ...)
}

#' @describeIn tools Extracts Chinese numerals from string. A wrapper around
#'   \code{\link[stringr]{str_extract_all}} from \code{stringr}.
#'
#' @return \code{extract_cnum} returns a list of character vectors containing
#'   the extracted Chinese numerals.
#'
#' @examples
#' extract_cnum("shisiyi ren")
#'
#' @export
#'
extract_cnum <- function(x, lang = default_cnum_lang(), mode = "casual", financial = FALSE, ...) {
  stringr::str_extract_all(x, return_regex(lang, mode, financial, FALSE), ...)
}
