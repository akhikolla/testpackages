#' @describeIn conversion Convert Arabic Numerals to Chinese Numerals.
#'
#' @return \code{num2c} returns a character vector.
#'
#' @examples
#' num2c(721)
#' num2c(-6)
#' num2c(3.14)
#' num2c(721, literal = TRUE)
#' num2c(1.45e12, financial = TRUE)
#' num2c(6.85e12, lang = "sc", mode = "casualPRC")
#' num2c(1.5e9, mode = "SIprefix", single = TRUE)
#'
#' @export
#'
num2c <- function(x, lang = default_cnum_lang(), mode = "casual", financial = FALSE, literal = FALSE, single = FALSE) {
  if (length(x) > 1)
    return(sapply(x, function(y) num2c(y, lang, mode, financial, literal, single)))

  if (!is.numeric(x))
    stop("`x` must be numeric.")

  conv_t <- conv_table(lang, mode, financial)
  if (x < 0) {
    neg_chr <- conv_t[["neg"]]
    x <- gsub("-", "", x)
    x <- as.numeric(x)
  } else
    neg_chr <- ""
  # single scale char
  if (single & x >= 11)
    paste0(neg_chr, integer2c_single(x, conv_t))
  else {
    if (x %% 1 == 0) {
      # integer
      x <- format(x, scientific = FALSE)
      if (literal)
        paste0(neg_chr, integer2c_literal(x, conv_t))
      else
        paste0(neg_chr, integer2c(x, conv_t))
    } else {
      # decimal
      n <- nchar(gsub("^.*\\.", "", x)) # count deciaml places
      x <- format(x, scientific = FALSE, nsmall = ifelse(n > 20, 20, n))
      int <- gsub("\\..*$", "", x)
      dec <- gsub("^..*\\.", "", x)
      if (literal)
        paste0(neg_chr, integer2c_literal(int, conv_t),
               conv_t[["dot"]], integer2c_literal(dec, conv_t))
      else
        paste0(neg_chr, integer2c(int, conv_t),
               conv_t[["dot"]], integer2c_literal(dec, conv_t))
    }
  }
}
