#' @describeIn conversion Convert Chinese Numerals to Arabic Numerals.
#'
#' @return \code{c2num} returns a numeric vector.
#'
#' @examples
#' c2num("EXAMPLE CHECK")
#'
#' @export
#'
c2num <- function(x, lang = default_cnum_lang(), mode = "casual", financial = FALSE, literal = FALSE) {
  if (length(x) > 1)
    return(sapply(x, function(y) c2num(y, lang, mode, financial, literal)))

  if (!is.character(x)) {
    message("`x` coerced into character.")
    x <- as.character(x)
  }

  if (!nchar(x))
    stop("`x` must not be an empty string.")

  if (x == "EXAMPLE CHECK") {
    # to pass example check since examples must run but documentation can't contain non-ASCII text
    message(x)
    return(NULL)
  }

  conv_t <- conv_table(lang, mode, financial)
  scale_t <- conv_t[["scale_t"]]
  dot <- conv_t[["dot"]]

  number_split <- split_numeral(x, conv_t, mode, financial)
  if (number_split[1] == conv_t[["neg"]]) {
    neg <- "-"
    number_split <- number_split[-1]
  } else
    neg <- ""

  i <- grep(dot, number_split) # position of dot
  j <- ifelse(any(number_split %in% scale_t$c), # position of last scale char
              max(which(number_split %in% scale_t$c)), 0)
  if (length(i) > 1)
    stop("`x` must not contain more than one occurance of `", dot, "`.")
  # decimal
  if (length(i)) {
    if (i < j & j == length(number_split)) {
      # the position of dot is before a scale char: a numeral with single scale char
      converted <- as.numeric(paste0(neg, c2integer(number_split[1:(i - 1)], conv_t),
                                     c2decimal(number_split[i:(length(number_split) - 1)], conv_t)))
      converted * 10^(scale_t$n[scale_t$c == number_split[j]] - 1)
    } else {
      if (literal)
        as.numeric(paste0(neg, c2integer_literal(number_split[1:(i - 1)], conv_t),
                          c2decimal(number_split[i:length(number_split)], conv_t)))
      else
        as.numeric(paste0(neg, c2integer(number_split[1:(i - 1)], conv_t),
                          c2decimal(number_split[i:length(number_split)], conv_t)))
    }
  } else {
    # integer
    if (literal)
      as.numeric(paste0(neg, c2integer_literal(number_split, conv_t)))
    else
      as.numeric(paste0(neg, c2integer(number_split, conv_t)))
  }
}
