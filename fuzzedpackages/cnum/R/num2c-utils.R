## Function to convert integer to Chinese numeral with single scale character
integer2c_single <- function(number, conv_t) {
  scale_t <- conv_t[["scale_t"]]

  if (number %% 1 == 0)
    n <- nchar(format(number, scientific = FALSE))
  else {
    int <- gsub("\\..*$", "", number)
    n <- nchar(format(int, scientific = FALSE))
  }
  nearest_scale <- scale_t$n[max(which(scale_t$n <= n))]
  new_number <- number / 10^(nearest_scale - 1)

  if (new_number %% 1 == 0) {
    new_number <- format(number, scientific = FALSE)
    paste0(integer2c(new_number, conv_t), scale_t$c[scale_t$n == nearest_scale])
  } else {
    n <- nchar(gsub("^.*\\.", "", new_number)) # count deciaml places
    new_number <- format(new_number, scientific = FALSE, nsmall = ifelse(n > 22, 22, n))
    int <- gsub("\\..*$", "", new_number)
    dec <- gsub("^..*\\.", "", new_number)
    paste0(integer2c(int, conv_t), conv_t[["dot"]],
           integer2c_literal(dec, conv_t), scale_t$c[scale_t$n == nearest_scale])
  }
}
