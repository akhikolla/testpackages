## Function to convert a Chinese numeral to integer then validate
c2integer <- function(number, conv_t) {
  converted <- c2integer_conv(number, conv_t)
  validate <- integer2c(format(converted, scientific = FALSE), conv_t)
  if (validate != paste0(number, collapse = ""))
    stop("`x` is not valid Chinese numeral.", call. = FALSE)
  converted
}

## Function to convert a Chinese numeral to decimal scale
c2decimal <- function(number, conv_t) {
  dot <- conv_t[["dot"]]
  number <- number[-grep(dot, number)]
  paste0(".", c2integer_literal(number, conv_t))
}

## Function to split up Chinese numerals
split_numeral <- function(number, conv_t, mode, financial) {
  chr_t <- conv_t[["chr_t"]]
  scale_t <- conv_t[["scale_t"]]
  zero <- conv_t[["zero"]]
  neg <- conv_t[["neg"]]
  dot <- conv_t[["dot"]]

  number_split <- strsplit(number, "")[[1]]
  # if multichar scale char exsists, check and re-merge, e.g. 万亿
  is_mchr <- nchar(scale_t$c) > 1
  if (any(is_mchr)) {
    mchr <- scale_t$c[is_mchr]
    if (any(stringr::str_detect(number, mchr))) {
      for (x in mchr) {
        index <- integer()
        # search through numerals for multichar
        i <- nchar(x) # starting index as length of multichar
        j <- i # hold value
        while (i <= length(number_split)) {
          # merge current and previous char
          merged <- paste0(number_split[(i - j + 1):i], collapse = "")
          if (merged == x)
            # if found multichar, save starting and ending index numbers
            index <- c(index, i - j + 1, i)
          i <- i + 1
        }
        if (!length(index)) next # skip if mulichar scale char not found
        i <- 1
        tmp <- character()
        while (i <= length(number_split)) {
          # merge multichar scale char while retaining order
          if (i %in% index) {
            # when i == starting index
            j <- index[grep(i, index) + 1] # ending index
            tmp <- c(tmp, paste0(number_split[i:j], collapse = "")) # merge
            i <- j + 1 # skip to element after ending index
          } else {
            tmp <- c(tmp, number_split[i])
            i <- i + 1
          }
        }
        number_split <- tmp
      }
    }
  }
  if (!all(number_split %in% c(chr_t$c, scale_t$c, zero, dot, neg))) {
    msg <- paste0("* `", number_split[!number_split %in% c(chr_t$c, scale_t$c, zero, dot, neg)],
                  "` is not a valid Chinese numeral character\n")
    if (financial)
      stop("`x` contains non-Chinese financial numerals in mode `", mode, "`:\n", msg, call. = FALSE)
    stop("`x` contains non-Chinese numerals in mode `", mode, "`:\n", msg, call. = FALSE)
  }
  number_split
}
