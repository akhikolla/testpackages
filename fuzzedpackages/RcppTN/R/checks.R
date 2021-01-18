checkInputs <- function(m, s, l, h) {
    if (!(length(m) == length(s) &
          length(m) == length(l) &
          length(m) == length(h)
          )
        ) {
        stop("Input vectors not all same length. Nothing done.")
    }
}

checkOutputs <- function(out) {
    if (any(is.na(out))) {
        warning("NAs returned. Check for invalid parameters.")
    }
}
