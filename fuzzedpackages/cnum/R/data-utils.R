## Function to return variable names of all character sets
return_vars <- function() {
  vars <- ls(envir = asNamespace("cnum"))
  i <- sapply(vars, function(x)
    is.character(get(x, envir = asNamespace("cnum"))), USE.NAMES = FALSE)
  vars[i]
}

## Function to return supported languages
return_langs <- function() {
  vars <- return_vars()
  vars <- vars[grepl("^scale_", vars)]
  vars <- gsub("^scale_[A-Za-z]+_([A-Za-z]+).*$", "\\1", vars)
  vars <- unique(vars)
  vars
}

## Function to return supported modes
return_modes <- function() {
  vars <- return_vars()
  vars <- vars[grepl("^scale_", vars)]
  vars <- gsub("^scale_([A-Za-z]+).*$", "\\1", vars)
  vars <- unique(vars)
  vars
}

## Fucntion to return a regular expression to match Chinse numerals
return_regex <- function(lang, mode, financial, only) {
  conv_t <- conv_table(lang, mode, financial)
  chr_c <- conv_t[["chr_t"]]$c
  scale_c <- conv_t[["scale_t"]]$c
  scale_c <- scale_c[scale_c != ""]
  zero <- conv_t[["zero"]]
  neg <- conv_t[["neg"]]
  dot <- conv_t[["dot"]]

  paste0(ifelse(only, "^", ""), neg, "?",
         "(", paste0(c(chr_c, scale_c, zero), collapse = "|"), ")+",
         dot, "+?", "(",
         paste0(c(chr_c, scale_c, zero), collapse = "|"), ")+",
         ifelse(only, "$|^", "|"), neg, "?",
         "(", paste0(c(chr_c, scale_c, zero), collapse = "|"), ")+")
}

## Function to return conversion table
conv_table <- function(lang, mode, financial) {
  if (!lang %in% return_langs())
    stop("unsupported language `", lang, "`.", call. = FALSE)

  if (!mode %in% return_modes())
    stop("unsupported mode `", mode, "`.", call. = FALSE)

  chr_var <- paste("chr", lang, sep = "_")
  scale_var <- paste("scale", mode, lang, sep = "_")
  interval_var <- paste("interval", mode, sep = "_")
  if (financial) {
    chr_var <- paste(chr_var, "f", sep = "_")
    scale_var <- paste(scale_var, "f", sep = "_")
  }
  zero_var <- paste("zero", lang, sep = "_")
  dot_var <- paste("dot", lang, sep = "_")
  neg_var <- paste("neg", lang, sep = "_")

  list(
    chr_t = data.frame(c = get(chr_var, envir = asNamespace("cnum")),
                       n = 1:9, stringsAsFactors = FALSE),
    scale_t = data.frame(c = get(scale_var, envir = asNamespace("cnum")),
                         n = get(interval_var, envir = asNamespace("cnum")),
                                 stringsAsFactors = FALSE),
    zero = get(zero_var, envir = asNamespace("cnum")),
    dot = get(dot_var, envir = asNamespace("cnum")),
    neg = get(neg_var, envir = asNamespace("cnum"))
  )
}
