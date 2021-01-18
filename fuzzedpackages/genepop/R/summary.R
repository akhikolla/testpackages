summary.GP_outfile <- function(object, ...) { ## object is a filename with added "GP_outfile" class
  splits <- strsplit(object, ".", fixed = TRUE)[[1L]]
  termin <- splits[length(splits)]
  switch(termin,
         "NUL" = .summary_GP_NUL(object, ...), ## i.e. outfile ....NUL for null alleles
         # Add more case here...
         NULL)
  invisible(object)
}

.summary_GP_NUL <- function(object, what=c("point est", "lower bound"), ...) {
  resu <- readLines(object)
  what <- what[1L]
  if (what=="point est") {
    pointests <- as.numeric(unlist(lapply(strsplit(resu[grep(" Null  ",resu)]," Null"),paste,collapse="")))
    return(pointests)
  } else if (what=="lower bound") {
    ruleslines <- grep("============",resu)
    nrules <- length(ruleslines)
    CIlines <- resu[(ruleslines[nrules-1L]+4L):(ruleslines[nrules]-1L)]
    infolines <- CIlines[grep("(No info for CI)", CIlines,invert = TRUE)]
    lowerbounds <- rep(NA_real_, length(CIlines))
    for (it in seq_along(CIlines)) {
      CIline <- CIlines[[it]]
      if ( ! length(grep("(No info for CI)",CIline))) {
        subline <- substr(CIline,start=33,stop=1000)
        lowerbounds[it] <- as.numeric(stringr::str_extract_all(subline,pattern="\\(?[0-9,.]+\\)?")[[1]][1])
      }
    }
    return(lowerbounds)
  }
}