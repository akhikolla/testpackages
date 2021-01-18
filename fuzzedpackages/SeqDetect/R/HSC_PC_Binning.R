HSC_PC_Binning <- function(min_value, max_value, bins, value_field) {
  if(!is.numeric(min_value))
    stop("HSC_PC binning: min_value must be of numeric value")
  if(!is.numeric(max_value))
    stop("HSC_PC binning: max_value must be of numeric value")
  if(!is.numeric(bins))
    stop("HSC_PC binning: number of bins must be of integer value")
  if(!is.character(value_field))
    stop("HSC_PC binning: value field name must be of character value")
  structure(list(minv=min_value,
                 maxv=max_value,
                 nbins=bins,
                 vf=value_field), 
            class = c("HSC_PC_Binning","HSC_PC"))
}

classify.HSC_PC_Binning <- function(x, stream, ...) {
  if(!is.data.frame(stream))
    stop("HSC_PC binning: stream must be of data frame type")
  if(!x$vf %in% colnames(stream))
    stop("HSC_PC binning: value field not found in the input stream")
  res <- data.frame(stream,stringsAsFactors=FALSE)
  res[,".clazz"] <- NA
  for(i in 1:nrow(stream)) {
    val <- as.numeric(stream[i,x$vf])
    if(val>x$maxv) {
      bin <- x$nbins
      warning("HSC_PC binning: value out of bounds > max_value... using top bin")
    } else if(val<x$minv) {
      bin <- 1
      warning("HSC_PC binning: value out of bounds < min_value... using bin 1")
    } else {
      n <- (val-x$minv)/(x$maxv-x$minv)
      nvalue <- n*x$nbins
      bin <- ceiling(nvalue)
      if(bin==0) bin <- 1
    }
    res[i,".clazz"] <- as.character(bin)
  }
  return(res)
}