HSC_PC_None <- function() {
  structure(list(),class = c("HSC_PC_None","HSC_PC"))
}

classify.HSC_PC_None <- function(x, stream, ...) {
  if(!".clazz" %in% colnames(stream))
    stop("HSC_PC none: .clazz field not found in the input stream")
  res <- data.frame(stream,stringsAsFactors=FALSE)
  return(res)
}