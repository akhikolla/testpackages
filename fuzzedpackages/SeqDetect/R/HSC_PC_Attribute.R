HSC_PC_Attribute <- function(field) {
  if(!is.character(field))
    stop("HSC_PC attribute: class field name must be of character value")
  structure(list(fieldv=field), 
            class = c("HSC_PC_Attribute","HSC_PC"))
}

classify.HSC_PC_Attribute <- function(x, stream, ...) {
  if(!is.data.frame(stream))
    stop("HSC_PC attribute: stream must be of data frame type")
  if(!x$fieldv %in% colnames(stream))
    stop(paste0("HSC_PC attribute: class field ",x$fieldv," not found in the input stream"))
  res <- data.frame(stream,stringsAsFactors=FALSE)
  res[,".clazz"] <- NA
  for(i in 1:nrow(stream))
    res[i,".clazz"] <- as.character(stream[i,x$fieldv])
  return(res)
}