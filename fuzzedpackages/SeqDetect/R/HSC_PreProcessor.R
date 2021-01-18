HSC_PP <- function(fields, timestamp_field, create_unique_key=FALSE, auto_id=FALSE) {
  structure(list(
    fields = fields,
    ts_field = timestamp_field,
    cuk = create_unique_key,
    aid = auto_id,
    id = 1
  ), class = c("HSC_PP"))
}

preprocess <- function(x, streams,  ...) UseMethod("preprocess")
preprocess.default <- function(x, streams, ...) {
  stop(gettextf("preprocess method is not implemented for class '%s'.",paste(class(x), collapse=", ")))
}

preprocess.HSC_PP <- function(x, streams, ...) {
  if(is.null(streams) || is.null(x$ts_field))
    stop("HSC_PP straight through: not all parameters are defined, some are NULL")
  if(!is.list(streams) || length(streams)==0)
    stop("HSC_PP straight through: streams parameter must be a list and must have at least one element")
  if(!is.vector(x$fields) || length(x$fields)==0)
    stop("HSC_PP straight through: fields parameter must be a vector and must have at least one element")
  if(!is.character(x$ts_field))
    stop("HSC_PP straight through: timestamp field parameter must be a string")
  for(i in 1:length(streams)) {
    if(!is.data.frame(streams[[i]]))
      stop(paste0("HSC_PP straight through: streams[[",i,"]] is not a data frame"))
  }
  res <- data.frame(stringsAsFactors=FALSE)
  if(!x$aid && !x$ts_field %in% x$fields)
    stop(paste0("HSC_PP straight through: fields do not contain the timestamp field '",x$ts_field,"'"))
  for(i in 1:length(streams)) {
    stream <- streams[[i]]
    stream_fields <- colnames(stream)
    if(!x$aid && !x$ts_field %in% stream_fields)
      stop(paste0("HSC_PP straight through: streams[[",i,"]] does not contain the timestamp field '",x$ts_field,"'"))
    stream_filtered <- stream[,base::intersect(x$fields,stream_fields)]
    stream_filtered <- data.frame(stream_filtered,stringsAsFactors=FALSE)
    stream_filtered[,base::setdiff(x$fields,stream_fields)] <- NA
    if(x$cuk) stream_filtered[,".key"] <- 1
    if(x$aid) {
      message("Generate auto identifiers...")
      for(j in 1:nrow(stream)) {
        stream_filtered[j,x$ts_field] <- x$id
        x$id <- x$id + 1
      }
    }
    res <- rbind(res,stream_filtered)
  }
  res <- res[order(res[,x$ts_field]),]
  return(list(obj=x,res=res))
}