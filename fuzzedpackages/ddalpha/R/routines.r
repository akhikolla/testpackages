is.numeric_data.frame <- function(x){
  if (is.data.frame(x) && all(sapply(x,base::is.numeric)))
    return (T)
  return (F)
}

is.numeric <- function(x){
  if (base::is.numeric(x))
    return (T)
  if (is.data.frame(x) && all(sapply(x,base::is.numeric)))
    return (T)
  return (F)
}

resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}

is.sorted <- function(x)  {
  return(!is.unsorted(x))
}

tryCatchCapture <- function(expr, warn = T, err = T) {
  val <- NULL
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, w$message)
    invokeRestart("muffleWarning")
  }
  myError <- NULL
  eHandler <- function(e) {
    myError <<- e$message
    NULL
  }
  if(warn && err){
    val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
    return(list(value = val, warnings = myWarnings, error=myError))
  }
  if(warn){
    val <- tryCatch(withCallingHandlers(expr, warning = wHandler))
    return(list(value = val, warnings = myWarnings))
  }
  if(err){
    val <- tryCatch(expr, error = eHandler)
    return(list(value = val, error=myError))
  }
  val <- expr
  return(list(value = val))
}
