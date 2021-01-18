.depthf <- function(fname, funcargs, ...){
  f <- try(match.fun(fname), silent = T)
  if (is.function(f)){
    args = list(...)
    fcnArgs <- names(formals(f))
    fcnArgs <- unlist(fcnArgs, use.names=FALSE)
    keep <- intersect(names(args), fcnArgs)
    unused <- setdiff(names(args), fcnArgs)
    args <- args[keep]
    
    args <- c(args, funcargs)
    res <- do.call(fname, args=args)
    
    if(length(unused)>0)
      warning("Unused by '", fname, "' arguments: ", paste(unused, collapse = ', '))
    
    #res <- f(x, data, ...)
    
    return(res)
  } else {
    warning("There is no depth function ", fname)
  }
}

depthf. <- function(datafA, datafB, notion = c("ABD", "BD", "fd1", "fd2", "hM", "hM2", "HR", "RP1", "RP2"), ...){
  
  if(is.null(notion))
    stop("Parameter 'notion' must be set")
  t <- notion
  try(t <- match.arg(notion), silent = T)
  
  fname = paste0("depthf.", t)
  funcargs = list(datafA = datafA, datafB = datafB)

  return(.depthf(fname, funcargs, ...))
}

#depthf.space. <- function(dataf, cardinalities, notion = c("ABD", "BD", "fd1", "fd2", "hM", "hM2", "HR", "RP1", "RP2"), ...){
#  
#  if(is.null(notion))
#    stop("Parameter 'notion' must be set")
#  t <- notion
#  try(t <- match.arg(notion), silent = T)
#  
#  # try to find a depth
#  fname = paste0("depth.space.", t)
#  funcargs = list(cardinalities = cardinalities, data = data)
#  return(.depth(fname, funcargs, ...))
#}




# d =  depth(data$train, data$train, exact = T)
