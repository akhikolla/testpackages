

convts <- function(x,ind=NULL,...)
  {
    if (is.null(ind))
      {
        ind <- as.character(seq(...))
      }
    
    names(x[[1]]) <- ind
    
    for (i in 2:3)
      {
        rownames(x[[i]]) <- ind
      }
    
    if (class(x)=="tvpreg") { names(x[[4]]) <- ind }
    if (class(x)=="mixest") { rownames(x[[4]]) <- ind }
    
    if (class(x)=="mixest")
      {
        if (!any(is.na(x[[5]]))) { names(x[[5]]) <- ind }
        if (!any(is.na(x[[6]]))) { rownames(x[[6]]) <- ind }
      }

    return(x)
  }

