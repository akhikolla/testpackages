.onUnload <- function(libpath) 
  {
    library.dynam.unload("fDMA",libpath)
  }