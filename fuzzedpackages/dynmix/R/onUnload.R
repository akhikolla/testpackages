.onUnload <- function(libpath) 
  {
    library.dynam.unload("dynmix",libpath)
  }