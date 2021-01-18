.onUnload <- function(libpath) {
  library.dynam.unload("benchr", libpath)
}
