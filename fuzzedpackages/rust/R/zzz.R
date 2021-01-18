.onUnload <- function (libpath) {
  library.dynam.unload("rust", libpath)
}
