
.onUnload <- function (libpath) {
  library.dynam.unload("kernelboot", libpath)
}
