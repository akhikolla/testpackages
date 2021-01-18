# Unload the library
.onUnload <- function (libpath) {
  library.dynam.unload("imptree", libpath)
}
