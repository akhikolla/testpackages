# http://r-pkgs.had.co.nz/src.html#c-best-practices

.onUnload <- function (libpath) {
  library.dynam.unload("FDX", libpath)
}