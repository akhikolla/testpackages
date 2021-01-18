.onLoad <- function(libname, pkgname) {
  op <- options()
  op.PUlasso <-list(PUlasso.skip_fitting=FALSE)
  toset <- !(names(op.PUlasso) %in% names(op))
  if(any(toset)) options(op.PUlasso[toset])
  invisible()
}

.onUnload <- function (libpath) {
  library.dynam.unload("PUlasso", libpath)
}



