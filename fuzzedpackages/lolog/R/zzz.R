



loadModule("lolog", TRUE)
.onLoad <- function(libname, pkgname) {
  .Call(`_lolog_initStats`)
}

.onUnload <- function(libpath) {
}



#library(roxygen2)
#roxygenize('lolog',roclets='rd')
