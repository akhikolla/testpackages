.onLoad <- function(libname, pkgname){
  for (cc in cnsts)
    assign(cc, get_cnst(cc), envir = parent.env(environment()))
}
