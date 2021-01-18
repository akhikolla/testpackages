.onAttach = function(libname, pkgname){
  pkg_ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields="Version")
  packageStartupMessage(paste(pkgname, pkg_ver))
}
