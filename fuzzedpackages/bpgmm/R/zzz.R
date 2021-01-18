.onAttch <- function(libname, pkgname) {
  packageStartupMessage(paste(
    "\nThis is bpgmm version",
    utils::packageVersion("bpgmm")
  ))
}
