.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("wCorr v", utils::packageDescription("wCorr")$Version, "\n"))
}
