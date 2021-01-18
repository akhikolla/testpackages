
# .onLoad <- function(libname, pkgname) {
#   Rcpp::loadModule("chunker_module", TRUE, loadNow = TRUE)
# }

.onUnload <- function (libpath) {
  library.dynam.unload("chunkR", libpath)
}

.onAttach <- function(...) {
  
  vers <- utils::packageDescription("chunkR", fields = "Version")
  
  textstart<- paste("\n", " <-- chunkR --<", "\n\n",
                    "  Version",  vers, "\n",
                    "  GitHub: chunkR_devel()", "\n",
                    "  Comments / suggestions / bug reports: learoser@gmail.com \n")
  
  packageStartupMessage(textstart)
}