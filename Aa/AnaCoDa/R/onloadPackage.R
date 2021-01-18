NAMESPACE <- environment()
# TODO: Test with exposure and what doesn't need to be exposed to R
Rcpp::loadModule("Test_mod", TRUE)
Rcpp::loadModule("Trace_mod", TRUE)
Rcpp::loadModule("CovarianceMatrix_mod", TRUE)
Rcpp::loadModule("MCMCAlgorithm_mod", TRUE)
Rcpp::loadModule("Model_mod", TRUE)
Rcpp::loadModule("Parameter_mod", TRUE)
Rcpp::loadModule("Genome_mod", TRUE)
Rcpp::loadModule("Gene_mod", TRUE)
Rcpp::loadModule("SequenceSummary_mod", TRUE)

#.onLoad <- function(libname, pkgname){
  #library.dynam("ribModel", pkgname, libname) 
#  invisible()
#} # End of .onLoad().

#.onUnload <- function(libpath){
#  library.dynam.unload("ribModel", libpath)
#  invisible()
#} # End of .onUnload().

#.onAttach <- function(libname, pkgname){
#  invisible()
#} # End of .onAttach().
