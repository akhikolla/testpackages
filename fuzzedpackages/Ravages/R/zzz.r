.onAttach <- function(libname, pkgname) {
  rnorm(1); # force RNG seed initialisation (not done when the RNG is called from C++ code)
  if(r_check_limit_cores())
    Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")
}

.onLoad <- function(libname, pkgname) {
  rnorm(1); # force RNG seed initialisation (not done when the RNG is called from C++ code)
  if(r_check_limit_cores())
    Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")
}


# note : Il se pourrait qu'au lieu de renvoyer "true" certains systèmes renvoient "true " ou quelque chose comme ça
# le plus simple semble etre de tester que c'est défini (à autre chose que "false" au cas où, admettons)...
r_check_limit_cores <- function() { 
  Rcheck <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  (nchar(Rcheck[1]) > 0) & (Rcheck != "false")
}

