.onAttach <- function(...){
  ## Retrieve Year Information
  date <- date()
  x <- regexpr("[0-9]{4}", date)
  this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
  
  # Retrieve Current Version
  this.version = packageVersion("RiemBase")
  
  ## Print on Screen
  packageStartupMessage("-------------------------------------------------------------------")
  packageStartupMessage("** RiemBase ")
  packageStartupMessage("**  - Functions and C++ Header Files for Computation on Manifolds.")
  packageStartupMessage("** ")
  packageStartupMessage("** Version    : ",this.version,"      (",this.year,")",sep="")
  packageStartupMessage("** Maintainer : Kisung You (kyoustat@gmail.com)")
  packageStartupMessage("**")
  packageStartupMessage("** Please share any bugs or suggestions with the maintainer.")
  packageStartupMessage("-------------------------------------------------------------------")
}

.onUnload <- function(libpath) {
  library.dynam.unload("RiemBase", libpath)
}