.onLoad <- function(libname, pkgname) {
  op <- options()
  op.rwdataplyr <- list(
    rwdataplyr.wy_month_tol = 6
  )
  toset <- !(names(op.rwdataplyr) %in% names(op))
  if(any(toset)) options(op.rwdataplyr[toset])
  
  invisible()
}

custom_words <- function()
{
  c("asis", "crss", "CRSS", "csv", "rdf", "rdfs", "RiverWare", "riverware", 
    "RiverSMART")
}

.onUnload <- function (libpath) {
  library.dynam.unload("RWDataPlyr", libpath)
}
