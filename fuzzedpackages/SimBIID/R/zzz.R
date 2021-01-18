.onUnload <- function(libpath) {
    library.dynam.unload("SimBIID", libpath)
}

globalVariables(c("Rcpp_object", "Output", "Parameter", "value", "nident",
                  "output", "V1", "uci", "Generation", "p1", "p2", "parnames",
                  "V2", "name", "pair", "output.x", "output.y", "lci", "dist", "med"))
