summary.PedMap <- function(object, ...){

  .Deprecated("GenomicTools.fileHandler::summary.PedMap", package="GenomicTools", msg="I/O Functions will be collected from now on in a new package GenomicTools.fileHandler")
    
  cat("PedMap Summary\n")
  cat("---------------\n")
  cat("# of imported SNPs       :",nrow(object$map),"\n")
  cat("# of imported samples    :",nrow(object$fam),"\n")
  cat("# of missing sites       :",object$meta$missing,"\n")
  cat("# of monomorphic sites   :",object$meta$mono,"\n")
  cat("# of multiallelic sites  :",object$meta$multiallelic,"\n")
  cat("Used ped file            :",object$meta$pedFile,"\n")
  invisible(object)
} 
