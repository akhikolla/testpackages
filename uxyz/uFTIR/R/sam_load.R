#' Sam load
#' 
#' The inverse function for sam_write
#'
#' @param dmdfile a connection to the .bin or .ban file (yes, is true).
#'
#' @return
#' An object of \code{\link[=SAM-class]{SAM}}
#' @export
#'
sam_load <- function(dmdfile){
  
  if(length(grep("\\.bin", dmdfile)) == 1){
    saname <- gsub("sub_", "sam_", dmdfile)
    suname <- gsub("sam_", "sub_", dmdfile)
    cuname <- gsub("sam_", "clu_", dmdfile)
    cuname <- gsub("sub_", "clu_", dmdfile)
    cuname <- gsub("\\.bin", "\\.ban", cuname)
  }
  
  if(length(grep("\\.ban", dmdfile)) == 1){
    cuname <- dmdfile
    dmdfile <- gsub("\\.ban", "\\.bin", dmdfile)
    saname <- gsub("clu_", "sam_", dmdfile)
    suname <- gsub("clu_", "sub_", dmdfile)
  }
  
  raw_sam <- csam_load(saname)
  substances <- csam_load(suname)
  clusters <- csam_load(cuname)
  
  new("SAM",
      raw_sam = raw_sam,
      substances = substances,
      clusters = clusters
  )
}