#' Write a sam object to disk
#' 
#' A wrapper to the C++ mosaic_sam_write
#'
#' @param tile The tile form which the sam originates
#' @param sam The sam that you want to write to disk
#'
#' @return
#' NULL
#' @export
#'
sam_write <- function(tile, sam){
  
  name <- gsub(".dmd$", "", tile@file)
  
  nsam <- paste(tile@path, "sam_", name, ".bin", sep ="")
  nsub <- paste(tile@path, "sub_", name, ".bin", sep ="")
  nclu <- paste(tile@path, "clu_", name, ".ban", sep ="")
  
  mosaic_sam_write(sam@raw_sam, nsam)
  mosaic_sam_write(sam@substances, nsub)
  mosaic_sam_write(sam@clusters, nclu)

  NULL
}