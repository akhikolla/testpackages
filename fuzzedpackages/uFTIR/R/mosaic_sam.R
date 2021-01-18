#' Mosaic SAM
#' 
#' @description The function performs the spectral angle mapper algorithm chunk by chunk for the mosaics files. It uses the parallel package, and by default all cores -1.
#' 
#' The function uses \code{\link{mosaic_chunk}} to cast a mosaic tile as a \code{\link[=Tile-class]{Tile}} object. Then performs a \code{\link{tile_base_corr}}, \code{\link[=TileRead.wavealign]{wavealing}}, and \code{\link{tile_sam}} to finaly write binary files that hold the Spectral Angle Mapper results in the destination folder (path slot of the \code{\link[=SpectralInfo-class]{SpectralInfo}} object). The files can be loaded back to R using \code{link{mosaic_compose}}.
#' 
#' The function is using the parallel package, as each tile_sam takes a while to complete.
#'
#' @param info \code{\link[=SpectralInfo-class]{SpectralInfo}} object.
#' @param sref \code{\link[=SpectralReference-class]{SpectralReference}} object.
#' @param derivative whether to apply the first (1) or second (2) derivative before sam. Default NULL.
#' @param base_corr TRUE/FALSE should \code{\link{tile_base_corr}} be call before processing each chunk?
#' @param FUN A function to be passed to \code{\link{preprocess}}.It is always applied as if 'data' were a \code{\link[=SpectralPack-class]{SpectralPack}} object.
#' @param n_cores The number of cores to parallelize the task. NULL means all cores -1.
#'
#' @return
#' TRUE
#' 
#' @export
#' @seealso 
#' For a single tile application see \code{\link{tile_sam}}.
#' @examples
#' x <- mosaic_info(base::system.file("extdata/mosaic.dmt", package = "uFTIR"))
#' mosaic_sam(x, primpke, n_cores = 1)
mosaic_sam <- function(info, sref, derivative = NULL, base_corr = TRUE, FUN = NULL, n_cores = NULL){
  
  fname <- gsub(info@path, "", info@file)
  fname <- gsub(".dmt$", "", fname)
  fname <- gsub("/", "", fname)
  
  dmdfiles <- list.files(pattern = fname, path = info@path)
  dmdfiles <- dmdfiles[grep(".dmd$", dmdfiles)]
  
  if(is.null(n_cores)){
    n_cores <- detectCores() - 1  
  }
  cl <- makeCluster(n_cores)
  
  clusterExport(cl, c("info", "sref", "derivative", "base_corr", "FUN"), envir=environment())
  
  parLapply(cl, dmdfiles,
            function(x){
              out <- mosaic_chunk(info=info, dmdfile=x)
              if(base_corr){
                out <- tile_base_corr(out)
              }
              out <- wavealign(out, sref)
              if(is.function(FUN)){
                out <- preprocess(out, FUN)
              }
              out <- tile_sam(out, derivative)
              
              name_out <- unlist(strsplit(x, "_"))
              name_out <- name_out[grep('[0-9]{4}', name_out)]
              name_out <- gsub('\\.[[:alnum:]]*$', '', name_out)
              name_out <- paste(name_out[1], name_out[2], sep = "_")
              
              sam_out_name <- paste("sam_", name_out, ".bin", sep = "")
              sub_out_name <- paste("sub_", name_out, ".bin", sep = "")
              
              sam_out_name <- paste(info@path, sam_out_name, sep = .Platform$file.sep)
              sub_out_name <- paste(info@path, sub_out_name, sep = .Platform$file.sep)
              
              mosaic_sam_write(out@raw_sam, sam_out_name)
              mosaic_sam_write(out@substances, sub_out_name)
            }
  )
  stopCluster(cl)
  
  TRUE
}