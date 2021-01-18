#' GDALPolygonize in R
#'
#' @description The function calls the C GDALPolygonize routine (using C++). It does overwrite. It uses always 8 point connectedness. It writes on disc a file in "ESRI Shapefile" format.
#' 
#' I wrote a separate package with the function at first. However, I did not want to submit such a small package to CRAN to resolve the 'polygonize' problem. Therefor, I included the source code of that package in uFTIR at a later stage. If you are looking for a package to polygonize a raster file (similar to GDAL polygonize.py) you can use the GDALPolygonize package directly. You can find it in [GitHub](https://github.com/fcorra/GDALPolygonize).
#' 
#' That's also why this function is not exported for you to use it.
#'
#' @param raster name/location of raster file to polygonize. Any driver included in GDAL should work.
#' @param folder name/location of the shape file to write. File extension should not be included.
#' @param layer name of the layer (.shp file) to write. Do not include de file extension.
#' @param field name of the field that will keep the raster values.
#' @param overwrite TRUE/FALSE to overwrite an existing file (folder).
#'
#' @return
#' integer 0
#' @examples
#' NULL
rgdal_polygonize <- function(raster, folder, layer, field, overwrite = TRUE){

  if(!file.exists(raster)){
    warning("\n Source file not found!")
    return(1L)
  }

  # check to overwrite
  dst_layer <- paste(paste(folder, layer, sep = .Platform$file.sep), ".shp", sep = "")
  if(file.exists(dst_layer) & overwrite){
    rmvec <- list.files(folder, pattern = layer, full.names = TRUE)
    file.remove(rmvec)

  } else if(file.exists(dst_layer) &! overwrite){
    warning("\nThe destination file lives already in your system,\nuse overwrite = TRUE if you want to replace it.")
    return(1L)
  }

  gdal_polygonize(raster, folder, layer, field)
  return(0L)
}
