#' @title Read map file
#' @description Read in the marker map  data. 
#' @param filename contains the name of the map file. The file name needs to be in quotes. If the file is not in the working directory, then the full path 
#' to the file is required. 
#' @param header   a logical value. When \code{TRUE}, the first row of the file contains the column headings. 
#' @param ...   arguments to be passed to read.table such as \code{skip}, \code{sep}. See \code{\link{read.table}} so the list 
#'             of arguments. 
#' @details
#' Association mapping, unlike classical linkage mapping, 
#' does not require a map to find marker-trait associations. 
#' So, reading in a map file is optional. 
#' If a map file is supplied, then the marker names from this file are used when reporting the findings from \code{\link{AM}}. 
#' If a map file is not supplied, then generic names M1, M2, ..., are assigned to the 
#' marker loci where the number refers to the  column number in the marker file. 
#' 
#' A space separated text file with column headings is assumed as the default input. The map file can have three or four 
#' columns. If the map file has three columns, then it is assumed that the three columns are the marker locus names, 
#' the chromosome number, and the map position (in any units). If the map file has four columns as with a 'PLINK map file, 
#' then the columns are assumed to be the marker locus names, the chromosome number, the map position in centimorgans, 
#' and the map position in base pairs. 
#'
#' Missing values are allowed but not in the first column of the file (i.e. the marker labels are not allowed to be missing). 
#' 
#' The order of the marker loci in this file is assumed to be  the same order as the loci in the marker data file.  
#'
#' The first column of the map file is assumed to contain the marker names. 
#'
#' @seealso \code{\link{ReadMarker}} and \code{\link{ReadPheno}}.
#' @return 
#' a data frame is returned of the map data. 
#'
#' @examples
#' # Read in  example map data from ./extdata/
#' 
#' # find the full location of the map data 
#' complete.name <- system.file('extdata', 'map.txt', package='Eagle')
#'   
#' # read in map data 
#' map_obj <- ReadMap(filename=complete.name) 
#'                                
#'# look at first few rows of the map file
#' head(map_obj)
#'
#'
ReadMap  <- function( filename = NULL, header=TRUE, ...)
{
 mapfile <- fullpath(filename)
 error.code <-  check.inputs(file_genotype=filename)
 if(error.code){
    message(" ReadMap has terminated with errors.")
   return(FALSE)
  }
  sep=""
  map <- try(read.table(mapfile, header=header, ...) )

  if (class(map) == "try-error"){
    message(" ReadMap has terminated with errors.")
     return(FALSE)
  }



  if (any(is.na(map[,1]))){
     message("  ")
     message("  ")
     message("  ")
     message(" ERROR: The map file contains missing marker locus labels. ")
     message("        The map file will not be used for analysis. Instead generic marker names will be assigned to the snp.")
     message("   ")
     message("        ReadMap has terminated with errors.")
     message(" ")
     return(FALSE)
  }

message("\n\n Loading map file ... \n\n")
message("                    Summary of Map File  \n")
message("                   ~~~~~~~~~~~~~~~~~~~~~~ \n")
message(" File name:                   ",  mapfile, "\n")
message(" Number of marker loci:       ", nrow(map) , "\n")
message(" Number of columns:           ", ncol(map), "\n")
message(" Number of chromosomes:       ", length(unique(map[[2]])), "\n\n")
message(" First 5 markers of the map file are \n")

if(nrow(map) > 5){
  mat <- as.matrix(map[1:5,])
  for(ii in 1:5){
  message(paste(mat[ii,], collapse="     "))
  }
} else {
  mat <- as.matrix(map)
  for(ii in 1:nrow(map) ){
  message(paste(mat[ii,], collapse="     "))
  }
}


message("The map file has been loaded.")


message("\n\n")

return(map)

}

