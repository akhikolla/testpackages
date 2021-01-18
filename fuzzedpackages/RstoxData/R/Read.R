##################################################
##################################################
#' Read biotic XML files
#' 
#' This function reads multiple biotic file to a list with a list of tables for each file.
#' 
#' @param FileNames     The paths of the biotic files.
#' @inheritParams general_arguments
#' 
#' @details
#' This function is awesome and does excellent stuff.
#' 
#' @return
#' An object of StoX data type BioticData: A list of a list of data.tables of the different levels of the input biotic files.
#' 
#' @examples
#' x <- 1
#' 
#' @seealso \code{\link[RstoxData]{readXmlFile}}.
#' 
#' @importFrom parallel makeCluster parLapply stopCluster mclapply
#' @export
#' 
ReadBiotic <- function(FileNames, NumberOfCores = 1L) {
	
	# Read BioticData possibly on several cores:
	BioticData <- lapplyOnCores(
		FileNames, 
		FUN = RstoxData::readXmlFile, 
		NumberOfCores = NumberOfCores
	)
	
	# Add names as the file names:
	names(BioticData) <- basename(FileNames)
	
	return(BioticData)
}



##################################################
##################################################
#' Read acoustic XML files
#' 
#' This function reads multiple acoustic file to a list with a list of tables for each file.
#' 
#' @param FileNames     The paths of the acoustic files.
#' @inheritParams general_arguments
#' 
#' @details
#' This function is awesome and does excellent stuff.
#' 
#' @return
#' An object of StoX data type AcousticData: A list of a list of data.tables of the different levels of the input acoustic files.
#' 
#' @examples
#' x <- 1
#' 
#' @seealso \code{\link[RstoxData]{readXmlFile}}.
#' 
#' @importFrom parallel makeCluster parLapply stopCluster mclapply
#' @export
#' 
ReadAcoustic <- function(FileNames, NumberOfCores = 1L) {
	
	# Read AcousticData possibly on several cores:
	AcousticData <- lapplyOnCores(
		FileNames, 
		FUN = RstoxData::readXmlFile, 
		NumberOfCores = NumberOfCores
	)
	
	# Add names as the file names:
	names(AcousticData) <- basename(FileNames)
	
	return(AcousticData)
}

