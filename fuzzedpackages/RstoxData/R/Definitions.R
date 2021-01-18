##################################################
##################################################
#' Definitions stored in the RstoxData environment
#' 
#' This function declares the RstoxData environment and writes vital definitions to it.
#' 
#' @return
#' A list of definitions.
#' 
#' @noRd
#' @seealso Use \code{\link{getRstoxDataDefinitions}} to get the definitions.
#' 
initiateRstoxData <- function(){
	
	# Define the number of digits (12) and the number of significant digits (6, used if values are very low) used by the Rstox packages:
	digits <- 12
	signifDigits <- 6
	
	# Define the time format used by Stox formats:
	StoxDateTimeFormat <- "%Y-%m-%d %H:%M:%OS"
	StoxTimeZone <- "UTC"
	
	# Get the path to the extdata folder:
	fpath <- system.file("extdata", package = "RstoxData")
	
	# Define formats that contain non-unique variables, i.e., columns with the same name in different tables:
	nonUniqueFormats <- c(
		"nmdbioticv1", 
		"nmdbioticv1.1", 
		"nmdbioticv1.2", 
		"nmdbioticv1.3", 
		"nmdbioticv1.4"
	)
	
	# StoxBioticKeys: 
	StoxBioticKeys <- c(
		"CruiseKey", 
		"StationKey", 
		"HaulKey", 
		"SpeciesCategoryKey", 
		"SampleKey", 
		"IndividualKey", 
		"SubIndividualKey"
	)
	# StoxBioticKeys: 
	StoxAcousticKeys <- c(
		"CruiseKey", 
		"LogKey", 
		"BeamKey", 
		"AcousticCategoryKey", 
		"ChannelReferenceKey", 
		"NASCKey"
	)
	
	targetAndSourceVariables <- list(
		target = "TargetVariable", 
		source = "SourceVariable"
	)
	
	# Define the columns required for VariableConversionTable:
	StoxBioticTranslationRequiredColumns <- c("VariableName", "Value", "NewValue")
	
	#### Assign to RstoxDataEnv and return the definitions: ####
	definitionsNames <- ls()
	definitions <- lapply(definitionsNames, get, pos = environment())
	names(definitions) <- definitionsNames
	
	#### Create the RstoxDataEnv environment, holding definitions on folder structure and all the projects. This environment cna be accesses using RstoxData:::RstoxDataEnv: ####
	assign("RstoxDataEnv", new.env(), parent.env(environment()))
	assign("definitions", definitions, envir=get("RstoxDataEnv"))
	
	#### Return the definitions: ####
	definitions
}


##################################################
##################################################
#' Get RstoxData definitions
#' 
#' This function gets vital definitions from the RstoxData environment.
#' 
#' @param name  An optional string vector denoting which definitions to extract.
#' @param ...   values overriding the values of definitions.
#' 
#' @return
#' A list of definitions.
#' 
#' @examples
#' getRstoxDataDefinitions()
#' 
#' @export
#' 
getRstoxDataDefinitions <- function(name = NULL, ...) {
	
	# Save the optional inputs for overriding the output:
	l <- list(...)
	
	# Get all or a subset of the definitions:
	definitions <- get("RstoxDataEnv")$definitions
	if(length(name)){
		definitions <- definitions[[name]]
	}
	
	l <- l[names(l) %in% names(definitions)]
	if(length(l)){
		definitions <- utils::modifyList(definitions, l)
	}
	
	definitions
}
