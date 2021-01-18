
# The general function for redefining StoxData:
RedefineData <- function(
	StoxData, RawData, 
	Redefinition = data.table::data.table(), 
	StoxDataFormat = c("Biotic", "Acoustic"), 
	NumberOfCores = 1L
) {
	
	StoxDataFormat <- match.arg(StoxDataFormat)
	
	# Add the requested variable:
	StoxData <- AddToStoxData(
		StoxData = StoxData, 
		RawData = RawData, 
		VariableNames = Redefinition$ReplaceBy, 
		NumberOfCores = NumberOfCores, 
		StoxDataFormat = StoxDataFormat
	)
	
	# Remove the old:
	lapply(StoxData, replaceAndDelete, VariableReplacement = Redefinition)
	
	return(StoxData)
}

# The general function for defining StoxData translation:
DefineDataTranslation <- function(
	processData, UseProcessData = FALSE, 
	DefinitionMethod = c("Table", "ResourceFile"), 
	Translation = data.table::data.table(), 
	FileName
) {
	
	# Return immediately if UseProcessData = TRUE:
	if(UseProcessData) {
		return(processData)
	}
	
	DefinitionMethod = match.arg(DefinitionMethod)
	
	if(DefinitionMethod == "ResourceFile") {
		# Get the conversion table:
		Translation <- readVariableConversion(
			processData = processData, 
			FileName = FileName, 
			UseProcessData = UseProcessData
		)
	}
	else if(DefinitionMethod != "Table"){
		stop("Invalid DefinitionMethod")
	}
	
	return(Translation)
}

# The general function for translating StoxData:
TranslateData <- function(
	StoxData, 
	TranslationDefinition = c("FunctionParameter", "FunctionInput"), 
	Translation = data.table::data.table(), 
	TranslationProcessData
) {
	
	TranslationDefinition <- match.arg(TranslationDefinition)
	
	if(TranslationDefinition == "FunctionInput") {
		Translation <- TranslationProcessData
	}
	else if(TranslationDefinition != "FunctionParameter"){
		stop("TranslationDefinition must be one of \"FunctionParameter\" and \"FunctionInput\"")
	}
	
	# Apply the translation:
	StoxData <- translateVariables(data = StoxData, Translation = Translation)
	
	return(StoxData)
}

# The general function for converting StoxData:
ConvertData <- function(
	StoxData, 
	ConversionFunction = c("Constant", "Addition", "Scaling", "AdditionAndScaling"), 
	GruopingVariables = c("SpeciesCategory"), 
	Conversion = data.table::data.table()
) {
	
	# Get the ConversionFunction input:
	ConversionFunction <- match.arg(ConversionFunction)
	
	# Check the Conversion for unique grouping variables:
	if(!all(GruopingVariables %in% names(Conversion))) {
		stop("All grouping variables must be present in the Conversion")
	}
	
	# Make a copy to allow for applying functions by reference:
	StoxDataCopy <- data.table::copy(StoxData)
	
	# Merge the data to allow for use of variables from different tables (e.g., length and lengthmeasurement for NMDBiotic 3.0)
	StoxDataCopyMerged <- mergeDataTables(StoxDataCopy, output.only.last = TRUE)
	
	# Merge in the Conversion:
	StoxDataCopyMerged <- mergeByIntersect(StoxDataCopyMerged, Conversion, all.x = TRUE)
	
	# Apply the conversion function for each row of the Conversion:
	ConversionList <- split(Conversion, seq_len(nrow(Conversion)))
	for(Conversion in ConversionList) {
		applyConversionFunction(
			data = StoxDataCopyMerged, 
			ConversionFunction = ConversionFunction, 
			TargetVariable = Conversion$TargetVariable, 
			SourceVariable = Conversion$SourceVariable, 
			RoundOffTo = Conversion$RoundOffTo
		)
	}
	
	# Extract the StoxData from the merged:
	StoxDataCopy <- getStoxDataFromMerged(
		StoxDataMerged = StoxDataCopyMerged, 
		StoxData = StoxDataCopy
	)
	
	
	

	return(StoxDataCopy)
}


getStoxDataFromMerged <- function(StoxDataMerged, StoxData) {
	toExtract <- lapply(StoxData, names)
	lapply(toExtract, function(x) unique(StoxDataMerged[, ..x]))
}


# Function to convert one or more variables of StoxData:
applyConversionFunction <- function(data, ConversionFunction, TargetVariable, SourceVariable, RoundOffTo) {
	
	# Get the conversion function:
	ConversionFunctionName <- paste("ConversionFunction", ConversionFunction, sep = "_")
	do.call(ConversionFunctionName, list(
		data, 
		TargetVariable = TargetVariable, 
		SourceVariable = SourceVariable, 
		RoundOffTo = RoundOffTo
		)
	)
	
	return(data)
}















###  # Function to convert one or more variables of StoxData:
###  applyConversionFunction <- function(data, ConversionFunction) {
###  	
###  	# Get the conversion function:
###  	ConversionFunctionName <- paste("ConversionFunction", ConversionFunction, sep = "_")
###  	
###  	# Get the unique target	and source variable (since the Conversion has been merged with the data):
###  	uniqueTargetAndSource <- getUniqueTargetAndSource(data)
###  	
###  	# Loop through the rows of uniqueTargetAndSource:
###  	rowIndicesOfTargetAndSource <- seq_len(nrow(uniqueTargetAndSource))
###  	for(row in rowIndicesOfTargetAndSource) {
###  		
###  		# Get the current TargetVariable:
###  		TargetVariable <- uniqueTargetAndSource$target[row]
###  		SourceVariable <- uniqueTargetAndSource$source[row]
###  		
###  		# Get the row indices of the data:
###  		#rowIndicesOfData <- 
###  		#	data[[TargetVariable]] == data[[TargetVariable]][row] & 
###  		#	data[[SourceVariable]] == data[[SourceVariable]][row]
###  		#
###  		
###  		# Apply the conversion function:
###  		#data[rowIndicesOfData, eval(TargetVariable) := do.call(ConversionFunctionName, list(.SD, SourceVariable = SourceVariable))]
###  		##data[, eval(TargetVariable) := do.call(ConversionFunctionName, list(.SD, SourceVariable = SourceVariable))]
###  		
###  		do.call(ConversionFunctionName, list(data, TargetVariable = TargetVariable, SourceVariable = SourceVariable))
###  	}
###  	
###  	return(data)
###  }

# The different available conversion functions, reflecting the methods listed as default in ConvertData (parameter ConversionFunction):
ConversionFunction_Constant <- function(data, TargetVariable, SourceVariable, RoundOffTo) {
	# Get the valid rows (those for which the parameters are defined):
	valid <- !is.na(data$Constant)
	# Apply the function:
	data[valid, eval(TargetVariable) := Constant]
	# Round off:
	roundOffValid(data = data, valid = valid, TargetVariable = TargetVariable, RoundOffTo = RoundOffTo)
	#data[valid, eval(TargetVariable) := RoundOff(get(TargetVariable), get(RoundOffTo))]
}

ConversionFunction_Addition <- function(data, TargetVariable, SourceVariable, RoundOffTo) {
	# Get the valid rows (those for which the parameters are defined):
	valid <- !is.na(data$Addition)
	# Apply the function:
	data[valid, eval(TargetVariable) := Addition + get(SourceVariable)]
	# Round off:
	roundOffValid(data = data, valid = valid, TargetVariable = TargetVariable, RoundOffTo = RoundOffTo)
}

ConversionFunction_Scaling <- function(data, TargetVariable, SourceVariable, RoundOffTo) {
	# Get the valid rows (those for which the parameters are defined):
	valid <- !is.na(data$Scaling)
	# Apply the function:
	data[valid, eval(TargetVariable) := Scaling * get(SourceVariable)]
	# Round off:
	roundOffValid(data = data, valid = valid, TargetVariable = TargetVariable, RoundOffTo = RoundOffTo)
}

ConversionFunction_AdditionAndScaling <- function(data, TargetVariable, SourceVariable, RoundOffTo) {
	# Get the valid rows (those for which the parameters are defined):
	valid <- !is.na(data$Addition) & !is.na(data$Scaling)
	# Apply the function:
	data[valid, eval(TargetVariable) := Addition + Scaling * get(SourceVariable)]
	# Round off:
	roundOffValid(data = data, valid = valid, TargetVariable = TargetVariable, RoundOffTo = RoundOffTo)
}


roundOffValid <- function(data, valid, TargetVariable, RoundOffTo) {
	# Round off either to the values of a column or to oa numeric:
	if(!RoundOffTo %in% names(data)) {
		RoundOffToNumeric <- as.numeric(RoundOffTo)
		if(!is.na(RoundOffToNumeric)) {
			# Round off to the RoundOffToNumeric by reference:
			#RoundOffTo <- RoundOffToNumeric
			data[valid, eval(TargetVariable) := roundOff(get(TargetVariable), eval(RoundOffToNumeric))]
		}
		else {
			stop("RoundOffTo must be a character string with either the name of column or a single numeric (coercable to numeric)")
		}
	}
	else {
		# Round off by reference:
		data[valid, eval(TargetVariable) := roundOff(get(TargetVariable), get(RoundOffTo))]	
	}
}

roundOff <- function(x, RoundOffTo) {
	round(x / RoundOffTo) * RoundOffTo
}


# Helper function to get unique conversions:
getUniqueTargetAndSource <- function(data) {
	# Get the defined names of the target and source columns:
	targetAndSourceVariables <- unlist(getRstoxDataDefinitions("targetAndSourceVariables"))
	#targetAndSourceVariablesPresent <- intersect(names(data), targetAndSourceVariables)
	# Uniquify and rename:
	output <- unique(data[, ..targetAndSourceVariables])
	setnames(output, c("target", "source"))
	# Remoev rows with all NAs:
	valid <- rowSums(is.na(output)) < ncol(output)
	output <- output[valid, ]
	
	return(output)
}

# Function for reading a conversion table:
readVariableConversion <- function(processData, FileName, UseProcessData = FALSE) {
	
	# Return immediately if UseProcessData = TRUE:
	if(UseProcessData) {
		return(processData)
	}
	
	conversion <- data.table::fread(FileName)
	
	return(conversion)
}

# Function to convert variables given a conversion table:
translateVariables <- function(data, Translation) {
	
	dataCopy <- data.table::copy(data)
	
	# Currently not defined
	requiredColumns <- getRstoxDataDefinitions("StoxBioticTranslationRequiredColumns")
	if(! all(requiredColumns %in% names(Translation))) {
		stop("The Translation must contain the columns ", paste(requiredColumns, collapse = ", "))
	}
	
	# Split into a list, thus treating only one row at the time. This is probably sloppy coding:
	translationList <- split(Translation, seq_len(nrow(Translation)))
	# Run the conversion for each row of the Translation:
	lapply(translationList, translateVariable, data = dataCopy)
	
	return(dataCopy[])
}

# Function to convert variables given one row of a conversion table:
translateVariable <- function(translationList, data) {
	lapplyToStoxData(data, translateOneTable, translationList = translationList)
}

# Function to apply to all tables of the input data, converting the variables:
translateOneTable <- function(x, translationList) {
	# Check that the table contains the variable to convert:
	if(translationList$VariableName %in% names(x)) {
		# Do nothing if the variable is a key:
		isKeys <- endsWith(translationList$VariableName, "Key")
		if(isKeys) {
			warning("StoX: The variable", translationList$VariableName, " is a key and cannot be modified ")
		}
		else {
			# Convert the class to the class of the existing value in the table:
			translationList <- convertClassToExisting(translationList, x)
			# Replace by the new value:
			x[, c(translationList$VariableName) := replace(
				x = get(translationList$VariableName), 
				list = get(translationList$VariableName) %in% translationList$Value, 
				values = translationList$NewValue)]
		}
	}
}

# Function to convert the class of the Value and NewValue of a translationList to the class of the existing value:
convertClassToExisting <- function(translationList, x) {
	# Convert the NewValue to the class of the existing value:
	existingClass <- class(x[[translationList$VariableName]])[1]
	newClass <- class(translationList$Value)[1]
	if(!identical(existingClass, newClass)) {
		class(translationList$Value) <- existingClass
		class(translationList$NewValue) <- existingClass
	}
	return(translationList)
}

# Function to apply a function to StoX data, which may be a list of tables or a list of lists of tables:
lapplyToStoxData <- function(x, fun, ...) {
	# Check the depth of the list, either with tables at the top level og with lists of tables:
	if(is.list(x[[1]]) && !data.table::is.data.table(x[[1]])) {
		lapply(x, function(y) lapply(y, fun, ...))
	}
	else {
		lapply(x, fun, ...)
	}
}

# Function to replace the existing column by the new, as stored in the VariableReplacement:
replaceAndDelete <- function(table, VariableReplacement) {
	present <- which(VariableReplacement$VariableName %in% names(table))
	if(any(present)) {
		# Delete the present column:
		table[, (VariableReplacement[present, VariableName]) := NULL]
		# ... and then rename the new to the old name:
		setnames(table, VariableReplacement[present, ReplaceBy], VariableReplacement[present, VariableName])
	}
}





##################################################
#' Redefine StoxBioticData variables by data from BioticData
#' 
#' This function redefines one or more columns of \code{\link{StoxBioticData}} by columns of \code{\link{BioticData}}.
#' 
#' @param StoxBioticData An input of \link{ModelData} object
#' @param BioticData An input of \link{ModelData} object
#' @param Redefinition A table of the columns "VariableName", representing the variable to redefine; and "RedefineBy", representing the variable from BioticData to replace by. 
#' 
#' @return
#' A \code{\link{StoxBioticData}} object.
#' 
#' @export
#' 
RedefineStoxBiotic <- function(
	StoxBioticData, BioticData, 
	Redefinition = data.table::data.table()
) {
	# Redefine StoxBioticData:
	RedefineData(
		StoxData = StoxBioticData, RawData = BioticData, 
		Redefinition = Redefinition, 
		StoxDataFormat = "Biotic"
	)
}


##################################################
#' Define StoxBioticData variables translation
#' 
#' This function defines the translation table used as input to \code{\link{TranslateStoxBiotic}} to translate values of one or more columns of \code{\link{StoxBioticData}} to new values given by a table or read from a CSV file.
#' 
#' @inheritParams general_arguments
#' @param DefinitionMethod  Character: A string naming the method to use, one of "Table" for defining the \code{Translation}, and "ResourceFile" for reading the table from the file given by \code{FileName}.
#' @param Translation A table of the columns "VariableName", representing the variable to translate; "Value", giving the values to translate; and "NewValue", giving the values to translate to.
#' @param FileName The csv file holding a table with the three variables listed for \code{Translation}.
#' 
#' @return
#' A \code{\link{StoxBioticTranslation}} object.
#' 
#' @export
#' 
DefineStoxBioticTranslation <- function(
	processData, UseProcessData = FALSE, 
	DefinitionMethod = c("Table", "ResourceFile"), 
	Translation = data.table::data.table(), 
	FileName
) {
	# Define translation for StoxBioticData:
	DefineDataTranslation(
		processData = processData, UseProcessData = UseProcessData, 
		DefinitionMethod = DefinitionMethod, 
		Translation = Translation, 
		FileName = FileName
	)
}

##################################################
#' Translate StoxBioticData
#' 
#' This function translates one or more columns of \code{\link{StoxBioticData}} to new values given by the table \code{Translation} or by the input \code{StoxBioticTranslation}.
#' 
#' @param StoxBioticData An input of \link{ModelData} object
#' @param TranslationDefinition  Character: A string naming the method to use for the translation, one of "FunctionParameter" for defining the \code{Translation}, and "FunctionInput" for using the table produced by the process given by the function input \code{StoxBioticTranslation}.
#' @param Translation A table of the columns "VariableName", representing the variable to translate; "Value", giving the values to translate; and "NewValue", giving the values to translate to.
#' @param StoxBioticTranslation The process from which to get the \code{\link{StoxBioticTranslation}} definition.
#' 
#' @return
#' A \code{\link{StoxBioticData}} object.
#' 
#' @export
#' 
TranslateStoxBiotic <- function(
	StoxBioticData, 
	TranslationDefinition = c("FunctionParameter", "FunctionInput"), 
	Translation = data.table::data.table(), 
	StoxBioticTranslation
) {
	# Translate StoxBioticData:
	TranslateData(
		StoxData = StoxBioticData, 
		TranslationDefinition = TranslationDefinition, 
		Translation = Translation, 
		TranslationProcessData = StoxBioticTranslation
	)
}


##################################################
#' Convert StoxBioticData
#' 
#' This function converts one or more columns of \code{\link{StoxBioticData}} by the function given by \code{ConversionFunction}.
#' 
#' @param StoxBioticData An input of \link{ModelData} object
#' @param ConversionFunction  Character: The function to convert by, one of "Constant", for replacing the specified columns by a constant value; "Addition", for adding to the columns; "Scaling", for multiplying by a factor; and "AdditionAndScaling", for both adding and multiplying.
#' @param GruopingVariables A vector of variables to specify in the \code{Conversion}. The parameters specified in the table are valid for the combination of the \code{GruopingVariables} in the data.
#' @param Conversion A table of the \code{GruopingVariables} and the columns "TargetVariable", "SourceVariable" and the parameters of the \code{ConversionFunction} (see details).
#' 
#' The parameters of the \code{ConversionFunction} are "Constant" for ConversionFunction "Constant", "Addition" for ConversionFunction"Addition", "Scaling" for ConversionFunction "Scaling", and "Addition" and "Scaling" for ConversionFunction "AdditionAndScaling".
#' 
#' @return
#' A \code{\link{StoxBioticData}} object.
#' 
#' @export
#' 
ConvertStoxBiotic <- function(
	StoxBioticData, 
	ConversionFunction = c("Constant", "Addition", "Scaling", "AdditionAndScaling"), 
	GruopingVariables = c("SpeciesCategory"), 
	Conversion = data.table::data.table()
) {
	# Convert StoxBioticData:
	ConvertData(
		StoxData = StoxBioticData, 
		ConversionFunction = ConversionFunction,
		GruopingVariables = GruopingVariables,
		Conversion = Conversion
	)
}









