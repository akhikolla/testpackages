#' Function specification for inclusion in StoX projects
#' @export
stoxFunctionAttributes <- list(
	# Read input biotic data:
	ReadBiotic = list(
		functionType = "modelData", 
		functionCategory = "baseline", 
		functionOutputDataType = "BioticData", 
		#functionParameterType = list(FileNames = "character"), 
		functionParameterFormat = list(FileNames = "filePaths"), 
		functionArgumentHierarchy = list()
	), 
	
	# Read input biotic data:
	ReadAcoustic = list(
		functionType = "modelData", 
		functionCategory = "baseline", 
		functionOutputDataType = "AcousticData", 
		#functionParameterType = list(FileNames = "character"), 
		functionParameterFormat = list(FileNames = "filePaths"), 
		functionArgumentHierarchy = list()
	), 
	
	# Convert AcousticData to StoxAcousticData:
	StoxAcoustic = list(
		functionType = "modelData", 
		functionCategory = "baseline", 
		functionOutputDataType = "StoxAcousticData", 
		functionArgumentHierarchy = list()
	),
	
	# Convert BioticData to StoxBioticData:
	StoxBiotic = list(
		functionType = "modelData", 
		functionCategory = "baseline", 
		functionOutputDataType = "StoxBioticData", 
		functionArgumentHierarchy = list()
	),
	
	# Convert LandingData to StoxLandingData:
	StoxLanding = list(
	 functionType = "modelData",
	 functionCategory = "baseline",
	 functionOutputDataType = "StoxLandingData",
	 #functionParameterType = list(StoxLanding = "character",
	 #                             appendColumns = "character",
	 #                             appendColumnsNames = "character"),
	 functionArgumentHierarchy = list()
	), 
	
	# Convert BioticData to StoxBioticData:
	FilterBiotic = list(
	    functionType = "modelData", 
	    functionCategory = "baseline", 
	    functionOutputDataType = "BioticData", 
	    functionParameterFormat = list(FilterExpression = "filterExpressionList")
	),
	
	# Convert BioticData to StoxBioticData:
	FilterAcoustic = list(
	    functionType = "modelData", 
	    functionCategory = "baseline", 
	    functionOutputDataType = "AcousticData", 
	    functionParameterFormat = list(FilterExpression = "filterExpressionList")
	),
	
	# Convert BioticData to StoxBioticData:
	FilterStoxBiotic = list(
	    functionType = "modelData", 
	    functionCategory = "baseline", 
	    functionOutputDataType = "StoxBioticData", 
	    functionParameterFormat = list(FilterExpression = "filterExpressionList")
	),
	
	# Convert BioticData to StoxBioticData:
	FilterStoxAcoustic = list(
	    functionType = "modelData", 
	    functionCategory = "baseline", 
	    functionOutputDataType = "StoxAcousticData", 
	    functionParameterFormat = list(FilterExpression = "filterExpressionList")
	),
	
	# Convert AcousticData to StoxAcousticData:
	MergeStoxAcoustic = list(
	    functionType = "modelData", 
	    functionCategory = "baseline", 
	    functionOutputDataType = "MergeStoxAcousticData", 
	    functionArgumentHierarchy = list()
	),
	
	# Convert BioticData to StoxBioticData:
	MergeStoxBiotic = list(
		functionType = "modelData", 
		functionCategory = "baseline", 
		functionOutputDataType = "MergeStoxBioticData", 
		functionArgumentHierarchy = list()
	),
	
	
	##### Define and Convert variables: #####
	# StoxBiotic:
	RedefineStoxBiotic = list(
		functionType = "modelData", 
		functionCategory = "baseline", 
		functionOutputDataType = "StoxBioticData", 
		functionParameterFormat = list(
			Redefinition = "redefinitionTable"
		)
	),
	
	DefineStoxBioticTranslation = list(
		functionType = "processData", 
		functionCategory = "baseline", 
		functionOutputDataType = "StoxBioticTranslation", 
		functionParameterFormat = list(
			Translation = "translationTable", 
			FileName = "filePath"
		), 
		functionArgumentHierarchy = list(
			DefinitionMethod = list(
				UseProcessData = FALSE
			), 
			# These two are joined with AND, and must both be fulfilled:
			Translation = list(
				DefinitionMethod = "Table", 
				UseProcessData = FALSE
			), 
			# These two are joined with AND, and must both be fulfilled:
			FileName = list(
				DefinitionMethod = "ResourceFile", 
				UseProcessData = FALSE
			)
		)
	),
	
	TranslateStoxBiotic = list(
		functionType = "modelData", 
		functionCategory = "baseline", 
		functionOutputDataType = "StoxBioticData", 
		functionParameterFormat = list(
			Translation = "translationTable"
		), 
		functionArgumentHierarchy = list(
			Translation = list(
				TranslationDefinition = "FunctionParameter"
			), 
			StoxBioticTranslation = list(
				TranslationDefinition = "FunctionInput"
			)
		)
	),
	
	ConvertStoxBiotic = list(
		functionType = "modelData", 
		functionCategory = "baseline", 
		functionOutputDataType = "StoxBioticData", 
		functionParameterFormat = list(
			Conversion = "conversionTable"
		)
	),
	
	AddToStoxBiotic = list(
		functionType = "modelData", 
		functionCategory = "baseline", 
		functionOutputDataType = "StoxBioticData", 
		functionParameterFormat = list(
			VariableNames = "variableNames_AddToStoxBiotic"
		)
	), 
	
	
	WriteICESAcoustic = list(
		functionType = "modelData", 
		functionCategory = "baseline", 
		functionOutputDataType = "ICESAcousticData"
	), 
	WriteICESBiotic = list(
		functionType = "modelData", 
		functionCategory = "baseline", 
		functionOutputDataType = "ICESBioticData"
	), 
	WriteICESDatras = list(
		functionType = "modelData", 
		functionCategory = "baseline", 
		functionOutputDataType = "ICESDatrasData"
	)
)

#' Define the process property formats:
#' 
#' @export
#' 
processPropertyFormats <- list(
	filePath = list(
		class = "single", 
		title = "The path to a single file"
	), 
	filePaths = list(
		class = "vector", 
		title = "The path to one or more files", 
		variableTypes <- "character"
	), 
	filterExpressionList = list(
		class = "list", 
		title = "A list of filter expressions, one for each table to filter on"
	), 
	
	redefinitionTable = list(
		class = "table", 
		title = "Define columns of StoX data by columns of the raw data", 
		columnNames = c(
			"VariableName", 
			"ReplaceBy"
		), 
		variableTypes = c(
			"character",
			"character"
		)
	), 
	translationTable = list(
		class = "table", 
		title = "Translate columns of StoX data", 
		columnNames = c(
			"VariableName", 
			"Value", 
			"NewValue"
		), 
		variableTypes = c(
			"character",
			"character",
			"character"
		)
	), 
	
	conversionTable = list(
		class = "table", 
		title = function(ConversionFunction = c("Constant", "Addition", "Scaling", "AdditionAndScaling")) {
			ConversionFunction <- match.arg(ConversionFunction)
			
			if(identical(ConversionFunction, "Constant")) {
				title <- "Replace variables of StoX data by a constant \"Constant\""
			}
			else if(identical(ConversionFunction, "Addition")) {
				title <- "Add the value \"Addition\" to variables of StoX data"
			}
			else if(identical(ConversionFunction, "Scaling")) {
				title <- "Multiply variables of StoX data by the value \"Scaling\""
			}
			else if(identical(ConversionFunction, "AdditionAndScaling")) {
				title <- "Multiply variables of StoX data by the value \"Scaling\" and add the value \"Addition\""
			}
			else {
				stop("Wrong ConversionFunction.")
			}
			
			return(title)
		}, 
		columnNames = function(ConversionFunction = c("Constant", "Addition", "Scaling", "AdditionAndScaling"), GruopingVariables = NULL) {
			ConversionFunction <- match.arg(ConversionFunction)
			
			if(identical(ConversionFunction, "Constant")) {
				parameters <- c("Constant")
			}
			else if(identical(ConversionFunction, "Addition")) {
				parameters <- "Addition"
			}
			else if(identical(ConversionFunction, "Scaling")) {
				parameters <- c("Scaling")
			}
			else if(identical(ConversionFunction, "AdditionAndScaling")) {
				parameters <- c("Addition", "Scaling")
			}
			else {
				stop("Wrong ConversionFunction.")
			}
			
			columnNames <- c(
				GruopingVariables, 
				c("TargetVariable", "SourceVariable"), 
				parameters, 
				"RoundOffTo"
			)
			
			return(columnNames)
		}, 
		variableTypes = function(ConversionFunction = c("Constant", "Addition", "Scaling", "AdditionAndScaling"), GruopingVariables = NULL) {
			ConversionFunction <- match.arg(ConversionFunction)
			
			if(identical(ConversionFunction, "Constant")) {
				types <- "double"
			}
			else if(identical(ConversionFunction, "Addition")) {
				types <- "double"
			}
			else if(identical(ConversionFunction, "Scaling")) {
				types <- "double"
			}
			else if(identical(ConversionFunction, "AdditionAndScaling")) {
				types <- c("double", "double")
			}
			else {
				stop("Wrong ConversionFunction.")
			}
			
			variableTypes <- c(
				rep("character", length(GruopingVariables)), 
				c("character", "character"), 
				types, 
				"character"
			)
			
			return(variableTypes)
		}
	), 
	
	variableNames_AddToStoxBiotic = list(
		class = "vector", 
		title = "One or more variables to add to the StoxBioticData from BioticData", 
		possibleValues = function(BioticData) {
			sort(unique(unlist(lapply(BioticData, function(x) lapply(x, names)))))
		}, 
		variableTypes <- "character"
	)
)





