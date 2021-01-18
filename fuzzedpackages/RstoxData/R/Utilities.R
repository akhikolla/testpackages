#' Merge list of data tables recursively
#'
#' @param data A list of data tables.
#' @param tableNames A character vector holding the names of the tables to merge.
#' @param output.only.last Only returns last merged table.
#' @param ... Extra parameters that will be passed into \code{\link[data.table]{merge}}.
#'
#' @return A merged data table.
#'
#' @export
#' 
mergeDataTables <- function(data, tableNames = NULL, output.only.last = FALSE, ...) {
	
	# Better use xsdObjects for getting header vars from XML data for merging
	## Get data type:
	plen <- NULL
	if(!is.null(data[["metadata"]])) {
		if(!exists("xsdObjects"))
			xsdObjects <- RstoxData::xsdObjects
		datatype <- unlist(data[["metadata"]][1, "useXsd"])
		plen <- xsdObjects[[paste0(datatype, ".xsd")]]$prefixLens
	}

    # Merge all tables by default:
	if(length(tableNames) == 0) {
		tableNames <- names(data)
	}
	# No merging if only one table given in 'tableNames':
	else if(length(tableNames) == 1)  {
		return(data)
	}
	
	# Make sure tableNames are ordered as in the data:
	dataNames <- names(data)
	tableNames <- dataNames[match(tableNames, dataNames)]
	tableNames <- tableNames[!is.na(tableNames)]
	
	# Merge
	for(ii in 2:length(tableNames)) {
		curr <- tableNames[ii]
		prev <- tableNames[(ii-1)]

		if(!is.null(plen) && !is.na(plen[prev]))
			vars <- names(data[[curr]])[1:plen[prev]]
		else
			vars <- intersect(names(data[[curr]]), names(data[[prev]]))

		# There can be duplicate names between two tables, see that we fix them by adding appropriate suffix before merging
		duplicates <- intersect(setdiff(names(data[[prev]]), vars), setdiff(names(data[[curr]]), vars))
		for(ddpl in duplicates) {
			message(paste("Duplicate columns in merging", prev, "and", curr,  ": ", ddpl, "->", paste0(ddpl, ".", curr)))
			setnames(data[[curr]], ddpl, paste0(ddpl, ".", curr))
		}
		
		data[[curr]] <- merge(data[[prev]], data[[curr]], by=vars, suffixes = suffixes, ...)
	}

	# If tableNamestableNames == "last", return the last table:
	if(output.only.last) {
		data <- data[[utils::tail(tableNames, 1)]]
	}
	
	return(data)
}

#' Merge two data tables by the intersect of the names
#'
#' @param x,y Data tables of class \code{\link[data.table]{data.table}}.
#' @param ... Various overrides.
#' @param msg Verbose message switch, default to \code{FALSE}.
#'
#' @return A merged data table.
#'
#' @export
#' 
mergeByIntersect <- function(x, y, ..., msg = FALSE) {
	# Cascading merge if a list of tables is given:
	if(length(x) > 1  && 
	   is.list(x)  &&  
	   !data.table::is.data.table(x)  && 
	   data.table::is.data.table(x[[1]])) {
		for(ind in seq(2, length(x))) {
			x[[ind]] <- mergeByIntersect(x[[ind - 1]], x[[ind]], ..., msg = msg)
		}
		output <- x[[ind]]
	}
	else {
		by <- intersect(names(x), names(y))
		if(msg) {
			message("Merging by ", paste(by, collapse = ", "))
		}
		if(length(by)) {
			output <- merge(x, y, by = by, ...)
		}
		else {
			stop("No intersect between the names of x and y")
		}
	}
	
	return(output)
}


#' Merge two data tables by StoX keys
#'
#' @param x,y Data tables of class \code{\link[data.table]{data.table}}.
#' @param StoxDataType Input data type. Text string of \code{StoxBiotic}
#' 			or \code{StoxAcoustic}.
#' @param toMergeFromY Specify key columns from \code{y}. \code{NULL} means
#' 			all similarly named columns from \code{x} and \code{y} will be
#' 			merged. Default to \code{NULL}.
#' @param replace Whether to replace the variables in the target.
#' 			Default to \code{FALSE}.
#' @param ... Extra parameters that will be passed into \code{\link[data.table]{merge}}.
#'
#' @return A merged data table.
#'
#' @export
#' 
mergeByStoxKeys <- function(x, y, StoxDataType, toMergeFromY = NULL, replace = FALSE, ...) {
	# Get the keys:
	#keys_x <- getKeys(x)
	#keys_y <- getKeys(y)
	#keys <- intersect(keys_x, keys_y)
    keys <- Reduce(intersect, 
        list(
            names(x), 
            names(y), 
            getStoxKeys(StoxDataType = StoxDataType)
        )
    )
	
	# Define the columns to merge:
	if(!length(toMergeFromY)) {
		toMergeFromY <- names(y)
	}
	# Make sure the toMergeFromY are present in y:
    toMergeFromY <- intersect(names(y), toMergeFromY)
    # Exclcude the keys:
	toMergeFromY <- setdiff(toMergeFromY, getStoxKeys(StoxDataType = StoxDataType))

	#  Replace the variable in the target:
	if(replace) {
		keep <- setdiff(names(x), toMergeFromY)
		x <- x[, ..keep]
	}
	
	# If there are any left, extract the keys and toMergeFromY:
	if(length(toMergeFromY)) {
		y <- y[, c(keys, toMergeFromY), with = FALSE]
		# Then merge:
		merge(x, y, by = keys, ...)
	}
	else {
		x
	}
}

#getKeys <- function(x, keystring = "Key", ignore.case = FALSE) {
#	namesx <- names(x)
#	namesx[endsWith(if(ignore.case) tolower(namesx) else namesx, if(ignore.case) tolower(keystring) else keystring#)]
#}

#' Get the keys of a StoX format
#' 
#' @param StoxDataType The name of the StoX format (only StoxBiotic implemented yet).
#' @param level The name of the level/table to get keys for.
#' @param keys.out Specification of what to return. One of "all", to return all keys of the level; "only.present", to return only the key of the \code{level}; and "all.but.present", to return all keys except the present key.
#'
#' @importFrom data.table key
#' @export
#' 
getStoxKeys <- function(StoxDataType = c("StoxBiotic", "StoxAcoustic"), level = NULL, keys.out = c("all", "only.present", "all.but.present")) {
	StoxDataType <- match.arg(StoxDataType)
	if(StoxDataType == "StoxBiotic") {
	    keys <- stoxBioticObject$convertTable[key == "Y", c("variable", "level")]
	    keys <- split(keys, by = "level")
	    keys <- lapply(keys, "[[", "variable")
	}
	else if(StoxDataType == "StoxAcoustic") {
		stop("Not yet implemented")
	}
	
	if(length(level)) {
		keys <- keys[[level]]
	}
	else {
	    keys <- unique(unlist(keys))
	    return(keys)
	}
	
	keys.out <- match.arg(keys.out)
	if(keys.out == "only.present") {
		keys <- utils::tail(keys, 1)
	}
	else if(keys.out == "all.but.present") {
		keys <- keys[-length(keys)]
	}
	
	return(keys)
}




# Detect OS
get_os <- function() {
	if (.Platform$OS.type == "windows") {
		"win"
	} else if (Sys.info()["sysname"] == "Darwin") {
		"mac"
	} else if (.Platform$OS.type == "unix") {
		"unix"
	} else {
		stop("Unknown OS")
	}
}

# Pick a suitable number of cores
#' @importFrom parallel detectCores
getCores <- function() {
	cores <- as.integer(getOption("mc.cores"))
	if (length(cores) == 0 || is.na(cores)) {
		cores <- parallel::detectCores()
		if (is.na(cores)) {
			return(1)
		} else {
			# Don't use too many cores in autodetect
			if (cores > 4)
				return(4)
			else
				return(cores)
		}
	} else {
		return(cores)
	}
}

#' Run a function on all elements of x on one or more cores
#'
#' @param x An object to apply \code{FUN} to.
#' @param FUN The function to apply.
#' @inheritParams general_arguments
#' @param ... Additional arguments to \code{FUN}.
#'
#' @return A list of outputs from \code{FUN}.
#'
#' @export
#' 
lapplyOnCores <- function(x, FUN, NumberOfCores = 1L, ...) {
	# Get the number of cores to use:
	if(length(NumberOfCores) == 0) {
		NumberOfCores <- getCores()
	}
	# Do not use more cores than the number of files:
	NumberOfCores <- min(length(x), NumberOfCores)
	
	# Simple Lapply if onle one core:
	if(NumberOfCores == 1) {
		out <- lapply(x, FUN, ...)
	}
	# Run in parallel on Windows and other platforms:
	else {
		# On Windows run special args to speed up:
		if(get_os() == "win") {
			cl <- parallel::makeCluster(NumberOfCores, rscript_args = c("--no-init-file", "--no-site-file", "--no-environ"))
			parallel::clusterEvalQ(cl, {
				library(RstoxData)
			})
			out <- parallel::parLapply(cl, x, FUN, ...)
			parallel::stopCluster(cl)
		} 
		else {
			out <- parallel::mclapply(x, FUN, mc.cores = NumberOfCores, ...)
		}
	}
	
	return(out)
}


#' Run a function on all elements of x on one or more cores
#'
#' @inheritParams lapplyOnCores
#' @param ...,MoreArgs,SIMPLIFY See \code{\link[base]{mapply}}.
#'
#' @return A list of outputs from \code{FUN}.
#'
#' @export
#' 
mapplyOnCores <- function(FUN, NumberOfCores = integer(), ..., MoreArgs = NULL, SIMPLIFY = FALSE) {
	# Get the number of cores to use:
	if(length(NumberOfCores) == 0) {
		NumberOfCores <- getCores()
	}
	# Do not use more cores than the number of files:
	lll <- list(...)
	if(length(lll)) {
		NumberOfCores <- min(length(lll[[1]]), NumberOfCores)
	}
	else {
		NumberOfCores <- 1
	}
	
	# Simple mapply if onle one core:
	if(NumberOfCores == 1) {
		out <- mapply(FUN, ..., MoreArgs = MoreArgs, SIMPLIFY = SIMPLIFY)
	}
	# Run in parallel on Windows and other platforms:
	else {
		# On Windows run special args to speed up:
		if(get_os() == "win") {
			cl <- parallel::makeCluster(NumberOfCores, rscript_args = c("--no-init-file", "--no-site-file", "--no-environ"))
			out <- parallel::clusterMap(cl, FUN, ..., MoreArgs = MoreArgs, SIMPLIFY = SIMPLIFY)
			parallel::stopCluster(cl)
		} 
		else {
			out <- parallel::mcmapply(FUN, mc.cores = NumberOfCores, ..., MoreArgs = MoreArgs, SIMPLIFY = SIMPLIFY)
		}
	}
	
	return(out)
}


#' Round off to number of digits
#'
#' @param x A list of \code{data.table}s or a single \code{data.table} object.
#'
#' @return A transformed object.
#'
#' @export
#' 
setRstoxPrecisionLevel <- function(x) {
	# Get the defines number of digits:
	digits <- getRstoxDataDefinitions("digits")
	signifDigits <- getRstoxDataDefinitions("signifDigits")
	
	# If a data.table run setPrecisionLevelOneDT() directly:
	if(data.table::is.data.table(x)) {
		setPrecisionLevelOneDT(x, digits = digits, signifDigits = signifDigits)
	}
	# If a list of data tables, loop through the list and set precision:
	else if(is.list(x)) {
		for(tableName in names(x)) {
			setPrecisionLevelOneDT(x[[tableName]], digits = digits, signifDigits = signifDigits)
		}
	}
}
# Function setting the precision of one data table:
setPrecisionLevelOneDT <- function(DT, digits, signifDigits) {
	# Detect numeric columns and round off to the specified number of digits:
	atNumeric <- sapply(DT, is.numeric)
	if(any(atNumeric)) {
		numericCols <- names(DT)[atNumeric]
		# DT[, (numericCols) := round(.SD, digits), .SDcols = numericCols]
		#DT[, (numericCols) := roundSignif(.SD, digits = ..digits, signifDigits = ..signifDigits), .SDcols = numericCols]
		for(numericCol in numericCols) {
			DT[, eval(numericCol) := roundSignif(get(numericCol), digits = ..digits, signifDigits = ..signifDigits)]
		}
	}
}


roundSignif <- function(x, digits = 12, signifDigits = NULL) {
	if(length(signifDigits)) {
		digits <- pmax(signifDigits - floor(log10(abs(x))) - 1, digits)
	}
	round(x, digits)
}

## Stolen from https://stackoverflow.com/questions/47190693/count-the-number-of-integer-digits:
#n_int_digits = function(x) {
#	result = floor(log10(abs(x)))
#	result[!is.finite(result)] = 0
#	result
#}


# Function to get the formats of StoX raw data:
getStoxRawDataFormat <- function(x, unlist = FALSE) {
	formats <- lapply(x, function(this) this$metadata$useXsd)
	names(x) <- names(x)
	if(unlist) {
		formats <- unlist(formats)
	}
	return(formats)
}
	
# Check that the formats are unique:
checkUniqueFormat <- function(x) {
	nonUniqueFormats <- getRstoxDataDefinitions("nonUniqueFormats")
	uniqueFormat <- !any(getStoxRawDataFormat(x, unlist = TRUE) %in% inapplicableFormats)
	return(uniqueFormat)
}


## Function to remove rows with duplicated keys in StoxBioticData:
#removeRowsOfDuplicatedKeysFromStoxBioticData <- function(StoxBioticData) {
#	StoxBioticKeys <- getRstoxDataDefinitions("StoxBioticKeys")
#	
#	# Run through the tables of the StoxBioticData and remove duplicate rows:
#	for(tableName in names(StoxBioticData)) {
#		# Get the names of the columns which are keys:
#		presentKeys <- intersect(names(StoxBioticData[[tableName]]), StoxBioticKeys)
#		# Find rows of duplicated keys:
#		duplicatedKeys <- duplicated(StoxBioticData[[tableName]][, ..presentKeys])
#		# Remove the rows with duplicated keys:
#		rowsToKeep <- !duplicatedKeys
#		if(any(duplicatedKeys)) {
#			warning("StoX: Removing ", sum(duplicatedKeys), " rows of duplicated keys.")
#			StoxBioticData[[tableName]] <- StoxBioticData[[tableName]][rowsToKeep, ]
#		}
#	}
#	
#	return(StoxBioticData)
#}


# Function to remove rows with duplicated keys in StoxBioticData:
#' @importFrom data.table .I
removeRowsOfDuplicatedKeys <- function(StoxData, stoxDataFormat = c("Biotic", "Acoustic")) {
	
	stoxDataFormat <- match.arg(stoxDataFormat)
	StoxKeys <- getRstoxDataDefinitions(paste0("Stox", stoxDataFormat, "Keys"))
	
	# Run through the tables of the StoxData and remove duplicate rows:
	for(tableName in names(StoxData)) {
		# Get the names of the columns which are keys:
		presentKeys <- intersect(names(StoxData[[tableName]]), StoxKeys)
		# Find rows of duplicated keys:
		duplicatedKeys <- duplicated(StoxData[[tableName]], by = presentKeys)
		# Remove the rows with duplicated keys:
		if(any(duplicatedKeys)) {
			# Get the rows with equall keys, and indicate this in a copy of the data, and write to a tempfile:
			allDuplicated <- duplicated(StoxData[[tableName]], by = presentKeys) | duplicated(StoxData[[tableName]], by = presentKeys, fromLast = TRUE)
			dupData <- data.table::copy(StoxData[[tableName]])
			dupData[, duplicated := ..allDuplicated]
			dupData[, rowIndex := .I]
			fileToWriteDupDataTo <- tempfile()
			data.table::fwrite(dupData, fileToWriteDupDataTo)
			
			warning("StoX: Removing ", sum(duplicatedKeys), " rows of duplicated keys from table ", tableName, ". To see the duplicated rows run the following in R: dat <- data.table::fread(\"", fileToWriteDupDataTo, "\")")
			#rowsToKeep <- !duplicatedKeys
			StoxData[[tableName]] <- StoxData[[tableName]][!duplicatedKeys, ]
		}
	}
	
	return(StoxData)
}



AddToStoxData <- function(
	StoxData, 
	RawData, 
	VariableNames = character(), 
	NumberOfCores = 1L, 
	StoxDataFormat = c("Biotic", "Acoustic")
) {
	
	if(length(VariableNames) == 0) {
		warning("StoX: No variables specified to extract. Returning data unchcanged")
		return(StoxData)
	}
	
	# Check the the BioticData are all from the same source (ICES/NMD):
	checkDataSource(RawData)
	
	# Convert from BioticData to the general sampling hierarchy:
	StoxDataFormat <- match.arg(StoxDataFormat)
	if(StoxDataFormat == "Biotic") {
		GeneralSamplingHierarchy <- BioticData2GeneralSamplingHierarchy(RawData, NumberOfCores = NumberOfCores)
		# Define a vector of the variables to extract:
		toExtract <- c(
			getRstoxDataDefinitions("StoxBioticKeys"), 
			VariableNames
		)
	}
	else if(StoxDataFormat == "Acoustic") {
		stop("Not yet implemented")
	}
	else {
		stop("Invalid StoxDataFormat")
	}
	
	# Extract the variables to add:
	toAdd <- lapply(GeneralSamplingHierarchy, function(x) lapply(x, extractVariables, var = toExtract))
	# Rbind for each StoxBiotic table:
	toAdd <- rbindlist_StoxFormat(toAdd)
	# Extract only those tables present in StoxBioticData:
	toAdd <- toAdd[names(StoxData)]
	# Keep only unique rows:
	toAdd <- lapply(toAdd, unique)
	
	# Merge with the present StoxBioticData:
	StoxData <- mapply(merge, StoxData, toAdd)
	
	return(StoxData)
}

# Function to extracct variables from a table:
extractVariables <- function(x, var) {
	varToExtract <- intersect(names(x), var)
	if(length(varToExtract)) {
		x[, ..varToExtract]
	}
	else {
		#warning("None of the variables present")
		data.table::data.table()
	}
}

checkDataSource <- function(BioticData) {
	# Function to match the metadata against data source strings:
	matchSource <- function(x, BioticData) {
		matched <- startsWith(sapply(lapply(BioticData, "[[", "metadata"), "[[", "useXsd"), x)
		output <- rep(NA, length(matched))
		output[matched] <- x
		return(output)
	}
	
	# Detect the data source:
	possibleDataSources <- c("nmd", "ices")
	detectedDataSources <- sapply(possibleDataSources, matchSource, BioticData = BioticData, simplify = FALSE)
	numberOfFormats <- sum(sapply(detectedDataSources, function(x) any(!is.na(x))))
	#detectedDataSources <- apply(detectedDataSources, 1, min, na.rm = TRUE)
	# Accept only BioticData from a single source:
	if(numberOfFormats > 1) {
		stop("The function AddToStoxBiotic can only be applied to BioticData where all files read are of the same data source (NMD or ICES)")
	}
	
	return(detectedDataSources)
}

