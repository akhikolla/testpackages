#' Run filter on any StoX related data source
#'
#' @name filterData
#'
#' @param inputData An input data. Can be a list of biotic data (StoX data type \code{\link{BioticData}}), list of acoustic data, StoxBiotic data, or StoxAcoustic data.
#' @param filterExpression Filter expression in list of strings. The name of the list and structures should be identical to the names of the input data list.
#' @param propagateDownwards Whether the filter action will propagate in the downwards direction. Default to TRUE.
#' @param propagateUpwards Whether the filter action will propagate in the upwards direction. Default to FALSE.
#'
#' @return An object of filtered data in the same format as the input data.
#'
#' @importFrom utils head
#' @importFrom data.table fsetdiff is.data.table %like% %flike% %ilike% %inrange% %chin% %between%
#' @export
#' 
filterData <- function(inputData, filterExpression, propagateDownwards = TRUE, propagateUpwards = FALSE) {
    
	`%notin%` <- Negate(`%in%`)

	processFilter <- function(filters) {
		# Assume each individual filters relation are the AND (&) operator 
		parsedFilterTxt <- paste(filters, collapse=" & ")
		parsedFilter <- parse(text = parsedFilterTxt)
		return(parsedFilter)
	}

	processFilterWrap <- function(x) {
		out <- lapply(x, processFilter)
		names(out) <- names(x)
		return(out)
	}

	applyFilter <- function(tableName, x, y, propDown, propUp) {
		ret <- list()

		filts <- x[[tableName]]

		if(!nrow( y[[tableName]])) {
			warning("StoX: Empty table ", tableName)
			return(ret)
		}
		
		test <- try(ret[[tableName]] <- y[[tableName]][eval(filts),], silent = TRUE)
		if(class(test)[1] == "try-error") {
			warning("StoX: Apply filter error:\n", test[1])
		} else {
			# 3. propagate down
			if(propDown) {
				start <- which(names(y) == tableName)
				range <- c(start:length(names(y)))
				# Don't propage if it's the only table
				if(length(range) > 1) {
					for(parent in head(range, -1)) {
						# Sometimes we need to hop (some) empty tables that were originally empty (e.g., NMD's prey)
						goUp <- 0
						while(nrow(y[[parent - goUp]]) == 0 && parent - goUp > 1) {
							goUp <- goUp + 1
						}
						# Find similar columns (assume those are keys)
						key <- intersect(names(y[[parent + 1]]), names(y[[parent - goUp]]))
						if(length(key) > 0) {
							# Find the not deleted keys after filtering
							deleted <- fsetdiff(unique(y[[parent - goUp]][, ..key]), unique(ret[[names(y)[parent - goUp]]][, ..key]))
							# Propagate to child
							ret[[names(y)[parent+1]]] <- y[[parent+1]][!deleted, on = names(deleted)]
						}
					}
				}
			}
			# 4. propagate up
			if(propUp) {
				start <- which(names(y) == tableName)
				range <- c(1:start)
				# Don't propage if it's the only table
				if(length(range) > 1) {
					for(parent in head(rev(range), -1)) {
						# Sometimes we need to hop (some) empty tables that were originally empty (e.g., NMD's prey)
						goDown <- 0
						while(nrow(y[[parent + goDown]]) == 0 && parent + goDown <= length(y)) {
							goDown <- goDown + 1
						}
						# Find similar columns (assume those are keys)
						key <- intersect(names(y[[parent - 1]]), names(y[[parent + goDown]]))
						if(length(key) > 0) {
							# Find the not deleted keys after filtering
							deleted <- fsetdiff(unique(y[[parent + goDown]][, ..key]), unique(ret[[names(y)[parent + goDown]]][, ..key]))
							# Propagate to parent
							ret[[names(y)[parent-1]]] <- y[[parent-1]][!deleted, on = names(deleted)]
						}
					}
				}
			}
		}
		return(ret)
	}

	applyFilterWrap <- function(fileName, x, y, propDown, propUp) {
		# Do per file filtering
		merged <- y[[fileName]]
		for (tName in intersect(names(merged), names(x[[fileName]]))) {
			out <- applyFilter(tName, x[[fileName]], merged, propDown, propUp)
			if(!length(out)) {
				warning("StoX: Filter on data from file \"", fileName, "\" returned empty table \"", tName, "\"")
			}
			merged <- replace(merged, intersect(names(out), names(merged)), out[intersect(names(out), names(merged))])
		}

		return(merged)
	}

	level <- 0
	
	# 1. Check Validity/Level of data
	if(!length(filterExpression)) {
		return(inputData)
	} else if(!is.list(filterExpression)) {
		warning("StoX: Invalid filter parameter (must be a list)!")
		return(NULL)
	} else if(!is.list(inputData) || !length(inputData)) {
		warning("StoX: Invalid or empty input data!")
		return(NULL)
	} else if(is.data.table(inputData[[1]])) {
		level <- 1
		if(!is.character(filterExpression[[1]])) {
			warning("StoX: Data/Filter level mismatch!")
			return(NULL)
		}
	} else if(is.data.table(inputData[[1]][[1]])) {
		level <- 2
		if(!is.character(filterExpression[[1]][[1]])) {
			warning("StoX: Data/Filter level mismatch!")
			return(NULL)
		}
	} else {
		warning("StoX: Something wrong with the input!")
		return(NULL)
	}

	# 2. Check and parse filters
	if(level == 1) {
		pFilters <- lapply(filterExpression, processFilter)
	} else {
		pFilters <- lapply(filterExpression, processFilterWrap)
	}
	names(pFilters) <- names(filterExpression)
	
	# 3. Apply filters
	if(level == 1) {
		merged <- inputData
		for (tName in intersect(names(merged), names(pFilters))) {
			out <- applyFilter(tName, pFilters, merged, propagateDownwards, propagateUpwards)
			merged <- replace(merged, intersect(names(out), names(merged)), out[intersect(names(out), names(merged))])
		}
	} else {
		ret <- lapply(names(pFilters), applyFilterWrap, pFilters, inputData, propagateDownwards, propagateUpwards)
		names(ret) <- names(pFilters)
		merged <- replace(inputData, intersect(names(ret), names(inputData)), ret[intersect(names(ret), names(inputData))])
	}

	return(merged)
}




expandFilterExpressionList <- function(FilterExpressionList, sep = "/") {
	
	# Error if not a list:
    if(!is.list(FilterExpressionList)) {
        #stop("FilterExpressionList must be a list")
        return(FilterExpressionList)
    }
	else if(length(FilterExpressionList) == 0) {
        return(FilterExpressionList)
	}
    # If the input list of expressions has 2 levels, return immediately:
    if(is.list(FilterExpressionList[[1]])) {
        return(FilterExpressionList)
    }
	
    # Get the file names and the table names:
    splited <- strsplit(names(FilterExpressionList), split = sep)
    fileNames <- sapply(splited, function(x) x[seq_len(length(x) - 1)])
    tableNames <- sapply(splited, utils::tail, 1)
    tableNames <- split(tableNames, fileNames)
    
    # Split the FilterExpression by the fileNames:
    FilterExpressionList <- split(FilterExpressionList, fileNames)
    
    names(FilterExpressionList) <- unique(fileNames)
    
    # Change the names of the individual tables:
    for(ind in seq_along(FilterExpressionList)) {
        names(FilterExpressionList[[ind]]) <- tableNames[[ind]]
    }
    
    return(FilterExpressionList)
}


#' Filter (raw) Biotic data
#' 
#' Filters \code{\link{BioticData}}.
#'
#' @param BioticData  Input \code{\link{BioticData}} data.
#' @param FilterExpression Filter expression given as a list of strings. The name of the list and structures should be identical to the names of the input data list. To extract or exclude missing values (NAs) use the operator \code{\%in\%} or the special operator \code{\%notin\%}, which is defined in RstoxData.
#' @param FilterUpwards Whether the filter action will propagate in the upwards direction. Default to FALSE.
#'
#' @return An object of filtered data in the same format as the input data.
#'
#' @export
#' 
FilterBiotic <- function(BioticData, FilterExpression, FilterUpwards = FALSE) {
    # For filtering directly on the input data, we need to split the list filter expression to one level for the file and one for the table:
    FilterExpression <- expandFilterExpressionList(FilterExpression)
    
    filterData(
        BioticData, 
        filterExpression = FilterExpression, 
        propagateDownwards = TRUE, 
        propagateUpwards = FilterUpwards
    )
}

#' Filter (raw) Acoustic data
#'
#' Filters \code{\link{AcousticData}}.
#' 
#' @param AcousticData  Input \code{\link{AcousticData}} data.
#' @param FilterExpression Filter expression in list of strings. The name of the list and structures should be identical to the names of the input data list.
#' @param FilterUpwards Whether the filter action will propagate in the upwards direction. Default to FALSE.
#'
#' @return An object of filtered data in the same format as the input data.
#'
#' @export
#' 
FilterAcoustic <- function(AcousticData, FilterExpression, FilterUpwards = FALSE) {
    # For filtering directly on the input data, we need to split the list filter expression to one level for the file and one for the table:
    FilterExpression <- expandFilterExpressionList(FilterExpression)
    
    filterData(
        AcousticData, 
        filterExpression = FilterExpression, 
        propagateDownwards = TRUE, 
        propagateUpwards = FilterUpwards
    )
}


#' Filter StoxBiotic data
#'
#' Filters \code{\link{StoxBioticData}}.
#' 
#' @param StoxBioticData  Input \code{\link{StoxBioticData}} data.
#' @param FilterExpression Filter expression in list of strings. The name of the list and structures should be identical to the names of the input data list.
#' @param FilterUpwards Whether the filter action will propagate in the upwards direction. Default to FALSE.
#'
#' @return An object of filtered data in the same format as the input data.
#'
#' @export
#' 
FilterStoxBiotic <- function(StoxBioticData, FilterExpression, FilterUpwards = FALSE) {
    filterData(
        StoxBioticData, 
        filterExpression = FilterExpression, 
        propagateDownwards = TRUE, 
        propagateUpwards = FilterUpwards
    )
}

#' Filter StoxAcoustic data
#'
#' Filters \code{\link{StoxAcousticData}}.
#' 
#' @param StoxAcousticData  Input \code{\link{StoxAcousticData}} data.
#' @param FilterExpression Filter expression in list of strings. The name of the list and structures should be identical to the names of the input data list.
#' @param FilterUpwards Whether the filter action will propagate in the upwards direction. Default to FALSE.
#'
#' @return An object of filtered data in the same format as the input data.
#'
#' @export
#' 
FilterStoxAcoustic <- function(StoxAcousticData, FilterExpression, FilterUpwards = FALSE) {
    filterData(
        StoxAcousticData, 
        filterExpression = FilterExpression, 
        propagateDownwards = TRUE, 
        propagateUpwards = FilterUpwards
    )
}
