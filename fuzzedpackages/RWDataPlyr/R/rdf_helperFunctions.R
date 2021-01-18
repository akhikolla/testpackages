
# ****************************************************************************
# 								rdf_slot_names
# ****************************************************************************
#' Returns all slots contained in an rdf file.
#' 
#' `rdf_slot_names()` returns a character vector of all slots contained within 
#' an rdf object.
#' 
#' @param rdf An `rdf` object. 
#' 
#' @return A character vector.
#' 
#' @examples
#' rdf_slot_names(keyRdf)
#' 
#' @seealso [read.rdf()]
#' 
#' @export
rdf_slot_names <- function(rdf)
{
  stopifnot(methods::is(rdf, "rdf"))
  
  names(rdf$runs[[1]]$objects)
}

#' @describeIn rdf_slot_names Deprecated version of `rdf_slot_names()`
#' @export
getSlotsInRdf <- function(rdf)
{
  .Deprecated("rdf_slot_names()")
  rdf_slot_names(rdf)
}

# ****************************************************************************
#                             rdf_get_slot
# ****************************************************************************
#' Get a slot out of an rdf object
#' 
#' `rdf_get_slot()` gets a slot from an `rdf` object and creates a matrix with
#' rows indexing through time and columns indexing over traces.
#' 
#' @param rdf An `rdf` object.
#' @param slot Character slot name that exists in the `rdf`.
#' 
#' @examples
#' pe <- rdf_get_slot(keyRdf, "Mead.Pool Elevation")
#' 
#' @return A matrix with traces as columns and timesteps as rows.
#' 
#' @export
rdf_get_slot <- function(rdf, slot)
{
  stopifnot(methods::is(rdf, "rdf"))
  stopifnot(length(slot) == 1)
  
  # check to see if the slot exists in the rdf, if it does not exit error out
  if (!(slot %in% rdf_slot_names(rdf)))
     stop(paste(slot,'not found in rdf:',deparse(substitute(rdf))))
     
	nn <- as.numeric(rdf$meta$number_of_runs)
	res <- do.call(
	  cbind, 
	  lapply(1:nn, function(xx) rdf$runs[[xx]]$objects[[slot]]$values)
	)
	
	attr(res, "timespan") <- rdf_get_timespan(rdf)
	
	res
}

#' @describeIn rdf_get_slot Deprecated version of `rdf_get_slot()`
#' @export
rdfSlotToMatrix <- function(rdf, slot)
{
  .Deprecated("rdf_get_slot()")
  rdf_get_slot(rdf, slot)
}

# ****************************************************************************
# **************************  getTimeSpan  *******************************
# ****************************************************************************
#' Returns the simulation timespan from an rdf 
#' 
#' `rdf_get_timespan()` gets the simulation timespan from an `rdf` object.
#' 
#' @param rdf An `rdf` object (likely from [read.rdf()]).
#' 
#' @return A named character vector with two elements. The first element, named 
#'   `"start"`, includes the start date of the simulation. The second element, 
#'   named `"end"`, includes the end date of the simulation.
#'   
#' @examples
#' rdf_get_timespan(keyRdf)
#' 
#' @export
rdf_get_timespan <- function(rdf)
{
  stopifnot(methods::is(rdf, "rdf"))
  
  c('start' = rdf$runs[[1]]$start, 'end' = rdf$runs[[1]]$end)
}

getTimeSpan <- function(rdf)
{
  .Deprecated("rdf_get_timespan()")
  rdf_get_timespan(rdf)
}

# ****************************************************************************
# **************************  rdfSlotToXTS  **********************************
# ****************************************************************************
#' Get one slot out of an rdf list and put it in an XTS object
#' 
#' \code{rdfSlotToXTS} Takes a list created by \code{\link{read.rdf}} and 
#' convert the nested slot values over the multiple traces into an XTS array 
#' indexing over traces.
#' 
#' @param rdf list returned by \code{\link{read.rdf}}
#' @param slot string of slot name that exists in \code{rdf} that will be 
#'   converted to a matrix
#' @return an XTS object with the selected slot data
#' @examples
#' pe <- RWDataPlyr:::rdfSlotToXTS(keyRdf, 'Mead.Pool Elevation')
#' 
#' @keywords internal

rdfSlotToXTS <- function(rdf, slot)
{
  # check to see if the slot exists in the rdf, if it does not exit
  if(!(slot %in% rdf_slot_names(rdf)))
    stop(paste(slot,'not found in rdf:',deparse(substitute(rdf))))
  # Get date-times from rdf
  tArray <- rdf$runs[[1]]$times
  # OPERATIONS IN ORDER OF EXECUTION
  # 1. rdfSlotToMatrix - read data for 'slot' string given 'rdf' file
  # 2. cbind - combine datetime and data arrays
  # 3. data.frame - define R dataframe for conversion to XTS
  # 4. read.zoo - convert dataframe to zoo matrix
  # 5. as.xts - convert zoo matrix to XTS
  # 6. Storage.mode() - convert char values in the XTS matrix to numeric
  rdfXTS <- xts::as.xts(zoo::read.zoo(data.frame(cbind(
    tArray,
    rdf_get_slot(rdf, slot)
  ))))
  storage.mode(rdfXTS) <- "numeric"
  runNames <- c()
  for (ithRun in c(1:as.numeric(rdf$meta$number_of_runs))){
    runNames <- c(runNames, paste('Trace',ithRun,sep=""))
  }
  names(rdfXTS) <- runNames
  rdfXTS
}
