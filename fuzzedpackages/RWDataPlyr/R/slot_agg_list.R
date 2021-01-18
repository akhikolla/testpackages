#' A class to control how RiverWare data are aggregated
#' 
#' `"slot_agg_list"` is a class that contains a set of RiverWare slots, which 
#' rdf file they are found in, and a set of keywords that are used to control 
#' how [getDataForAllScens()] aggregates RiverWare data. 
#' 
#' @details 
#' The `slot_agg_list` class, contains a list that includes: which rdf file to 
#' find each slot in, how to aggregate and process the slot data, and any 
#' thresholds or scaling factors. The function can either read in a csv file or 
#' start from an N x 4 or N x 5 string matrix (the 5th column is optional).
#' 
#' The csv file and the matrix should be in the form of an Nx4 or Nx5 matrix.  
#' Each row is a single slot, aggregation, and threshold combination. If you 
#' want to compare a single slot value to multiple thresholds, it needs to have 
#' one row for each threshold. The first column is the rdf the slot is found in. 
#' The second column is the slot name. The third column is the aggregation 
#' method that will be applied to the slot (see below for a list of the 
#' aggregation methods).  The fourth column is a scaling factor or threshold to 
#' compare the slot data to. The fifth column is an optional column; if 
#' specified, the 5th column will be used for the variable name for the 
#' data.frame created by [getDataForAllScens()]. If it is not specified
#' the variable name will be created by concatenating the slot, aggregation 
#' method, and threshold/scaling factor using '_' to separate the columns. 
#' Below is an example table. All values should be strings except for `NA`, 
#' if specified as a matrix in R.
#' 
#' \tabular{ccccc}{
#' \strong{rdf} \tab \strong{Slot} \tab \strong{Aggregation Method} \tab 
#' \strong{Threshold or Scaling Factor} \tab \strong{Variable Name (optional)}\cr
#' 'KeySlots.rdf' \tab 'Mead.Pool Elevation' \tab 'EOCY' \tab NA \tab Mead EOCY Elevation\cr
#' 'KeySlots.rdf' \tab 'Mead.Pool Elevation' \tab 'AnnMinLTE' \tab '1100' \tab Mead < 1,100\cr
#' 'KeySlots.rdf' \tab 'Mead.Pool Elevation' \tab 'AnnMinLTE' \tab '1060' \tab Mead < 1,060\cr
#' 'Other.rdf' \tab 'Powell.Outflow' \tab 'AnnualSum' \tab '0.001' \tab Powell Annual Release\cr
#' }
#' 
#' The above table lists each slot, the rdf the slot is saved in, the summary 
#' function, the threshold to be used to scale the data by or compare the data 
#' to, and an optionally specified variable name. The threshold and scaling 
#' factors are described in more detail below. For example, the first row will 
#' result in compiling all end-of-December values for Mead's pool elevation.  
#' The data will not be scaled, and [getDataForAllScens()] will look in 
#' KeySlots.rdf for the "Mead.Pool Elevation" slot. The second row will find 
#' the annual minimum Mead pool elevation and see if it is less than or equal 
#' to 1,100' feet in the second line and less than or equal to 1,060' feet in 
#' the third row. 
#' 
#' To scale the data by a value less than 1, use decimals rather than fractions, 
#' as shown in the fourth row. If the Variable Name column was not specified, 
#' the variable name for the first row would be `Mead.Pool Elevation_EOCY_1`
#' `NA` is replaced with a 1 when constructing the variable names.
#' 
#' See the **Aggregation Methods** section for available aggregation methods.
#' 
#' @section Aggregation Methods:
#' 
#' The available aggregation methods are as follows. The behavior of the 
#' "Threshold or scaling factor" are described and a bold **"Threshold"** or 
#' **"Scaled"** indicates which is used by the aggregation method. For scaling 
#' factors, a value of `NA` will not scale the data.
#' 
#' \describe{
#'  \item{'AnnMin'}{Returns the minimum annual \strong{scaled} value.}
#'  \item{'AnnMax'}{Returns the maximum annual \strong{scaled} value.}
#'  \item{'AnnualSum'}{Returns the annual \strong{scaled} sum.}
#'  \item{'AnnMinLTE'}{Checks to see if the annual minimum value is less than or 
#'  equal to a \strong{threshold.} Returns 1 if it is less than or equal to 
#'  the \strong{threshold} and 0 otherwise.}
#'  \item{'AnnualRaw'}{Returns the annual \strong{scaled} data. This aggregation 
#'  method should only be used if the rdf file contains only annual data. For 
#'  rdf files that include monthly data and only an annual value is desired, the 
#'  EOCY aggregation method should be used. This differs from the Monthly 
#'  aggregation method, only in the timestep naming.}
#'  \item{'BOCY'}{Beginning-of-calendar year values are reported and 
#'  \strong{scaled}. Any values that are NaNs are changed to 0s.}
#'  \item{'EOCY'}{End-of-calendar year values are reported and \strong{scaled}. 
#'  Any values that are NaNs are changed to 0s.}
#'  \item{'EOCYGTE'}{Checks to see if the end-of-calendar year values are 
#'  greater than or equal to a \strong{threshold}. Returns 1 if it is greater 
#'  than or equal to the \strong{threshold} and 0 otherwise.}
#'  \item{'EOCYLTE'}{Checks to see if the end-of-calendar year values are less 
#'  than or equal to a \strong{threshold}. Returns 1 if it is less than or 
#'  equal to the \strong{threshold} and 0 otherwise.}
#'  \item{'EOWY'}{End-of-water year values are reported and \strong{scaled}. 
#'  Any values that are NaNs are changed to 0s.}
#'  \item{'Monthly'}{Returns the monthly \strong{scaled} data.}
#'  \item{'WYMaxLTE'}{Checks to see if the maximum water year value is less than 
#'  or equal to a \strong{threshold.} Returns 1 if it is less than or equal to 
#'  the \strong{threshold} and 0 otherwise. This can be used to determine if an 
#'  entire water year is below a \strong{threshold}. The water year is defined 
#'  as October through September of the next year. For the first year, only 
#'  January through September are evaluated as RiverWare does not typically 
#'  export pre-simulation data.}
#'  \item{'WYMinLTE'}{Checks to see if the minimum water year value is less than 
#'  or equal to a \strong{threshold.} Returns 1 if it is less than or equal to 
#'  the \strong{threshold} and 0 otherwise. The water year is defined as October 
#'  through September of the next year. For the first year, only January through 
#'  September are evaluated as RiverWare does not typically export 
#'  pre-simulation data.}
#' }
#' 
#' @param x Either an Nx4 character matrix or a character with an absolute 
#'   or relative path to a csv file.
#'   
#' @return A  `slot_agg_list`  object.
#' 
#' @examples
#' # read in a csv file that contains the data
#' slot_agg_list(
#'   system.file('extdata','SlotAggTable.csv',package = 'RWDataPlyr')
#' )
#' 
#' # or specify as a matrix
#' slot_agg_matrix <- matrix(
#'   c("KeySlots.rdf", "Powell.Outflow", "AnnualSum", ".001", "powellAnnRel", 
#'   "KeySlots.rdf", "Mead.Pool Elevatoin", "AnnMinLTE", "1050", "meadLt1050"), 
#'   nrow = 2, 
#'   byrow = TRUE
#' )
#' slot_agg_list(slot_agg_matrix)
#' 
#' @seealso [getDataForAllScens()]
#'
#' @export

slot_agg_list <- function(x)
{
  structure(
    create_slot_agg_list(x),
    class = "slot_agg_list"
  )
}

#' @export
print.slot_agg_list <- function(x, ...)
{
  printSlot <- function(x, i, rdfI){
    thresh <- ifelse(
      is.na(x[[rdfI]]$annualize[2, i]), 
      "NA", 
      x[[rdfI]]$annualize[2, i]
    )
    
    var <- ifelse(
      is.na(x[[rdfI]]$varNames[i]) || x[[rdfI]]$varNames[i] == "", 
      "default constructed", 
      x[[rdfI]]$varNames[i]
    )
    
    cat("  ", x[[rdfI]]$annualize[1, i], "(", x[[rdfI]]$slots[i], ",", thresh, 
        ")",  "->", var ,"\n")
  }
  
  printRdf <- function(x, i){
    cat(x[[i]]$rdf, ":\n", sep = "")
    nslots <- seq_len(length(x[[i]]$slots))
    lapply(nslots, function(nn) printSlot(x, nn, i))
  }
  
  nRdf <- seq_len(length(x))
  lapply(nRdf, function(i) printRdf(x, i))
  invisible(x)
}

#' @export
summary.slot_agg_list <- function(object, ...)
{
  cat(length(object[[1]]$slots), "slots to aggregate from", length(object), 
      "rdf(s).")
}

#' Test if the object is a slot_agg_list
#'
#' @param x An object
#' 
#' @return `TRUE` if the object inherits from the `slot_agg_list` class.
#' @export
is_slot_agg_list <- function(x) inherits(x, "slot_agg_list")

#' @rdname is_slot_agg_list
#' @export
is.slot_agg_list <- is_slot_agg_list
