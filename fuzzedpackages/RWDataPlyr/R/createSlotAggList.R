
#' Creates a list for use by \code{\link{getDataForAllScens}}.
#' 
#' Deprecated: please use [slot_agg_list()] instead, which returns the same 
#' list, but now as a "slot_agg_list" object.
#' 
#' @param iData Either an N x 4 character matrix or a character with an absolute 
#'   or relative path to a csv file.
#' 
#' @export
createSlotAggList <- function(iData)
{
  .Deprecated("slot_agg")
  create_slot_agg_list(iData)
}

#' creates the slot_agg object; see ?slot_agg.
#' @noRd

create_slot_agg_list <- function(iData)
{
  if (!is.matrix(iData)) {
    if (length(iData) > 1) {
      if (length(iData) %% 4 == 0) {
        warning(
          "Attempting to convert `iData` to a N x 4 matrix. Results may be unexpected.\n",
          "Probably better to stop and pass a matrix to `slot_agg_list()`.",
          call. = FALSE
        )
        iData <- matrix(iData, ncol = 4, byrow = TRUE)
      } else if (length(iData) %% 5 == 0){
        warning(
          "Attempting to convert `iData` to a N x 5 matrix. Results may be unexpected.\n",
          "Probably better to stop and pass a matrix to `slot_agg_list()`.",
          call. = FALSE
        )
        iData <- matrix(iData, ncol = 5, byrow = TRUE)
      } else {
        stop(
          "`iData` is not a matrix, nor can it be converted to an Nx4 or Nx5 matrix",
          call. = FALSE
        )
      }
    } else if (!file.exists(iData)) {
      stop(paste(iData,'does not exist.'), call. = FALSE)
    } else {
      iData <- as.matrix(utils::read.csv(iData, header = FALSE))
    }
  } else {
    # it is a matrix
    # if it is a matrix, make sure it has 4 or 5 columns
    if(ncol(iData) != 4 & ncol(iData) != 5) {
      stop(
        "`iData` is a matrix with ", ncol(iData), 
        " columns. There should either be 4 or 5 columns.",
        call. = FALSE
      )
    }
  }
  
  # check and see if the "monthly" aggregation method exists, if it does, it 
  # should be the only aggregation method
  # this is contained in column 3
  if ("Monthly" %in% iData[,3] & !all(iData[,3] == "Monthly")) {
    stop(
      "The \"Monthly\" aggregation method cannot currently be mixed with other aggregation methods.\n",
      "Please create a seperate slot aggregation list with only the monthly data.",
      call. = FALSE
    )
  }
  
  # make sure that all of the slot agg methods are valid
  if (!all(iData[,3] %in% slotAggMethods())) {
    tmp <- iData[!(iData[,3] %in% slotAggMethods()),3]
    stop(
      paste0("Invalid aggregation methods:\n    ", 
      paste(tmp, collapse = ", "), "\n  ",
      paste("Fix the", length(tmp), "aggregation method(s) and try again.")),
      call. = FALSE
    )
  }
  
  #check that all variable names are unique
  if (ncol(iData) == 5 && length(unique(iData[,5])) < nrow(iData)) {
    stop(
      "All variable names passed to `slot_agg_list()` must be unique.",
      call. = FALSE
    )
  }
  
  # check and see if alternative variable names have been added
  altNames <- ncol(iData) == 5
  
  sl <- list() #slot list
  # create one entry for each unique rdf
  rdfs <- levels(as.factor(iData[,1]))
  
  for (i in seq_len(length(rdfs))) {
    tmp <- matrix(iData[which(rdfs[i]==iData[,1]),],ncol = dim(iData)[2])
    sl[[i]] <- list()
    sl[[i]]$rdf <- rdfs[i]
    sl[[i]]$slots <- tmp[,2]
    sl[[i]]$annualize <- matrix(t(tmp[,3:4]),nrow = 2)
    if (altNames) {
      sl[[i]]$varNames <- tmp[,5]      
    } else {
      sl[[i]]$varNames <- rep(NA,nrow(tmp))
    }
  }
  
  sl
}

#' @keywords internal 
slotAggMethods <- function() {
  c("AnnMin", "AnnMax", "AnnualSum", "AnnMinLTE", "AnnualRaw", "BOCY", 
    "EOCY", "EOCYGTE", "EOCYLTE", "EOWY", "Monthly", "WYMaxLTE", "WYMinLTE")
}
