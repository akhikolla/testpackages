
#' Aggregate the slot data.
#' 
#' \code{processSlots} gets slot data from a rdf list and aggregates it as 
#' specified.
#' 
#' @param slotsAnnualize A string vector with three entries.  
#'   `slotsAnnualize[1]` is the slot to process. `slotsAnnualize[2]` is the 
#'   aggregation method to use. `slotsAnnualize[3]` is the threshold or scaling 
#'   factor to use. `slotsAnnualize[4]` is the variable name to use. If 
#'   `slotsAnnualize[4]` is `NA`, then the variable is constructed as 
#'   `slotsAnnualize[1]_slotsAnnualize[2]_slotsAnnualize[3]`.
#' 
#' @param rdf The rdf list returned by [read.rdf()] to get the slot data from.  
#' 
#' @param rdfName String of the rdf name.
#' 
#' @return A data frame table with the aggregated slot data.
#' 
#' @keywords internal
#' @noRd

processSlots <- function(slotsAnnualize, rdf, rdfName, findAllSlots)
{
  ann <- slotsAnnualize[2]
	thresh <- as.numeric(slotsAnnualize[3])
	# can use thresh as a scale also.  If it is not specified, then multiply by 1.
	thresh[is.na(thresh)] <- 1 
	slot <- slotsAnnualize[1]

	if (!(slot %in% rdf_slot_names(rdf))) {
	  if (findAllSlots) {
		  stop(paste("slot:", slot, "not found in rdf:", rdfName))
	  } else {
	    # Trace Year                     Variable   Value
	    # construct a df indicating the slot couldn't be found, and return it
	    zz <- data.frame(
	      Trace = -99,
	      Year = -99,
	      Variable = ifelse(
	        is.na(slotsAnnualize[4]),
	        paste(slotsAnnualize[1],ann,thresh,sep = '_'),
	        slotsAnnualize[4]
	      ),
	      Value = -99
	    )
	    return(zz)
	  }
	}
	slot <- rdf_get_slot(rdf, slot)
	
	startData <- strsplit(rdf$runs[[1]]$start, '-')[[1]] # start year
	endData <- strsplit(rdf$runs[[1]]$end, '-')[[1]] # end year

	yy <- seq(as.numeric(startData[1]), as.numeric(endData[1]), 1)
	
	tsUnit <- rdf$runs[[1]]$time_step_unit # should either be 'year' or 'month'
	if (!(tsUnit %in% c('month','year'))) {
	  stop(
	    'rdf: ', rdfName,
	     ' contains data that is on a timestep other than year or month.\n',
	     'Currently, RWDataPlyr can only handle monthly and annual rdf data.',
	    call. = FALSE
	  )
	}
	
	if (tsUnit == 'year' & ann != 'AnnualRaw') {
	  # data is annual, so none of the aggregation methods besides annualRaw 
	  # make sense
	  warning(
	    "rdf contains annual data, but the aggregation method is not 'AnnualRaw'.\n",
	    "Processing using 'AnnualRaw' instead.\n",
	    "Edit the slotAggList and call `getDataForAllScens()` again, if necessary.",
	    call. = FALSE
	  )
	  
	  ann <- "AnnualRaw"
	}
	
	# XXX
	# Need to add other summerization methods to this area
	# XXX
	# now summarize in some way
	if(ann == 'AnnMin'){
		slot <- apply(slot, 2, trace_min_ann) # minimum annual value
		rownames(slot) <- yy
	} else if(ann == 'EOWY'){
	  slot <- slot[seq(9, nrow(slot), 12),,drop = FALSE] # 9 is september
	  slot[is.nan(slot)] <- 0
	  slot <- slot * thresh
	  rownames(slot) <- yy
	} else if(ann == 'EOCY'){
	  slot <- slot[seq(12, nrow(slot), 12),,drop = FALSE] 
		slot[is.nan(slot)] <- 0
		slot <- slot * thresh
		rownames(slot) <- yy
	} else if(ann == 'BOCY'){
	  slot <- slot[seq(1, nrow(slot), 12),,drop = FALSE] 
	  slot[is.nan(slot)] <- 0
	  slot <- slot * thresh
	  rownames(slot) <- yy
	} else if(ann == 'AnnMax'){
		slot <- apply(slot, 2, trace_max_ann) # maximum annual value
		slot <- slot * thresh
		rownames(slot) <- yy
	} else if(ann == 'AnnualSum'){
		slot <- rwslot_annual_sum(slot,thresh)
		rownames(slot) <- yy
	} else if(ann == 'AnnMinLTE'){
		slot <- apply(slot, 2, trace_min_ann) # minimum annual value
		slot <- (slot <= thresh) * 1 # convert to numeric
		rownames(slot) <- yy
	} else if(ann == 'Monthly'){
		rownames(slot) <- as.character(
		  zoo::as.yearmon(yy[1] + seq(0, (length(yy) * 12) - 1) / 12)
		)
		slot <- slot*thresh
	} else if(ann == 'WYMinLTE'){
		slot <- rbind(slot[1,],slot[1,],slot[1,],slot)
		slot <- slot[1:(nrow(slot)-3),, drop = FALSE]
		slot <- apply(slot, 2, trace_min_ann) # minimum annual value
		slot <- (slot <= thresh) * 1 # convert to numeric
		rownames(slot) <- yy
	} else if(ann == 'WYMaxLTE'){
	  slot <- rbind(slot[1,],slot[1,],slot[1,],slot)
	  slot <- slot[1:(nrow(slot)-3),, drop = FALSE]
	  slot <- apply(slot, 2, trace_max_ann) # minimum annual value
	  slot <- (slot <= thresh) * 1 # convert to numeric
	  rownames(slot) <- yy
	}	else if(ann == 'EOCYLTE'){
		slot <- slot[seq(12, nrow(slot), 12),,drop = FALSE]
		slot[is.nan(slot)] <- 0
		slot <- (slot <= thresh) * 1 # convert to numeric
		rownames(slot) <- yy
	} else if(ann == 'EOCYGTE'){
		slot <- slot[seq(12, nrow(slot), 12),, drop = FALSE]
		slot[is.nan(slot)] <- 0
		slot <- (slot >= thresh) * 1 # convert to numeric
		rownames(slot) <- yy
	} else if(ann == 'AnnualRaw'){
		if(tsUnit == 'month'){
		  # data is monthly, so will use EOCY
		  warning(
		    "User specified aggregation is 'AnnualRaw', but the rdf contains monthly data.\n",
		    "Will use 'EOCY' aggregation instead.\n",
        "If other aggregation method is desired, edit the slotAggList and call `getDataForAllScens()` again.",
		    call. = FALSE
		  )
		  
		  slot <- slot[seq(12, nrow(slot), 12),, drop = FALSE] 
		  slot[is.nan(slot)] <- 0
		  slot <- slot * thresh
		  rownames(slot) <- yy
		} else{
	    # data is annual
		  rownames(slot) <- yy
		  slot <- slot*thresh
		}
	} else{
		stop(paste0("'",ann, "'", " is an invalid aggregation method.\n",
		            "  Fix the slot aggregation list and try again."))
	}
	
	
	colnames(slot) <- seq_len(ncol(slot))

	if(ann != 'Monthly'){
		slot <- tidyr::gather(
		  tibble::rownames_to_column(as.data.frame(slot), var = "Year"), 
		  Trace, 
		  Value, 
		  -Year
		) %>%
		  dplyr::mutate(
		    Year = as.numeric(Year),
		    Trace = as.integer(Trace),
		    Variable = dplyr::if_else(
		      is.na(slotsAnnualize[4]),
		      paste(slotsAnnualize[1],ann,thresh,sep = '_'),
		      slotsAnnualize[4]
		    )
		  ) %>%
		  dplyr::select(Trace, Year, Variable, Value)
		
	} else{
		slot <- tidyr::gather(
		  tibble::rownames_to_column(as.data.frame(slot), var = "Month"), 
		  Trace, 
		  Value, 
		  -Month
		) %>%
		  dplyr::mutate(
		    Year = as.numeric(simplify2array(strsplit(Month, ' '))[2,]),
		    Month = month.name[match(
		      simplify2array(strsplit(Month, " "))[1,], 
		      month.abb
		    )],
		    Trace = as.integer(Trace),
		    Variable = dplyr::if_else(
		      is.na(slotsAnnualize[4]),
		      paste(slotsAnnualize[1],ann,thresh,sep = '_'),
		      slotsAnnualize[4]
		    )
		  ) %>%
		  dplyr::select(Trace, Month, Year, Variable, Value)
		
	}
	
	slot
}

#' Get and aggregate data from a single rdf file.
#' 
#' `getSlots()` gets all of the slots contained in a single rdf file and 
#' aggregates them as specified by the summary functions in `slotAggList`. 
#' 
#' @param slotAggList The slot aggregation list. A list containing the slots 
#'   that will be imported and aggregated, the aggregation method(s) to use, 
#'   and the rdf files that contain the slots. See [slot_agg_list()].
#' @param scenPath A relative or absolute path to the scenario folder.
#' 
#' @keywords internal
#' @noRd

getSlots <- function(slotAggList, scenPath, findAllSlots)
{
  rdf <- slotAggList$rdf
  rdf <- read.rdf(paste(scenPath,'/',rdf,sep = ''))
  
  if(slotAggList$slots[1] == 'all'){
	  # if slots is all, then need to create the slotAggList
	  # after reading in all the slot names
    slots <- rdf_slot_names(rdf)
    nSlots <- length(slots)
    if(rdf$runs[[1]]$time_step_unit == 'month'){
      aggMeth <- 'Monthly'
    } else if (rdf$runs[[1]]$time_step_unit == 'year'){
      aggMeth <- 'AnnualRaw'
    } else{
      stop(
        paste('The', slotAggList$rdf, 
              'contains data of an unexpected timestep.'),
        call. = FALSE
      )
    }
    
    slotAggList <- slot_agg_list(cbind(
      rep(slotAggList$rdf,nSlots), slots, rep(aggMeth, nSlots), rep(NA, nSlots)
    ))
    
    # go in one level into the list as that is what happens when
    # this function is called if using the normal slotAggList structure
    slotAggList <- slotAggList[[1]] 
	}
  
  slotsAnnualize <- rbind(
    slotAggList$slots, 
    slotAggList$annualize, 
    slotAggList$varNames
  )

	allSlots <- apply(
	  slotsAnnualize, 
	  2, 
	  processSlots, 
	  rdf, 
	  slotAggList$rdf, 
	  findAllSlots
	 )
	allSlots <- do.call(rbind, lapply(allSlots, function(X) X))
	allSlots
}

#' Get and aggregate data from rdf file(s) for one scenario.
#' 
#' `getAndProcessAllSlots()` gets data for a single scenario.  The slots 
#' from each rdf are processed and aggregated together.
#' 
#' @param scenPath A relative or absolute path to the scenario folder.
#' 
#' @inheritParams getDataForAllScens
#' 
#' @seealso \code{\link{getDataForAllScens}}
#' 
#' @keywords internal
#' @noRd

getAndProcessAllSlots <- function(scenPath, slotAggList, findAllSlots)
{
  sPath <- scenPath[1]
	sName <- scenPath[2]
	zz <- lapply(slotAggList, getSlots, sPath, findAllSlots)

	allRes <- do.call(rbind, lapply(zz, function(X) X))
	nn <- colnames(allRes)

	allRes$Scenario <- rep(sName, nrow(allRes))
	allRes <- subset(allRes, select=c('Scenario', nn))
	
	allRes
}

#' Get and aggregate data from an rdf file(s)
#' 
#' `getDataForAllScens()` gets slot data from multiple rdf files and/or multiple 
#' scenarios, aggregates it, and saves it as a data.frame. The slot data can be 
#' aggregated in multiple ways (see [slot_agg_list]). 
#' 
#' @param scenFolders A string vector containing the folder names (scenarios) 
#'   that the rdf files are saved in.
#' 
#' @param scenNames A string vector containing the scenario names.  This should 
#'   be the same length as `scenFolders`. The scenario names are used as 
#'   attributes to the data in the `Scenario` column.
#' 
#' @param slotAggList The slot aggregation list. Either an object of class 
#'   [slot_agg_list] or a "special" list with the keyword `"all"`. If, it is
#'   a [slot_agg_list], see that documentation for how to control the 
#'   aggregation methods used in this function. If all of the slots in an 
#'   entire rdf are desired, use a list of lists with each entry containing an 
#'   rdf file and the keyword `"all"` for the slots, e.g., 
#'   `list(list(rdf = 'KeySlots.rdf',slots = 'all'))`. If this option is used, 
#'   the function will return raw monthly, or annual data, i.e., no aggregation 
#'   methods will be applied to the data in the rdf file. 
#' 
#' @param scenPath An absolute or relative path to the folder containing 
#'   `scenFolders`.
#' 
#' @param oFile If not `NULL`, then an absolute or relative path with the file 
#'   name of the location the table will be saved to. Valid file types are 
#'   .csv, .txt, or .feather.
#' 
#' @param findAllSlots Boolean; if `TRUE` (default), then the function will
#'   abort if it cannot find a particular slot. If \code{FALSE}, then the 
#'   function will continue, even if a slot cannot be found. If a slot is not 
#'   found, then the function will return `-99` for the Trace, Year, and Value.
#'   
#' @param retFile Deprecated. Data are always returned invisibly.
#' 
#' @return A data.frame returned invisibly.
#' 
#' @examples 
#' # get a specified set of slots and apply some aggregation method to them
#' # get the data from two scenarios
#' scenFolders <- c('ISM1988_2014,2007Dems,IG,Most', 
#'   'ISM1988_2014,2007Dems,IG,2002') 
#' # slotAggTable.csv lists the slots to obtain, and the aggregation method to 
#' # apply to them
#' slotAggList <- slot_agg_list(
#'   system.file('extdata','SlotAggTable.csv',package = 'RWDataPlyr')
#' )
#' scenPath <- system.file('extdata','Scenario/',package = 'RWDataPlyr')
#' # expect Deprecated warning
#' testthat::expect_warning(
#'   keyData <- getDataForAllScens(
#'     scenFolders, 
#'     scenNames = scenFolders, 
#'     slotAggList = slotAggList, 
#'     scenPath = scenPath
#'   )
#' )
#' 
#' # get all of the data from the KeySlots rdf file
#' scenFolders <- scenFolders[1] # only one scenario
#' slotAggList <- list(list(rdf = 'KeySlots.rdf', slots = 'all'))
#' # will return monthly data for all slots in KeySlots.rdf
#' # expect Deprecated warning
#' testthat::expect_warning(
#'   allData <- getDataForAllScens(
#'     scenFolders, 
#'     scenNames = scenFolders, 
#'     slotAggList = slotAggList, 
#'     scenPath = scenPath
#'   )
#' )
#' 
#' @seealso [slot_agg_list()]
#' 
#' @export
#' 
getDataForAllScens <- function(scenFolders, scenNames, slotAggList, scenPath, 
                               oFile = NULL, retFile = NULL, findAllSlots = TRUE)
{
  .Deprecated(
    "`rw_scen_aggregate()`",
    msg = paste(
      "`getDataForAllScens()` is deprecated.",
      "Use `rw_scen_aggregate()` instead.",
      "`rw_scen_aggregate()` provides a more user friendly way of specifying",
      "and cusomizing the aggregation of RiverWare data.",
      "`getDataForAllScens()` will be removed in a future release.",
      sep = "\n"
    )
  )
  
  # determine file type to save data as:
  if (!is.null(oFile)) {
    fExt <- tools::file_ext(oFile)
    if (!(fExt %in% c('txt', 'csv', 'feather'))) {
      stop(paste0('oFile has an invalid file exention.\n',
                  'getDataForAllScens does not know how to handle ".', fExt,
                  '" extensions.'))
    }
  }
  
  if (!missing(retFile)) {
    warning(
      "In `getDataForAllScens()`, `retFile` is deprecated.\n",
      "Data are always invisibly returned.",
      call. = FALSE
    )
  }
  
	scenPath <- file.path(scenPath, scenFolders)
	scen <- cbind(scenPath, scenNames)
	zz <- apply(scen, 1, getAndProcessAllSlots, slotAggList, findAllSlots)
	zz <- do.call(rbind, lapply(zz, function(X) X))
	
	if (!is.null(oFile))
	  write_rw_data(zz, oFile)
	
	invisible(zz)
}

#' Write out csv, txt, or a feather file.
#' @noRd
write_rw_data <- function(zz, oFile)
{
  fExt <- tools::file_ext(oFile)
  if(fExt == 'txt'){
    data.table::fwrite(zz, file = oFile, row.names = FALSE, sep = '\t')
  } else if(fExt == 'csv'){
    data.table::fwrite(zz, oFile, row.names = FALSE, sep = ",")
  } else if(fExt == 'feather'){
    feather::write_feather(zz, oFile)
  }
  
  invisible(zz)
}
