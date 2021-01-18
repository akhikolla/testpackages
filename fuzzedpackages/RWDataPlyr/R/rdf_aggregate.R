
#' Aggregate RiverWare output for one or more scenarios
#' 
#' Process the user specified `rwd_agg` object for one or more scenarios to 
#' aggregate and summarize RiverWare output data.
#' 
#' `rdf_aggregate()` aggregates a single scenario of data by processing a 
#' [rwd_agg] object. 
#' 
#' In both cases, the user specifies the [rwd_agg], which 
#' determines the slots that are aggregated, and how they are aggregated. See
#' [rwd_agg] for more details on how it should be specified. 
#' 
#' See the **Directory Structure** section for how to specify `scenarios`, 
#' `scen_dir`, and `rdf_dir`.
#' 
#' @param agg A [rwd_agg] object specifying the rdfs, slots, and 
#'   aggregation methods to use.
#'   
#' @param rdf_dir The top level directory that contains the rdf files. See
#'   **Directory Structure**.
#'   
#' @param keep_cols Either boolean, or a character vector of column names to 
#'   keep in the returned tibble. The values of `keep_cols` work as follows:
#'   * `FALSE` (default) only includes the defaults columns:  
#'   `TraceNumber`, `ObjectSlot`, and `Value`. `Scenario` is also returned if 
#'   `scenario` is specified.
#'   * `TRUE`, all columns are returned.
#'   * A character vector, e.g., `c("ObjectName", "Units")`, allows the user to 
#'   include other columns that are not always required, in addition to the 
#'   "default" set of columns. If any of the values in `keep_cols` are not 
#'   found, a warning will post, but all other columns will be returned.
#'    
#' @inheritParams rdf_to_rwtbl
#' 
#' @param nans_are Either "0" or "error". If "0", then `NaN`s in the rwtbl are
#'   treated as 0s. If "error", then any `NaN`s will cause an error in this 
#'   function.
#'   
#' @param find_all_slots Boolean; if `TRUE` (default), then the function will 
#'   abort if it cannot find a particular slot. If `FALSE`, then the function 
#'   will continue, even if a slot cannot be found. If a slot is not found, 
#'   then the function will return `-99` for the Trace, and `NaN` for Year, and 
#'   Value.
#'   
#' @param cpp Boolean; if `TRUE` (default), then use [rdf_to_rwtbl2], which 
#'   relies on C++, otherwise, use original [rdf_to_rwtbl] function. 
#'   
#' @param verbose Boolean; if `TRUE` (default), then print out status of 
#'   processing the scenario(s) and the slots in each scenario.
#'   
#' @examples 
#' # rdf_aggregate() ----------
#' 
#' rdfPath <- system.file(
#'   "extdata/Scenario/ISM1988_2014,2007Dems,IG,Most", 
#'   package = "RWDataPlyr"
#' )
#' 
#' rwa <- read_rwd_agg(
#'   system.file(
#'     "extdata/rwd_agg_files/passing_aggs.csv", 
#'     package = "RWDataPlyr"
#'   )
#' )
#' 
#' x <- rdf_aggregate(rwa[1,], rdf_dir = rdfPath, scenario = "Most")
#' 
#' @export

rdf_aggregate <- function(agg, 
                          rdf_dir = ".",
                          scenario = NULL,
                          keep_cols = FALSE,
                          nans_are = "0",
                          find_all_slots = TRUE,
                          cpp = TRUE,
                          verbose = TRUE)
{
  if (!is_rwd_agg(agg))
    stop("`agg` passed to `rdf_aggregate()` is not a `rwd_agg`")
  
  nans_are <- match.arg(nans_are, choices = c("0", "error"))
  
  # check that rdf_dir is a valid directory
  if (!dir.exists(rdf_dir))
    stop("`rdf_dir` is not a valid directory")
  
  # get unique rdf files
  rdfs <- unique(agg$file)
  rdf_files <- file.path(rdf_dir, rdfs)
  rdfs_exist <- file.exists(rdf_files)
  
  # verify rdfs exist
  if (!any(rdfs_exist)) {
    stop(
      "The following rdfs were not found in ", normalizePath(rdf_dir), ":\n",
      toString(rdfs[!rdfs_exist]), 
      call. = FALSE
    )
  }
  
  rdf_len <- seq_len(length(rdfs))
  
  if (verbose) rdf_agg_msg(agg)
  
  if (cpp) {
    rwtblsmmry <- lapply(
      rdf_len,
      function(x) {
        rwtbl <- rdf_to_rwtbl2(
          rdf_files[x],
          scenario = scenario,
          keep_cols = keep_cols,
          add_ym = TRUE
        ) %>%
          check_nans(nans_are, rdf_file = rdf_files[x])
       
        tmp_sam <- agg[agg$file == rdfs[x],]
        
        rwtbl_apply_sam(rwtbl, tmp_sam, find_all_slots)
      }
    )
  } else {
  
    rwtblsmmry <- lapply(
      rdf_len,
      function(x){
        # call rwtbl_apply_sam for each unique rdf
        # seperate sam into one sam for each rdf;
        # read the rdf, then apply the sam to that rdf
        
        rwtbl <- rdf_to_rwtbl(
          read.rdf(rdf_files[x]), 
          scenario = scenario, 
          keep_cols = keep_cols, 
          add_ym = TRUE
        ) %>%
          check_nans(nans_are, rdf_file = rdf_files[x])
        
        tmp_sam <- agg[agg$file == rdfs[x],]
        
        rwtbl_apply_sam(rwtbl, tmp_sam, find_all_slots)
      }
    )
  }
  
  rwtbl_atts <- lapply(rdf_len, function(x) rwtbl_get_atts(rwtblsmmry[[x]]))
  names(rwtbl_atts) <- rdfs
  
  rwtblsmmry <- dplyr::bind_rows(rwtblsmmry)
  
  cols <- colnames(rwtblsmmry)
  cols <- cols[!(cols %in% c("Variable", "Value"))]
  rwtblsmmry <- rwtblsmmry %>% 
    dplyr::select(dplyr::one_of(cols, "Variable", "Value"))
  
  scen_folder <- data.frame(
    "scenario" = ifelse(is.null(scenario), NA_character_, scenario), 
    "folder" = normalizePath(rdf_dir),
    stringsAsFactors = FALSE
  )
  
  # save the sam as an attribute
  structure(
    rwtblsmmry,
    "rwd_agg" = agg,
    "rdf_atts" = rwtbl_atts,
    "scen_folder" = scen_folder
  )
}

#' Apply all of the operations from a rwd_agg to a single rdf file
#' @noRd

rwtbl_apply_sam <- function(rwtbl, agg, find_all_slots)
{
  # if rwd_agg uses the "all" keyword, need to construct the full rwd_agg
  # need to determine if only the all key word exists, or if there is one 
  # row that is all, but then there are summary rows
  if ("all" %in% agg$slot) {
    if (nrow(agg) == 1) {
      agg <- rwd_agg_build_all(rwtbl, agg$file)
    } else {
      # must have other summary rows, so remove all row and combine with the 
      # "all" rows
      agg <- rbind(
        agg[agg$slot != "all",],
        rwd_agg_build_all(
          rwtbl,
          agg$file[agg$slot == "all"]
        )
      )
    }
  }
  
  if (find_all_slots) {
    tmp_slots <- agg$slot %in% rwtbl_slot_names(rwtbl)
    if (!(all(tmp_slots))) {
      missing_slots <- agg$slot[!(tmp_slots)]
      stop(
        "`find_all_slots` is `TRUE`, and the following slots were not found in the ",
        agg$file[1], " file:\n",
        paste(paste("   ", missing_slots), collapse = "\n"), 
        call. = FALSE
      )
    }
  }
  
  sam_rows <- seq_len(nrow(agg))
  rwtblsmmry <- lapply(
    sam_rows, 
    function(x) rwtbl_apply_sar(rwtbl, agg[x,])
  )
  
  rwtblsmmry <- dplyr::bind_rows(rwtblsmmry)
  
  attributes(rwtblsmmry) <- c(attributes(rwtblsmmry), rwtbl_get_atts(rwtbl))
  
  rwtblsmmry
}

#' Apply a single row's aggregation (sar) from an rwd_agg object
#' @noRd
rwtbl_apply_sar <- function(rwtbl, slot_agg_row)
{
  # if the slot cannot be found, return -99, and NaN, NaN, since if we have 
  # gotten here, `find_all_slots`, must be `FALSE`
  
  if (slot_agg_row$slot %in% rwtbl_slot_names(rwtbl)) {
  zz <- apply_period(rwtbl, slot_agg_row) %>%
    apply_summary(slot_agg_row) %>%
    apply_eval(slot_agg_row) %>%
    add_month_to_annual() %>%
    add_var_drop_objectslot(slot_agg_row)
  } else{
    zz <- build_missing_slot_values(slot_agg_row) %>%
      add_month_to_annual() %>%
      add_var_drop_objectslot(slot_agg_row)
  }
  
  zz
}

#' Add in a `Variable` columns and drop the `ObjectSlot` column
#' @noRd

add_var_drop_objectslot <- function(rwtbl, slot_agg_row)
{
  tmp_groups <- dplyr::group_vars(rwtbl)
  tmp_groups <- tmp_groups[tmp_groups != "ObjectSlot"]
  
  rwtbl$Variable <- slot_agg_row$variable
  
  rwtbl %>%
    dplyr::group_by_at(tmp_groups) %>%
    dplyr::select(-dplyr::matches("ObjectSlot"))
}

#' If the data had dropped the Month column, add in December so that it is 
#' considered annual data.
#' 
#' @noRd

add_month_to_annual <- function(rwtbl)
{
  cols <- colnames(rwtbl)
  dec_fun <- function(x) "December"
  if (!"Month" %in% cols){
    rwtbl <- rwtbl %>%
      dplyr::mutate_at("Year", list("Month" = dec_fun)) %>%
      dplyr::select_at(c("Year", "Month", cols[cols != "Year"]))
  }
  rwtbl
}

#' Check if there are any NaNs in the `rwtbl`, and either convert to 0s or 
#' throw an error.
#' 
#' @inheritParams rdf_aggregate
#' @param rdf_file The rdf file name as a character.
#' 
#' @noRd
check_nans <- function(rwtbl, nans_are, rdf_file)
{
  nans_are <- match.arg(nans_are, choices = c("0", "error"))
  
  if (any(is.nan(rwtbl$Value))) {
    if (nans_are == "error") {
      slots <- rwtbl %>%
        dplyr::filter_at("Value", dplyr::all_vars(is.nan(.)))
      
      slots <- unique(slots$ObjectSlot)
      nSlots <- length(slots)
      slotM <- ""
      if (length(slots) > 10) {
        slots <- slots[1:10]
        slotM <- "...\n(Only the first 10 slots containing `NaN`s are printed.)"
      }
      
      stop(
        "`NaN`s were found in ", rdf_file, 
        " and `nans_are` treated as an error.\n",
        "`NaN`s were found in ", nSlots, " slots:\n",
        paste(slots, collapse = "\n"),
        "\n", slotM,
        call. = FALSE
      )
      
    } else {
      # convert any NaNs to 0
      rwtbl <- rwtbl %>%
        dplyr::mutate_at("Value", nan_to_zero)
    }
  }
  
  rwtbl
}

#' Convert NaNs to 0
#' @noRd
nan_to_zero <- function(x) {
  x[is.nan(x)] <- 0
  x
}

#' Build a single row tibble that contains the appropriate values for slots
#' that cannot be found
#' @noRd
build_missing_slot_values <- function(slot_agg_row) 
{
  tibble::tibble(
    Year = NaN,
    TraceNumber = -99,
    ObjectSlot = slot_agg_row$slot,
    Value = NaN
  )
}

rdf_agg_msg <- function(agg)
{
  rdfs <- unique(agg$file)
  
  lapply(seq_len(length(rdfs)), function(x) {
    rwa2 <- agg[agg$file == rdfs[x],]
    if("all" %in% rwa2$slot){
      message("   Processing all slots in ", rdfs[x])
    } else{
      slot <- if(nrow(rwa2) == 1) "slot" else "slots"
      message("   Processing ", nrow(rwa2), " ", slot, " in ", rdfs[x])
    }
  })
  
  invisible(agg)
}
