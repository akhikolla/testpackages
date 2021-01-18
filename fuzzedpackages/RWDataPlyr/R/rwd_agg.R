
#' Class to specify how to aggregate RiverWare data
#' 
#' `rwd_agg()` creates a RiverWare data aggregator (`rwd_agg`) object, which 
#' lets users specify how specific RiverWare slots should be aggregated.
#' 
#' `rwd_agg` objects can be created in three ways:
#' 1. By providing a data.frame, with the following expected columns `file`, 
#'   `slot`, `period`, `summary`, `eval`, `t_s`, and `variable`. Each row in the
#'   data.frame should include all of the information for how each individual
#'   slot will be aggregated. See the *Aggregation Specification* section below
#'   for details on how to specify each column. 
#' 2. By providing a vector of rdf files. If specified in this manor, all of the
#'    slots in each rdf file will be read in to a `rwtbl`, but will not be 
#'    aggregated/summarized.
#'    
#'    In this case, the `variable` names are automatically constructed from the
#'    `ObjectSlot` names. The variable names are constructed as the all lower 
#'    case version of the object_slot name. If the slot name is "Pool Elevation", 
#'    it is shortened to "pe", otherwise the full object and slot name are used. 
#'    If there are any spaces, they are replaced with underscores.
#'    
#'  3. By reading in a csv file with [read_rwd_agg()]. This csv file must have 
#'     the correct column names and meet other requirements described below. To
#'     ensure it has the correct column names, [rwd_agg_template()] can be used
#'     to create a blank csv file for the user to fill in. 
#'   
#' @section Aggregation Specification:
#' 
#' In order to specify how each slot should be aggregated, each column should 
#' include specific keywords, which are described below. It is up to the user
#' to specify which rdf file contains each slot. In a general case, the user
#' specifies the `slot` that is found in a specific rdf file (`file`). A 
#' `summary` function is applied to a subset `period` of the `slot`, and then
#' compared (`eval`) to a threshold (`t_s`) and saved as the `variable`.
#' 
#' - *file:* specifies the rdf file that contains the slot.
#' 
#' - *slot:* the full RiverWare slot name, i.e., "Object.Slot".
#' 
#' - *period:* the period that the slot should be summarized over. This should 
#'   either be a function name, a full month name (found in [month.name]), or 
#'   the keyword "asis". 
#'   
#'     - *function name:* Specifying a function name allows for pre-specified 
#'       or custom functions to group together several months in the *period*. 
#'       This package provides the following functions: `cy()`, `wy()`, 
#'       `eocy()`, and `eowy()`. `cy()` indicates the data will be summarized 
#'       over the calendar year, i.e., January - December, while `wy()` 
#'       summarizes over the water year, i.e., October - September. `eocy()` 
#'       selects the end of the calendar year, and `eowy()` selects the end of 
#'       the water year. When specified in the `slot_agg` object, leave off the 
#'       parenthesis, i.e., only specify the function name. If `wy()` is 
#'       specified, the function will remove data for any water years that have 
#'       less than 7 months of data. This "tolerance" is specified by the 
#'       `rwdataplyr.wy_month_tol` option, and can be modified by updating this 
#'       option to another number. For standard monthly data that starts in 
#'       January and ends in December, this results in keeping the first water 
#'       year, since it includes 9 months of data, and removing the final water 
#'       year, since it includes only three months of data. Setting this option 
#'       to 0 will result in keeping any water year data that has at least one 
#'       month of data; setting this option to 11, ensures that there must be 
#'       a full water year of data for that year to be kept.
#'   
#'       This can also be a user specified custom function; see the
#'       *Custom Period Functions* section for details on constructing the custom
#'       functions.
#'     
#'     - *full month name:* When the full month name is specified, data will 
#'       be filtered to only include data for that particular month. To select
#'       multiple months of data, use a function as described above. If the
#'       month specified is not found in [month.name], an error will occur.
#'     
#'     - *asis:* If the keyword "asis" is specified, the data is returned for 
#'       its native timestep, i.e, monthly data will return monthly data and 
#'       annual data will return annual.
#'     
#' - *summary:* the summary function that should be applied to the period 
#'   specified as a function name, or `NA`. If the `period` specified is "asis"
#'   or returns only one month, e.g., `eocy()`, then the summary should be `NA`.
#'   The summary function should only return one value; for that reason, most
#'   of the `Summary` [S4groupGeneric]s work. Notably, `range()` will not
#'   since it returns two values. There is no reason that a custom function
#'   will not work here, but it has not been tested. 
#'   
#' - *eval:* the comparison operator to use (see the `Compare` 
#'   [S4groupGeneric]s). If no comparison is desired, then `NA` should be used. 
#'   If `eval` is specified the value returned from applying the `summary` to
#'   the `period` will be compared to the threshold specified by `t_s`. The 
#'   results of the comparison are returned as 0 and 1 instead of `TRUE` and 
#'   `FALSE`.
#'   
#' - *t_s:* either the threshold to be compared to if `eval` is not `NA` or a 
#'   value to scale the result by, e.g,. 0.001 to convert from acre-ft to 
#'   thousand acre-ft. `NA` can also be specified to not scale the data.
#'   
#' - *variable:* the variable name that will be used to identify the results
#'   of applying the period, summary, comparison/scaling to. All variable names
#'   should be unique.
#' 
#' For example, to determine if the minimum water year elevation at Lake Powell
#' is below elevation 3550 feet, the following would be specified:
#' ```
#' data.frame(
#'   file = "KeySlots.rdf",
#'   slot = "Powell.Pool Elevation",
#'   period = "wy",
#'   summary = "min",
#'   eval = "<",
#'   t_s = 3550,
#'   variable = "powellLt3550",
#'   stringsAsFactors = FALSE
#' )
#' ```
#'   
#' @section Custom Period Functions:
#' 
#' Users can specify custom period functions to make it easier to group months
#' together in custom ways. For example a function could return all of the 
#' summer months, or the more complicated case groups months across different
#' calendar years together. In fact, `wy()` is an example of a function that 
#' does this; another example might be grouping December - February together 
#' for winter months. 
#' 
#' The custom period function should return a list with three elements:
#' - `fun` - a function that will modify a rwtbl and properly determine the 
#'   new `Year`s based on the custom period.
#' - `filter_months` - the months that should be grouped together.
#' - `group_tbl` - how to group the returned rwtbl; likely either `c("Year")` or
#'   `c("Year", "Month")`
#'   
#' See the "RWDataPlyr Workflow" vignette for example implementations of both 
#' the summer and winter custom functions described above.
#' 
#' @param x A data.frame with required column names and valid entries; see 
#'   *Details* and *Aggregation Specification* sections.
#' @param rdfs A vector of rdf names; see *Details* and 
#'   *Aggregation Specification* sections.
#'   
#' @examples 
#' # determine if Powell's minimum water year elevation is < 3550'
#' rwd_agg(
#'   data.frame(
#'     file = "KeySlots.rdf",
#'     slot = "Powell.Pool Elevation",
#'     period = "wy",
#'     summary = "min",
#'     eval = "<",
#'     t_s = 3550,
#'     variable = "powellLt3550",
#'     stringsAsFactors = FALSE
#'   )
#' )
#' 
#' # get all the monthly slots in KeySlots.rdf
#' rwd_agg(rdfs = "KeySlots.rdf")
#' 
#' @seealso [rwd_agg_template()], [read_rwd_agg()]
#'
#' @export

rwd_agg <- function(x = NULL, rdfs = NULL)
{
  if (!missing(x) & !missing(rdfs))
    stop("When creating a `rwd_agg`, specify either `x` or `rdfs`, not both.")
 
  if (missing(x)){
    file_col <- rdfs
    slot_col <- rep("all", length(rdfs))
    x <- data.frame(
      "file" = file_col,
      "slot" = slot_col,
      "period" = NA,
      "summary" = NA,
      "eval" = NA,
      "t_s" = NA,
      "variable" = NA,
      stringsAsFactors = FALSE
    )
  }
  
  validate_rwd_agg(new_rwd_agg(x))
}

new_rwd_agg <- function(x)
{
  stopifnot(is.data.frame(x))
  structure(x, class = c("rwd_agg", "data.frame"))
}

validate_rwd_agg <- function(x)
{
  if (!is.data.frame(x))
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  
  cols <- c("file", "slot", "period", "summary", "eval", "t_s", "variable")
  if (ncol(x) != 7 || !all(colnames(x) == cols)) {
    stop(
      "The `colnames(x)` must be exactly: ", paste(cols, collapse = ", "), 
      call. = FALSE
    )
  }
  
  # all columns should not be factors
  if ("factor" %in% simplify2array(lapply(1:7, function(cc) class(x[[cc]])))) {
    stop("No columns should be factors.", call. = FALSE)
  }
  
  # check valid file extensions (for now only rdfs)
  if (!all(tools::file_ext(x$file) %in% rwd_ext)) {
    stop(
      "All `file` extensions should be: ", paste(rwd_ext, collapse = ","),
      call. = FALSE
    )
  }
  
  # all variables should be unique if the "all" keyword isn't used
  tmp <- x[x$slot != "all",]
  if (length(unique(tmp$variable)) != nrow(tmp)) {
    stop("All `variable`s should be unique.", call. = FALSE)
  }
  
  # if period is "asis" or "eocy", then summary must be NA
  check_period_asis_eocy(x)
  
  # checks: if CY or WY, summary NA probably doesn't make sense
  check_period_wy_cy(x)
  
  # check the eval and t_s columns
  lapply(seq_len(nrow(x)), function(rr) check_eval_and_t_s(x[rr,]))
  
  x
}

#' Returns all of the file extensions that rwd_agg and the aggregate functions
#' can handle. Will eventually include .csv and .nc.
#' @noRd

rwd_ext <- c("rdf")

#' Test if the object is a rwd_agg
#'
#' @param x An object
#' 
#' @return `TRUE` if the object inherits from the `rwd_agg` class.
#' @export
is_rwd_agg <- function(x)
{
  inherits(x, "rwd_agg")
}

#' @rdname is_rwd_agg
#' @export
is.rwd_agg <- is_rwd_agg

#' Given a rdf_file name and the rwtbl of that rdf, create a `rwd_agg` object
#' that contains all of the slots found in that rdf.
#' @noRd
rwd_agg_build_all <- function(rwtbl, rdf_file)
{
  slots <- as.character(levels(as.factor(rwtbl$ObjectSlot)))
  slot_names <- obj_slot_name(slots)
  
  rwd_agg(data.frame(
    file = rdf_file,
    slot = slots,
    period = "asis",
    summary = NA,
    eval = NA,
    t_s = NA,
    variable = slot_names,
    stringsAsFactors = FALSE
  ))
}

#' Take an Object.Slot and construct a valid variable name. In general, this 
#' will replace "." with "_" and convert to all lowercase. It will also shorten
#' commonly used slot names. Right now, "Pool Elevation" is the only shortened
#' slot name, and it is set to pe.
#' @noRd
obj_slot_name <- function(slots)
{
  obj_slot <- simplify2array(strsplit(slots, ".", fixed = TRUE))
  obj_slot[2,obj_slot[2,] == "Pool Elevation"] <- "pe"
  # replace : and " " with _
  obj_slot <- gsub(":", "_", obj_slot)
  obj_slot <- gsub(" ", "_", obj_slot)
  # convert to all lower case
  obj_slot <- tolower(obj_slot)
  
  paste(obj_slot[1,], obj_slot[2,], sep = "_")
}
