
#' Convert an rdf to a tibble
#' 
#' `rdf_to_rwtbl()` converts an rdf list to a tibble.
#' 
#' The rdf object is converted to a data frame, and then converted to a 
#' [tibble::tibble()]. All of the `meta` entries in the rdf object
#' are stored as attributes in the returned tibble. These attributes are:
#' `mrm_config_name`, `owner`, `description`, `create_date`, and `n_traces`.
#' 
#' If the rdf contains a scalar slot(s), the scalar slot value(s) will be 
#' repeated for every timestep.
#' 
#' @param rdf An rdf object (from [read_rdf()]).
#' @param scenario An optional parameter, that if it is not `NULL` or `NA` 
#'   (default) will be added to the tibble as another variable. Coerced to a 
#'   character if it is not already a character.
#' @param keep_cols Either boolean, or a character vector of column names to 
#'   keep in the returned tibble. The values of `keep_cols` work as follows:
#'   * `FALSE` (default) only includes the defaults columns: `Timestep`, 
#'   `TraceNumber`, `ObjectSlot`, and `Value`. `Scenario` is also returned if 
#'   `scenario` is specified.
#'   * `TRUE`, all columns are returned.
#'   * A character vector, e.g., `c("ObjectName", "Units")`, allows the user to 
#'   include other columns that are not always required, in addition to the 
#'   "default" set of columns. If any of the values in `keep_cols` are not 
#'   found, a warning will post, but all other columns will be returned. 
#' @param add_ym Boolean that controls whether or not `Year` and `Month` columns
#'   are added to the returned tibble. If `TRUE` (default), they will be added, 
#'   and if `FALSE` they will not be added. They are constructed from the dates 
#'   in the `Timestep` column.
#' 
#' @return A tbl_df with additional attributes from the rdf object.
#' 
#' @examples 
#' rdftbl <- rdf_to_rwtbl(keyRdf)
#' # same as previous, except you do not want "Year" and "Month" columns
#' rdftbl <- rdf_to_rwtbl(keyRdf, add_ym = FALSE)
#' # but you do want to keep the object name seperately:
#' rdftbl <- rdf_to_rwtbl(keyRdf, add_ym = FALSE, keep_cols = "Object")
#' rdftbl <- rdf_to_rwtbl(sysRdf, scenario = "ISM1988_2014,2007Dems,IG,2002")
#' 
#' # rdf_to_rwtbl2 wants a file path instead of an rdf object
#' rdfPath <- system.file(
#'   "extdata/Scenario/ISM1988_2014,2007Dems,IG,Most/KeySlots.rdf", 
#'   package = "RWDataPlyr"
#' )
#' rdftbl <- rdf_to_rwtbl2(rdfPath)
#' 
#' @seealso [read_rdf()]
#' 
#' @export

rdf_to_rwtbl <- function(rdf, scenario = NULL, keep_cols = FALSE, add_ym = TRUE)
{
  .Deprecated(
    msg = paste0(
      "`rdf_to_rwtbl()` is deprecated and will be removed in a future release.\n",
      "Consider using `rdf_to_rwtbl2()` instead.\n",
      "See help(\"rdf_to_rwtbl()\"), as it uses different arguments."
    )
  )
  stopifnot(is_rdf(rdf))
  
  check_rdf_to_rwtbl_args(scenario, keep_cols, add_ym, "rdf_to_rwtbl")
  
  # rdf[["meta"]] contains meta data
  # rdf[["runs"]] will be the length of rdf[[1]]$number_of_runs
  atts <- rdf$meta
  nRun <- as.integer(rdf$meta$number_of_runs)
  
  tbl <- do.call(
    rbind, 
    lapply(1:nRun, function(x) rdf_trace_to_tbl(rdf[[2]][[x]], traceNum = x))
  )
  
  if ((is.logical(keep_cols) && !keep_cols) || is.character(keep_cols))
    tbl <- select_rdftbl_cols(tbl, keep_cols)
  
  if (!is.null(scenario)) {
    if (!is.character(scenario))
      scenario <- as.character(scenario)
    tbl$Scenario <- scenario
  }
  
  if (add_ym)
    tbl <- add_ym_to_rdftbl(tbl)
    
  structure(
    tibble::as_tibble(tbl),
    "mrm_config_name" = atts$name,
    "owner" = atts$owner,
    "description" = atts$description,
    "create_date" = atts$create_date,
    "n_traces" = nRun
  )
}

#' `rdf_to_rwtbl2()` converts an rdf file into a tibble, but is 
#' faster than `rdf_to_rwtbl()` since it uses c++. It also reads the rdf file
#' in, while `rdf_to_rwtbl()` needs an `rdf` object.
#' 
#' @param file The relative or absolute rdf filename.
#'  
#' @rdname rdf_to_rwtbl
#' 
#' @export
rdf_to_rwtbl2 <- function(file, scenario = NA_character_, keep_cols = FALSE, 
                          add_ym = TRUE)
{
  if (! is.character(file) & length(file) != 1) {
    stop(
      "In `rdf_to_rwtbl2()`, `file` should be a single entry character vector,",
      "  i.e., a file path to one rdf file.",
      call. = FALSE
    )
  }
  
  rdf_vec <- read_rdf(file, rdf = FALSE)
  
  check_rdf_to_rwtbl_args(scenario, keep_cols, add_ym, "rdf_to_rwtbl2")

  add_cols <- c("ObjectName", "SlotName", "ObjectType" ,"Unit", 
                "RulesetFileName", "InputDMIName")

  if (is.logical(keep_cols)) {
    if (keep_cols) {
      keep_cols <- c(req_rwtbl_cols(), add_cols)
    } else {
      keep_cols <- req_rwtbl_cols()
    }
  } else {
    # combine keep_cols with the required columns
    keep_cols <- keep_cols[!(keep_cols %in% req_rwtbl_cols())]
    keep_cols <- c(req_rwtbl_cols(), keep_cols)
  }
  
  # check that keep_cols are all expected column names
  missing_cols <- keep_cols[!(keep_cols %in% c(req_rwtbl_cols(), add_cols))]
  if (length(missing_cols) > 0)
    warning(
      "The following columns specified by 'keep_cols' were not found in the rwtbl:\n",
      "    ", toString(missing_cols)
    )
  
  if (is.null(scenario)) {
    scenario <- NA_character_
  }
  
  if (! is.character(scenario))
    scenario <- as.character(scenario)
  
  t1 <- rdf_to_rwtbl_cpp(
    rdf_vec, 
    scenario = scenario, 
    keep_cols = keep_cols,
    add_ym = add_ym
  ) %>% 
    tibble::as_tibble()
  
  t1
}


#' create a data frame from one of the objects in one trace of the rdf list
#' 
#' @noRd

rdf_object_to_tbl <- function(rdfObject, timeSteps)
{
  ot <- rdfObject$object_type
  on <- rdfObject$object_name
  sn <- rdfObject$slot_name
  os <- paste(on, sn, sep = ".")
  units <- rdfObject$units
  scale <- as.numeric(rdfObject$scale)
  value <- rdfObject$values * scale
  
  tbl <- data.frame("Timestep" = timeSteps) %>%
    dplyr::mutate(
      ObjectName = on,
      SlotName = sn,
      ObjectType = ot,
      ObjectSlot = os,
      Value = value,
      Unit = units
    )
  
  tbl
}


#' for a single run/trace, take the list and create a data frame from it

#' "loop" through all objects within 1 trace of data
#' @noRd

rdf_trace_to_tbl <- function(rdfTrace, traceNum) 
{
  trace <- as.integer(rdfTrace$idx_sequential)
  timeSteps <- rdfTrace$times
  ruleSet <- rdfTrace$rule_set
  slotSet <- rdfTrace$slot_set # the input DMI on the input tab
  
  # ***
  # maybe check if time_step_unit is year, and if so, change the timestep label
  # maybe change to yearmon here?
  # ***
  
  # for all objects, call getObjectData for rdfTrace$objects
  tbl <- do.call(
    rbind, 
    lapply(rdfTrace$objects, rdf_object_to_tbl, timeSteps = timeSteps)
  ) %>%
    dplyr::mutate(
      TraceNumber = traceNum, 
      RulesetFileName = ruleSet, 
      InputDMIName = slotSet
    ) %>%
    dplyr::mutate_at(.vars = "Timestep", .funs = as.character)
  
  tbl
}

#' Standard (required) Column Names for RDF Tables that will always be returned
#' by `rdf_to_rwtbl()`
#' 
#' @noRd

req_rwtbl_cols <- function()
{
  c("Timestep","TraceNumber", "ObjectSlot", "Value")
}

#' Select a subset of the columns in the rdftbl, as specified by `keep_cols`
#' 
#' @noRd

select_rdftbl_cols <- function(rdftbl, keep_cols)
{
  if (is.logical(keep_cols))
  {
    select_cols <- req_rwtbl_cols()
  } else {
    # check that the user specified keep_cols are actually columns in rdftbl
    colsInTbl <- (keep_cols %in% colnames(rdftbl))
    missing_cols <- keep_cols[!colsInTbl]
    if (length(missing_cols) > 0)
      warning(
        "The following columns specified by 'keep_cols' were not found in the rwtbl:\n",
        "    ", toString(missing_cols)
      )
    select_cols <- c(req_rwtbl_cols(), keep_cols[colsInTbl])
  }
  
  rdftbl %>%
    dplyr::select_at(.vars = select_cols)
}

#' Adds "Year", and "Month" columns to the tbl, which are constructed from 
#' the "Timestep" column.
#' @noRd

add_ym_to_rdftbl <- function(tbl)
{
  other_cols <- colnames(tbl)[colnames(tbl) != "Timestep"]
  tbl <- tbl %>%
    dplyr::mutate_at("Timestep", list("ym" = zoo::as.yearmon)) %>%
    dplyr::mutate_at(
      "ym", 
      list("Year" = ym_get_year, "Month" = ym_get_month_str)) %>%
    dplyr::select_at(c("Timestep", "Year", "Month", other_cols))
}

rwtbl_get_atts <- function(rwtbl)
{
  attributes(rwtbl)[match(
    c("mrm_config_name", "owner", "description", "create_date", "n_traces"), 
    names(attributes(rwtbl))
  )]
}


check_rdf_to_rwtbl_args <- function(scenario, keep_cols, add_ym, call_func)
{
  if (is.logical(keep_cols)) {
    if (length(keep_cols) != 1 || is.na(keep_cols)) {
      stop(
        "In `", call_func, "`(), if `keep_cols`, is logical, it should have ",
        "a length of 1, and only be TRUE or FALSE.",
        call. = FALSE
      )
    }
  } else if (is.character(keep_cols)) {
    if ( ! (length(keep_cols) > 0)) {
      stop(
        "In `", call_func, "`(), if `keep_cols`, is a character, it should ",
        "have a length of 1 or more.",
        call. = FALSE
      )
    }
  } else {
    stop(
      "In `", call_func, "`(), `keep_cols` should either be a character or ",
      "boolean.",
      call. = FALSE
    )
  }
  
  if (length(add_ym) != 1 || !is.logical(add_ym) || is.na(add_ym)) {
    stop(
      "In `", call_func, "`(), `add_ym` should be a length 1 boolean.",
      call. = FALSE
    )
  }
  
  # this function only reads in one scenario, so should only have one 
  # scenario name
  if (!is.null(scenario)) {
    if (length(scenario) != 1) {
      stop(
        "In `", call_func, "`(), `scenario`, should either be NULL, or have a ",
        "length of 1.",
        call. = FALSE
      )
    }
  }
}