
#' @details 
#' `rw_scen_aggregate()` aggregates multiple scenarios of data. It processes the 
#' [rwd_agg] object (`agg`) for each single scenario, and then binds all of the
#' individual scenario data together into a single `tbl_df`.
#' 
#' @inheritParams rdf_aggregate
#' @param scenarios A character vector of scenario folders. This is usually a 
#'   vector of folder names, where each folder name contains one scenario worth
#'   of data. `scenarios` can be named or unnamed. The names are used as the 
#'   scenario name in the returned `tbl_df`. Scenario names can also be 
#'   specified through the `scen_names` argument. If `scen_names` is specified, 
#'   `scenarios` should not already have names. If `scen_names` is not specified
#'   and, `scenarios` is not already named, then the scenario folders will also
#'   be used as the scenario names. See **Directory Structure**.
#' @param scen_dir File path to the directory that contains the scenario 
#'   folders. **Directory Structure**.
#' @param file Optionally save the `tbl_df` of aggregated scenario data as a 
#'   .txt, .csv, or .feather file. If `file` is specified, then the data are 
#'   saved in the specified output format.
#' @param scen_names An alternative way to specify scenario names. 
#' 
#' @section Directory Structure:
#' 
#' RiverWare and RiverSMART typically write data into an expected directory 
#' structure. The below shows an example directory structure and corresponding
#' variable names for `rw_scen_aggregate()` and `rdf_aggregate()`. (In the
#' example below, C:/user/crss/CRSS.Jan2017/Scenario is the more complete 
#' directory setup for the data included in `system.file("extdata/Scenario/")`.)
#' 
#' \preformatted{
#' C:/user/crss
#' |
#' |- CRSS.Jan2017
#' |    - model
#' |    - ruleset
#' |    - Scenario
#' |         - ISM1988_2014,2007Dems,IG,Most
#' |         - ISM1988_2014,2007Dems,IG,2002 
#' |    - ...
#' |- CRSS.Jan2018
#' |    - model
#' |    - ... (same general setup as CRSS.Jan2017)
#' }
#' 
#' To get one scenario's data, `rdf_aggregate()` can be called with `rdf_dir`
#' set to "C:/user/crss/CRSS.Jan2017/Scenario/ISM1988_2014,2007Dems,IG,Most".
#' (`scenario` can optionally be specified to git a scenario name.)
#' 
#' To aggregate multiple scenarios of data together, `rw_scen_aggregate()` 
#' should be called with `scen_dir` set to "C:/user/CRSS/CRSS.Jan2017/Scenario" 
#' and `scenarios` set to 
#' `c("ISM1988_2014,2007Dems,IG,Most", "ISM1988_2014,2007Dems,IG,2002")`. 
#' (Optionally, `scenarios` can be named, or `scen_names` specified to use 
#' scenario names that are different from the above scenario folders.)
#' 
#' Finally, to aggregate scenario data from both CRSS.Jan2017 and CRSS.Jan2018,
#' `rw_scen_aggregate()` should be called with `scen_dir` set to 
#' "C:/users/crss/". `scenarios` can then be set to 
#' `c("CRSS.Jan2017/Scenario/ISM1988_2014,2007Dems,IG,Most","CRSS.Jan2018/Scenario/ISM1988_2014,2007Dems,IG,Most")`,
#' assuming the same scenario exists in both folders. In this case it is 
#' advisable to also specify `scen_names` or name `scenarios`.
#' 
#' @return A `tbl_df` containing all aggregated and summarized data for all of
#'   the specified `scenarios`.
#'   
#' @examples 
#' # rw_scen_aggregate() ----------
#' 
#' scens <- c("ISM1988_2014,2007Dems,IG,2002", "ISM1988_2014,2007Dems,IG,Most")
#' scenNames <- c("2002", "Most")
#' namedScens <- scens
#' names(namedScens) <- scenNames
#' 
#' scenPath <- system.file("extdata/Scenario", package = "RWDataPlyr")
#' 
#' rwa <- read_rwd_agg(
#'   system.file(
#'     "extdata/rwd_agg_files/passing_aggs.csv", 
#'     package = "RWDataPlyr"
#'   )
#' )
#' 
#' x <- rw_scen_aggregate(namedScens, agg = rwa[1,], scen_dir = scenPath)
#' 
#' # y will be identical to x
#' 
#' y <- rw_scen_aggregate(
#'   scens, 
#'   agg = rwa[1,], 
#'   scen_dir = scenPath, 
#'   scen_names = scenNames
#' )
#' 
#' identical(x, y) # is TRUE
#' 
#' @rdname rdf_aggregate
#' 
#' @export
rw_scen_aggregate <- function(scenarios, 
                              agg, 
                              scen_dir = ".", 
                              nans_are = "0",
                              keep_cols = FALSE,
                              file = NULL, 
                              scen_names = NULL,
                              find_all_slots = TRUE,
                              cpp = TRUE,
                              verbose = TRUE)
{
  # check all UI --------------------
  if (!is.rwd_agg(agg)) {
    stop("In `rw_scen_aggregate()`, `agg` must be a `rwd_agg` object.")
  }
  
  if (length(scen_dir) != 1) {
    stop("`scen_dir` should only have a length of 1.", call. = FALSE)
  }
  scen_dir <- normalizePath(scen_dir, mustWork = TRUE)
  
  nans_are <- match.arg(nans_are, choices = c("0", "error"))
  
  check_scen_rdf_paths(scenarios, scen_dir, agg)
  
  # check that file is correct type if it's specified
  if (!missing(file)) {
    check_rw_agg_file(file)
  }

  scenarios <- get_scen_names(scenarios, scen_names)
  
  if (verbose) rw_scen_agg_msg(scenarios)
  # aggregate all scenarios -------------
  
  nScen <- seq_len(length(scenarios))
  
  rwscenagg <- lapply(nScen, function(x) {
    if (verbose) rw_scen_agg_msg(scenarios, x)
    
    rdf_aggregate(
      agg, 
      rdf_dir = file.path(scen_dir, scenarios[x]),
      scenario = names(scenarios)[x],
      keep_cols = keep_cols,
      nans_are = nans_are,
      find_all_slots = find_all_slots,
      cpp = cpp,
      verbose = verbose
    )
  })
  rwd_agg <- attr(rwscenagg[[1]], "rwd_agg")
  rdf_atts <- lapply(nScen, function(x) attr(rwscenagg[[x]], "rdf_atts"))
  names(rdf_atts) <- names(scenarios)
  scen_folder <- lapply(nScen, function(x) attr(rwscenagg[[x]], "scen_folder"))
  scen_folder <- dplyr::bind_rows(scen_folder)
    
  rwscenagg <- dplyr::bind_rows(rwscenagg)
  
  if (!missing(file)) {
    write_rw_data(rwscenagg, file)
  }
  
  structure(
    rwscenagg,
    "rwd_agg" = rwd_agg,
    "rdf_atts" = rdf_atts,
    "scen_folder" = scen_folder
  )
}

#' Since scenario names can be specified multiple ways, `get_scen_names()` will
#' check the different ways and error/return scenarios with proper names
#' @noRd
get_scen_names <- function(scenarios, scen_names)
{
  # if scenarios have names, then scen_names should not be specified
  if (!is.null(names(scenarios)) && !is.null(scen_names)) {
    stop(
      "In `rw_scen_aggregate()`, `scenarios` have `names()`, so `scen_names` should not be specified.",
      call. = FALSE
    )
  }
  
  if (is.null(names(scenarios)) && !is.null(scen_names)) {
    if (length(scenarios) != length(scen_names)) {
      stop(
        "In `rw_scen_aggregate()`, `scenarios` and `scen_names` must have the same length.",
        call. = FALSE
      )
    }
    names(scenarios) <- scen_names
  }
  
  if (!is.null(names(scenarios))) {
    if (any(names(scenarios) %in% "")) {
      repI <- which(names(scenarios) %in% "")
      names(scenarios)[repI] <- scenarios[repI]
    }
  } else if (is.null(names(scenarios)) && is.null(scen_names)) {
    names(scenarios) <- scenarios
  }
  
  scenarios
}

#' Check that all scenario paths, and rdf files in the scenario paths can be
#' found, prior to processing data.
#' @noRd
check_scen_rdf_paths <- function(scenarios, scen_dir, agg)
{
  # check all scenario paths
  scen_paths <- normalizePath(file.path(scen_dir, scenarios), mustWork = FALSE)
  if (! all(dir.exists(scen_paths))) {
    stop(
      "The following scenario directories do not exist:\n",
      paste(scen_paths[!dir.exists(scen_paths)], collapse = "\n"),
      call. = FALSE
    )
  }
  
  # check that all rdfs exist
  rdfs <- expand.grid(scen_paths, unique(agg$file))
  rdfs <- normalizePath(apply(rdfs, 1, paste, collapse = "/"), mustWork = FALSE)
  rdfs_exist <- file.exists(rdfs)
  if (!all(rdfs_exist)) {
    stop(
      "The following rdf files do not exist:\n",
      paste(rdfs[!rdfs_exist], collapse = "\n"),
      call. = FALSE
    )
  }
}

#' For the output file, verify it is correctly specified.
#' @noRd
check_rw_agg_file <- function(file)
{
  if (length(file) != 1) {
    stop("In `rw_scen_aggregate()`, `file` should have a length of 1.")
  }
  
  if (!(tools::file_ext(file) %in% c("csv", "feather", "txt"))) {
    stop(
      "In `rw_scen_aggregate()`, `file` should have a .csv, .feather, or .txt extension.", 
      call. = FALSE
    )
  }
  
  if (!dir.exists(dirname(file))) {
    stop("In `rw_scen_aggregate()`, `file` should point to a valid location.")
  }
  
  invisible(file)
}

rw_scen_agg_msg <- function(scenarios, n = NA_integer_)
{
  if (is.na(n)) {
    scen <- if (length(scenarios) == 1) "scenario" else "scenarios"
    message("Processing ", length(scenarios), " total ", scen, ".\n",
            "------------------")
  } else {
    message("** Starting scenario ", n, " of ", length(scenarios), ": ", 
            scenarios[n])
  }
  
  invisible(scenarios)
}
