
#' Given user specified `summary` input found in `slot_agg_row`, summarise the 
#' `rwtbl`, dropping unused columns `Timestep` and possibly `Month`. If no 
#' summary is specified (`NA`), then the function groups and drops the columns
#' as summary would, but does not modify `Value`.
#' @noRd

apply_summary <- function(rwtbl, slot_agg_row)
{
  if (length(dplyr::groups(rwtbl)) == 0 || is.null(dplyr::groups(rwtbl)))
    stop("rwtbl should already have groups when `apply_summary()` is called")
  
  # slot_agg_row should either be NA, or a string for an existing "summary" 
  # type function
  
  # should group by all columns, except Timestep and Value;
  # only group by month if it is already a grouping variable
  cur_groups <- dplyr::group_vars(rwtbl)
  cols <- colnames(rwtbl)
  cols <- cols[!(cols %in% c("Timestep", "Year", "Month", "Value"))]
  rwtbl <- dplyr::group_by_at(rwtbl, c(cur_groups, cols))
  
  if (!is.na(slot_agg_row$summary)){
    rwtbl <- summary_summarise(rwtbl, slot_agg_row$summary)
  } else {
    # drop the columns that aren't returned if you are summarising the tbl
    drop_cols <- colnames(rwtbl)[!(
      colnames(rwtbl) %in% c(dplyr::group_vars(rwtbl), "Value")
    )]
    
    rwtbl <- dplyr::select(rwtbl, -dplyr::one_of(drop_cols)) %>%
      # drop the last grouping variable to match output if you do summarise
      dplyr::group_by_at(c(cur_groups, utils::head(cols, -1)))
  }
  
  rwtbl
}

#' Checks that the summary function `sam_summary` exists and meets other 
#' requirements, e.g., returns only one value and accepts only one vector as
#' its arguement. Then summarises `rwtbl`.
#' @noRd

summary_summarise <- function(rwtbl, sam_summary)
{
  if (!exists(sam_summary, mode = "function"))
    stop("specified `summary`: ", sam_summary, 
         " does not match existing functions.\n",
         "   Please see ?rwd_agg for help.", call. = FALSE)
  
  smry_fun <- tryCatch(
    eval(parse(text = sam_summary)), 
    error = function(c) -1
  )

  check_summary_function(smry_fun, sam_summary)
  
  dplyr::summarise_at(rwtbl, "Value", .funs = list(~smry_fun(.)))
}

#' Checks that the summary function `sam_summary` meets other 
#' requirements, e.g., returns only one value and accepts only one vector as
#' its arguement.
#' 
#' @noRd

check_summary_function <- function(smry_fun, sam_summary)
{
  ftxt <- paste0("`", sam_summary, "()` ")
  
  if (is.numeric(smry_fun) && smry_fun == -1)
    stop(ftxt, "exists, but could not be evaluated.", call. = FALSE)
  
  s1 <- tryCatch(smry_fun(5), error = function(c) "error")
  s2 <- tryCatch(smry_fun(0:12), error = function(c) "error")
  s3 <- tryCatch(smry_fun(-24), error = function(c) "error")
  s4 <- tryCatch(smry_fun(-34:1), error = function(c) "error")
  sL <- list(s1, s2, s3, s4)
  
  areErrors <- simplify2array(lapply(1:4, function(x) {
    (is.character(sL[[x]]) && sL[[x]] == "error")
  }))
  
  if (any(areErrors))
    stop(ftxt, "resulted in an error for the simple test cases.\n",
         "Ensure that it only requires one vector as its arguement.")
  
  areErrors <- simplify2array(lapply(1:4, function(x) length(sL[[x]]) != 1))
  
  if (any(areErrors))
    stop(ftxt, "returns more than 1 value for a vector")
  
  invisible(smry_fun)
}
