
#' Given user specified `period` input found in `slot_agg_row`, filter the 
#' `Month` in `rwtbl` to only include those specified by `period`. 
#' @noRd

apply_period <- function(rwtbl, slot_agg_row)
{
  # check that it has Year and Month columns
  
  # filter based on slot
  rwtbl <- dplyr::filter_at(
    rwtbl, 
    "ObjectSlot", 
    dplyr::any_vars(. == slot_agg_row$slot)
  ) %>%
    # filter and group for period
    period_filter_group(period = slot_agg_row$period)
  
  rwtbl
}

#' Call the function specified by `period`, and filter `Month` based on it
#' 
#' Also group by the pre-specified or user specified grouping.
#' 
#' @noRd

period_filter_group <- function(rwtbl, period)
{
  if (!(period %in% c(month.name, "asis"))) {
    # convert to appropriate function, and then chuck that function exists
    if (!exists(period, mode = "function"))
        stop("specified `period`: ", period, 
             " does not match expected values or existing functions.\n",
             "   Please see ?rwd_agg for help.", call. = FALSE)
    
    period_filter <- tryCatch(
      eval(parse(text = paste0(period, "()"))), 
      error = function(c) -1
    )
    
    check_period_filter(period_filter, period)
    
  } else if (period %in% month.name) {
    period_filter <- period
  }
  
  if (period == "asis") {
    rwtbl <- rwtbl %>%
      dplyr::group_by_at(c("Year", "Month"))
  } else if (is_custom_period_fun(period_filter)) {
    # list(fun = wy_convert, filter_months = month.name, group_tbl = c("Year"))
    rwtbl <- period_filter$fun(rwtbl) %>%
      dplyr::filter_at(
        "Month", 
        dplyr::any_vars(. %in% period_filter$filter_months)
      ) %>%
      dplyr::group_by_at(period_filter$group_tbl)
  } else {
    rwtbl <- rwtbl %>%
      dplyr::filter_at("Month", dplyr::any_vars(. %in% period_filter)) %>%
      dplyr::group_by_at("Year")
  }
  
  rwtbl
}

#' After evaluating `period` (`period_filter`), ensure it meets all requirements
#' @noRd

check_period_filter <- function(period_filter, period)
{
  ftxt <- paste0("Function `", period, "()` ")
  
  if (is.numeric(period_filter) && period_filter == -1)
    stop(ftxt, "exists, but could not be evaluated.\n",
         "  It should not require any arguments.", call. = FALSE)
  
  if (is.character(period_filter) && !all(period_filter %in% month.name))
    stop(ftxt, "is returning characters that are not found in `month.name`",
         call. = FALSE)
  
  if (!is.character(period_filter) && !is_custom_period_fun(period_filter)) {
    stop(ftxt, "must return either a function or character vector.", 
         call. = FALSE)
  }
}

#' If the period is "asis" or eocy, then summary should be NA
#' @noRd

check_period_asis_eocy <- function(rwd_agg)
{
  r2 <- rwd_agg[rwd_agg$period == "asis",]
  if (nrow(r2) > 0 && nrow(r2[!is.na(r2$summary),]) > 0) {
    stop(
      "If the `period` is specified as 'asis', then the `summary` must be `NA`.",
      call. = FALSE
    )
  }
  
  r2 <- rwd_agg[rwd_agg$period == "eocy",]
  if (nrow(r2) > 0 && nrow(r2[!is.na(r2$summary),]) > 0) {
    stop(
      "If the `period` is specified as 'eocy', then the `summary` must be `NA`.",
      call. = FALSE
    )
  }
  invisible(rwd_agg)
}

#' check to see if the period is wy or cy; if it is, summary should not be NA
#' @noRd

check_period_wy_cy <- function(rwd_agg)
{
  r2 <- rwd_agg[rwd_agg$period %in% c("cy", "wy"),]
  if (nrow(r2) > 0 && nrow(r2[is.na(r2$summary),]) > 0) {
    stop(
      "If the `period` is 'CY', or 'WY', then the `summary` should not be `NA`.",
      call. = FALSE
    )
  }
  
  invisible(rwd_agg)
}
