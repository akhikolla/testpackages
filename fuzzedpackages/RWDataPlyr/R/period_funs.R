
eocy <- function() "December"

eowy <- function() "September"

cy <- function() month.name

wy <- function()
{
  wy_convert <- function(rwtbl)
  {
    tmp <- rwtbl %>%
      dplyr::mutate_at(
        "Timestep", 
        .funs = list("Year" = zoo::as.yearmon)
      ) %>%
      dplyr::mutate_at("Year", .funs = list(ym_get_wateryear))
    
    # drop if WY contains less than 6 months of data for the year
    cols <- names(tmp)
    cols <- cols[!(cols %in% c("Timestep", "Month", "Value"))]
    keep_yrs <- unique(
      (tmp %>% 
         dplyr::group_by_at(cols) %>% 
         dplyr::tally() %>% 
         dplyr::filter_at(
           "n", 
           dplyr::any_vars(. > getOption("rwdataplyr.wy_month_tol"))
         )
      )$Year
    )
    
    dplyr::filter_at(tmp, "Year", dplyr::any_vars(. %in% keep_yrs))
  }
  
  list(fun = wy_convert, filter_months = month.name, group_tbl = c("Year"))
}

is_custom_period_fun <- function(x)
{

  is.list(x) && length(x) == 3 && 
    all(names(x) %in% c("fun", "filter_months", "group_tbl")) &&
    is.function(x$fun) && all(x$filter_months %in% month.name)
}
