
#' Given the `eval` and `t_s` columns, evaluate the comparison threshold or
#' multiply by the scalilng factor if now comparison threshold exists. If 
#' using a comparison threshold, will return either 0 or 1 instead of FALSE or
#' TRUE.
#' 
#' @noRd

apply_eval <- function(rwtbl, slot_agg_row)
{
  check_eval_and_t_s(slot_agg_row)
  
  eval_fun <- slot_agg_row$eval
  t_s <- slot_agg_row$t_s
  
  if (is.na(eval_fun) & !is.na(t_s))
    eval_fun <- "*"
  
  if (!is.na(eval_fun) & !is.na(t_s)){
    rwtbl <- rwtbl %>%
      dplyr::mutate_at(
        "Value", 
        list(~as.numeric(eval(parse(text = paste(., eval_fun, t_s)))))
      )
  }
  
  rwtbl
}

#' Ensure that `eval` and `t_s` columns in a slot agg matrix are valid
#'
#' @noRd

check_eval_and_t_s <- function(slot_agg_row)
{
  # eval column should be NA or "<", "<=", ">", ">=", "!=", "=="
  
  eval_col <- slot_agg_row$eval
  
  if (!is.na(eval_col) & !(eval_col %in% methods::getGroupMembers("Compare"))) {
    stop(
      "'", eval_col, "' is not a valid `eval` value.\n",
      "The `eval` column in the slot agg matrix should either be\n",
      "`NA` or one of the 'Compare' S4 group generics. See ?S4groupGeneric.", 
      call. = FALSE
    )
  }
  
  # if eval is NA, then t_s can either be na or a numeric; 
  # if eval is a Compare generic, then t_s must be a numeric
  
  t_s <- slot_agg_row$t_s
  
  if (!is.na(t_s)) {
    t_s <- tryCatch(
      as.numeric(t_s), 
      error = function(c) NaN, 
      warning = function(c) NaN
    )
  }

  if (is.nan(t_s))
    stop(
      "'", slot_agg_row$t_s, "' is not a valid `t_s` value.\n",
      "The `t_s` column in the slot agg matrix should either be\n",
      "`NA` or a numerical value.", 
      call. = FALSE
    )
  
  if (!is.na(eval_col) & is.na(t_s))
    stop("When the `eval` column is a comparison function, the `t_s` column\n",
         "must be a numerical value.", call. = FALSE)
  
  invisible(slot_agg_row)
}
