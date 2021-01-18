forcelayout <- function(schedule, webinteract = TRUE, ttime = 0){
  # if(!is_tibble(schedule)) {
  #   schedule <- tibble(schedul)eR
  # }
  if(webinteract == TRUE){
    rcpp_forcelayout(schedule = schedule, path.package("DataViz"))
  }
  else
  {
    r_forcelayout(schedule = schedule, ttime = ttime)
  }
  return(invisible())
}


