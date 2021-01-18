
# @param agmeenpgi_dat a \code{list} equivalent to the output of \code{synth_data_inc}
# @param race_vec A vector containing counts of the total population 
# as well as counts disagregating by gender and race.
# should be: \code{unlist(<synth_data>$estimates$pop_by_race[<row i>,])}
synth_data_race <- function(agmeenpgi_dat, race_vec) {
  dat <- agmeenpgi_dat[[1]]
  race_vec <- race_vec[!grepl("other", names(race_vec))] # remove other as category
  
  # 2. create age buckets on which to condition
  ag_list <- split(dat, dat$gender)
  
  # 3. Apply race
  race_levels <- c("black, afr amer", "native amer", "asian", "pacific islander", "two or more races", 
                   "white alone", "hispanic, latino")
  
  race_vec_m <- race_vec[grepl("m_", names(race_vec))]
  race_vec_f <- race_vec[grepl("f_", names(race_vec))]
  
  ag_list[[1]] <- race_lapply(ag_list[[1]], v= race_vec_m, levels= race_levels)
  ag_list[[2]] <- race_lapply(ag_list[[2]], v= race_vec_f, levels= race_levels)
  ag_list <- do.call("rbind", ag_list)
  
  # final cleaning, etc, return
  ag_list <- factor_return(ag_list, prob_name= "p")
  if (!is.micro_synthetic(ag_list)) class(ag_list) <- c(class(ag_list), "micro_synthetic")
  return(ag_list)
}


# helper function. 
race_lapply <- function(l, v, levels) {
  if (sum(v) > 0) comp <- v / sum(v)
  
  st <- data.frame(pct= comp, levels= factor(levels, levels= levels))
  st <- base::split(st, 1:nrow(st))
  
  dat <- replicate(length(levels), l, simplify = FALSE)
  dat <- do.call("rbind", mapply(add_synth_attr_level, dat= dat, prob_name= "p", attr= st,
                                 attr_name= "race", SIMPLIFY = FALSE))
  return(dat)
}

