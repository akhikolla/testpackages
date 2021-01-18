
# @param agmee_dat a \code{list} equivalent to the output of \code{synth_data_emp}
# @param nativity_vec A vector containing counts of the total population 
# (Universe:  Population 25 years and over in the United States)
# as well as counts disagregating geographic mobility by educational attainment.
# should be: \code{unlist(<synth_data>$estimates$geo_mob_edu[<row i>,])}
synth_data_geomob <- function(agmeenp_dat, geomob_vec) {
  dat <- agmeenp_dat[[1]]
  # remove total counts
  geomob_vec <- geomob_vec[!grepl("_all", names(geomob_vec))]
  
  # 1. create hash table for mapping educational attainment
  ht <- data.frame(prior_dat= agmeenp_dat[[2]],
                   new_dat= c(rep("lt_hs", 2), "high_sch", rep("some_col",2), "bachelors", "graduate"),
                   stringsAsFactors = FALSE)
  # 2. create age buckets on which to condition
  ag_list <- split(dat, dat$edu_attain)
  
  # 3. Apply nativity status
  geo_levels <- c("same house", "same county", "same state", "diff state", "moved from abroad")
  
  ag_list <- do.call("rbind", lapply(ag_list, geo_lapply, ht= ht, 
                                     v= geomob_vec, levels= geo_levels))
  
  ag_list <- factor_return(ag_list, prob_name= "p")
  return(list(ag_list))
}


# helper function 
geo_lapply <- function(l, ht, v, levels) {
  if (nrow(l) < 1) return(l)
  l_comp <- ht[,2][which(l$edu_attain[1] == ht[,1])]
  comp <- v[which(grepl(l_comp, names(v)))]
  if (sum(comp) > 0) comp <- (comp / sum(comp)) 
  
  st <- data.frame(pct= comp, levels= factor(levels, levels= levels))
  st <- base::split(st, 1:nrow(st))
  
  dat <- replicate(length(levels), l, simplify = FALSE)
  dat <- do.call("rbind", mapply(add_synth_attr_level, dat= dat, prob_name= "p", attr= st,
                                 attr_name= "geog_mobility", SIMPLIFY = FALSE))
  return(dat)
}

