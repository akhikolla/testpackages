

# @param agmee_dat a \code{list} equivalent to the output of \code{synth_data_emp}
# @param inc_nat_vec A vector containing counts of the total population 
# (Universe:  Population 15 years and over in the United States)
# as well as counts disagregating income level by nativity status.
# should be: \code{unlist(<synth_data>$estimates$by_inc_12_mo[<row i>,])}
synth_data_inc <- function(agmeenpg_dat, inc_nat_vec) {
  dat <- agmeenpg_dat[[1]]
  # remove total counts
  inc_nat_vec <- inc_nat_vec[!grepl("cnts_", names(inc_nat_vec))]
  
  # 1. create hash table for mapping educational attainment
  ht <- data.frame(prior_dat= c("born_other_state","born_out_us", "born_state_residence", "foreigner"),
                   new_dat= c("cit_born_other_st", "cit_born_out_us", "cit_born_st_res", "foreign_born"),
                   stringsAsFactors = FALSE)
  
  # 2. create buckets on which to condition
  ag_list <- split(dat, dat$pov_status)
  ag_list[[1]] <- split(ag_list[[1]], ag_list[[1]]$nativity) # those at_above_pov_level
  ag_list[[2]] <- split(ag_list[[2]], ag_list[[2]]$nativity) # those below_pov_level
  
  # subset counts for incomes below pov level (2014 ~= $12k)
  inc_levels_pov <- c("no_income", "1_lt10k", "10k_lt15k")
  inc_levels_nonpov <- c("10k_lt15k", "15k_lt25k", "25k_lt35k", "35k_lt50k", 
                         "50k_lt65k", "65k_lt75k", "gt75k")
  
  pov_inc_vec <- inc_nat_vec[grepl(paste(c("no_inc", "1_lt10k", "10k_lt15k"), collapse="|"),
                                   names(inc_nat_vec))]
  nonpov_inc_vec <- inc_nat_vec[grepl(paste(inc_levels_nonpov, collapse="|"),
                                      names(inc_nat_vec))]
  
  # 3. Apply income status
  ag_list[[1]] <- do.call("rbind", lapply(ag_list[[1]], inc_lapply, ht= ht, 
                                          v= nonpov_inc_vec, levels= inc_levels_nonpov))
  ag_list[[2]] <- do.call("rbind", lapply(ag_list[[2]], inc_lapply, ht= ht, 
                                     v= pov_inc_vec, levels= inc_levels_pov))
  # normalize probabilities to account for % < pov , pct >= pov
  pov_rates <- tapply(dat$p, dat$pov_status, sum, na.rm=TRUE)
  pov_rates <- ifelse(is.na(pov_rates), 0, pov_rates)
  ag_names <- names(ag_list)
  if (sum(ag_list[[1]]$p, na.rm=TRUE) < pov_rates[which(names(pov_rates) == names(ag_list)[1])]) {
    ag_list[[1]]$p <- ag_list[[1]]$p * (pov_rates[1] / sum(ag_list[[1]]$p, na.rm=TRUE)) }
  if (sum(ag_list[[2]]$p, na.rm=TRUE) < pov_rates[which(names(pov_rates) == names(ag_list)[2])]) {
    ag_list[[2]]$p <- ag_list[[2]]$p * (pov_rates[2] / sum(ag_list[[2]]$p, na.rm=TRUE)) }
  
  
  ag_list <- do.call("rbind", ag_list)
  ag_list <- factor_return(ag_list, prob_name= "p")
  return(list(ag_list))
}


# helper function 
inc_lapply <- function(l, ht, v, levels) { # need to normalize by poverty rates w/in groups
  if (nrow(l) < 1) return(l)
  levels(l$nativity) <- c("born_other_state","born_out_us", "born_state_residence", "foreigner")
  l_comp <- ht[,2][which(l$nativity[1] == ht[,1])]
  comp <- v[which(grepl(l_comp, names(v)))]
  if (sum(comp) > 0) comp <- (comp / sum(comp)) 
  
  st <- data.frame(pct= comp, levels= factor(levels, levels= levels))
  st <- base::split(st, 1:nrow(st))
  
  dat <- replicate(length(levels), l, simplify = FALSE)
  dat <- do.call("rbind", mapply(add_synth_attr_level, dat= dat, prob_name= "p", attr= st,
                                 attr_name= "ind_income", SIMPLIFY = FALSE))
  return(dat)
}

