
# @param agmeen_dat a \code{list} equivalent to the output of \code{synth_data_nativ}
# @param pov_ge_vec A vector containing counts of the total population 
# (civilian individuals 16+ years for whom poverty status is determined)
# as well as counts breaking out the counts by gender and employment status. 
# should be: \code{unlist(<synth_data>$estimates$pov_status1[<row i>,])}
# NOTE: employment status is determined by Age/gender
synth_data_pov <- function(agmeen_dat, pov_ge_vec, total_pop) {
  
  # 0. check population on which poverty reported vs total population
  # if (pop / pov_pop) > 3) | (pop / pov_pov > 2 & pov_pov < 100) --> make 100% >= pov line
  pov_pop <- pov_ge_vec[1]
  if (total_pop / pov_pop > 3 | (total_pop / pov_pop > 2 & pov_pop < 100)) {
    message("Data on poverty status is underreported. Marking all individauls as >= poverty line.")
    dat <- add_synth_attr_level(dat= agmeen_dat[[1]], prob_name= "p", attr_name= "pov_status",
                                attr= list(pct= 1.0, level="at_above_pov_level"))
    dat$pov_status <- factor(dat$pov_status, levels= c("below_pov_level", "at_above_pov_level"))
    return(list(dat, levels(dat$edu_attain)))
  }
  
  dat <- agmeen_dat[[1]]
  # 1. create age hash table, marginalize pov by employment status / gender
  age_ht <- data.frame(age_dat= agmeen_dat[[2]],
                       age_new= c("15_17", "15_17", "18_24", rep("25_34", 2), rep("35_44",2), 
                                  rep("45_54", 2), rep("55_64", 2), rep("65_75", 2), rep("75up", 3)),
                       stringsAsFactors = FALSE)
  
  # marginalize pov by gender/employment status by emp status (sep by gender)
  pv_m <- pov_ge_vec[which(substr(names(pov_ge_vec), 1,1) == "m")]
  pv_f <- pov_ge_vec[which(substr(names(pov_ge_vec), 1,1) == "f")]
  
  emp_lvls <- levels(agmeen_dat[[1]]$emp_status)
  emp_marg_m <- marginalize_emp_status(pv_m, emp_lvls)
  emp_marg_f <- marginalize_emp_status(pv_f, emp_lvls)
  
  # 2. create buckets on which to condition -- gender & employment status
  ag_list <- split(dat, dat$gender)
  ag_list[[1]] <- split(ag_list[[1]], ag_list[[1]]$emp_status)
  ag_list[[2]] <- split(ag_list[[2]], ag_list[[2]]$emp_status)
  
  # 3. Apply poverty status
  pov_levels <- c("below_pov_level", "at_above_pov_level")
  
  ag_list[[1]] <- do.call("rbind", lapply(ag_list[[1]], pov_lapply, emp_marg= emp_marg_m,
                                          levels= pov_levels))
  ag_list[[2]] <- do.call("rbind", lapply(ag_list[[2]], pov_lapply, emp_marg= emp_marg_f,
                                          levels= pov_levels))
  
  dat <- do.call("rbind", ag_list)
  dat <- factor_return(dat, prob_name= "p")
  return(list(dat, levels(dat$edu_attain)))
}


# helper function for synth_data_pov. 
# @param l a list of data.frames -- internal to synth_data_pov. In this case either ag_list[[1]]
# or ag_list[[2]]
# @param levels the levels of the new variable to be added (poverty in this case)
# @param emp_marg the poverty rates by employment type
pov_lapply <- function(l, levels, emp_marg) {
  if (nrow(l) < 1) return(l)
  emp_levels <- names(table(l$emp_status))[table(l$emp_status) > 0]
  comp <- emp_marg[which(names(emp_marg) %in% emp_levels)]
  comp <- c(comp, 1-comp)
  
  st <- data.frame(pct= comp, levels= factor(levels, levels= levels))
  st <- base::split(st, 1:nrow(st))
  
  dat <- replicate(length(emp_levels), l, simplify = FALSE)
  dat <- do.call("rbind", mapply(add_synth_attr_level, dat= dat, prob_name= "p", attr= st,
                                 attr_name= "pov_status", SIMPLIFY = FALSE))
  return(dat)
}

# after first splitting the vector of poverty status by gender/employment status 
# this function takes a vector of poverty statuses by employment status (joint density)
# for a single gender and produces the marginal density by of poverty for each 
# employment status
# the result is a vector of poverty status by 
marginalize_emp_status <- function(pov_vec, emp_levels) {
  n <- length(emp_levels)
  out <- vector("numeric", length= n); names(out) <- emp_levels
  for (i in 1:n) { # marginalize by looping through employment statuses
    l <- emp_levels[i]; l2 <- emp_levels[-i]
    if (l == "employed") { # needed to separate employed from unemployed
      ind <- which(grepl(l, names(pov_vec)) & !grepl(paste(l2, collapse="|"), names(pov_vec)))
    } else {
      ind <- which(grepl(l, names(pov_vec)))
    }
    
    out[i] <- pov_vec[ind][1] / sum(pov_vec[ind]) # returns % lt pov
  }
  return(out)
}
