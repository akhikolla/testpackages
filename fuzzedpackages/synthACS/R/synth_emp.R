
# @param agme_dat a \code{list} equivalent to the output of \code{synth_data_edu}
# @param emp_status_vec A vector containing counts of the total, male, and female population 
# (age 16+) as well as counts breaking out the gender counts by age and education status. 
# should be: \code{unlist(<synth_data>$estimates$emp_status[<row i>,])}
synth_data_emp <- function(agme_dat, emp_status_vec) {
  dat <- agme_dat[[1]]
  
  # 1. create hash table of age/gender ages to employment status ages
  age_ht <- data.frame(age_gen= agme_dat[[2]],
                       emp= c(NA, "16_19", "20_24", "25_29", "30_34", rep("35_44",2), rep("45_54",2),
                              "55_59", "60_64", "65_69", "70_74", rep("75up", 3)), 
                       stringsAsFactors = FALSE)
  
  # 2. create age/gender buckets on which to condition
  ag_list <- split(dat, dat$gender)
  ag_list[[1]] <- split(ag_list[[1]], ag_list[[1]]$age)
  ag_list[[2]] <- split(ag_list[[2]], ag_list[[2]]$age)

  emp_m <- emp_status_vec[which(substr(names(emp_status_vec), 1,1) == "m")]
  emp_f <- emp_status_vec[which(substr(names(emp_status_vec), 1,1) == "f")]
  
  # 3. Apply employment status
  # assumptions: under 15 == not in labor force
  # assumptions: 15-17 best represented by 16-19
  # assumptions: 18-24 best represented by employment status of 20-24 (vs 16-19)
  emp_levels <- c("employed", "unemployed", "not_in_labor_force")
  
  ag_list[[1]] <- do.call("rbind", 
                           lapply(ag_list[[1]], emp_lapply, ht= age_ht, 
                                  emp_v= emp_m, levels= emp_levels))
  ag_list[[2]] <- do.call("rbind", 
                           lapply(ag_list[[2]], emp_lapply, ht= age_ht, 
                                  emp_v= emp_f, levels= emp_levels))
  
  dat <- do.call("rbind", ag_list)
  dat <- factor_return(dat, prob_name= "p")
  return(list(dat, levels(dat$age)))
}
  
  
# helper function for synth_data_emp. 
emp_lapply <- function(l, ht, emp_v, levels) {
  if (is.na(l$age[1])) # error catch, break/next not allowed in lapply
    return(data.frame(age= "under15", gender= "Male", 
                      marital_status= "never_mar", edu_attain= "lt_hs",
                      emp_status= "employed", p= 0)) 
  else if (l$age[1] == "under15") {
    return(data.frame(age=l$age, gender= l$gender, 
                      marital_status= l$marital_status, edu_attain= l$edu_attain,
                      emp_status= levels[3],
                      p= l$p))
  } else {
    l_age_comp <- ht[,2][which(l$age[1] == ht[,1])]
    emp_comp <- emp_v[which(grepl(l_age_comp, names(emp_v)))]
    if (sum(emp_comp) > 0) emp_comp <- (emp_comp / sum(emp_comp)) 
    
    st <- data.frame(pct= emp_comp, levels= factor(levels, levels= levels))
    st <- base::split(st, 1:nrow(st))
    
    dat <- replicate(length(levels), l, simplify = FALSE)
    dat <- do.call("rbind", mapply(add_synth_attr_level, dat= dat, prob_name= "p", attr= st,
                                   attr_name= "emp_status", SIMPLIFY = FALSE))
    return(dat)
  }
}
