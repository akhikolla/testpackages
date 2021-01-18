
# @param agm_dat a \code{list} equivalent to the output of \code{synth_data_mar}
# @param edu_vec A vector containing counts of the total, male, and female population (age 18+)
# as well as counts breaking out the gender counts by age and education status. 
# should be: \code{unlist(<synth_data>$estimates$edu[<row i>,])}
synth_data_edu <- function(agm_dat, edu_vec) {
  dat <- agm_dat[[1]]
  
  # 1. create hash table of age/gender ages to education ages
  age_ht <- data.frame(
              age_gen= c("under15", "15_17", "18_24", "25_29", "30_34","35_39", "40_44", 
                         "45_49", "50_54", "55_59", "60_64", "65_69","70_74", "75_79", "80_84", "85up"),
              edu= c(NA, NA, "18_24", rep("25_34",2), rep("35_44", 2), rep("45_64", 4), rep("65up",5)),
              stringsAsFactors = FALSE)
  
  # 2. create age/gender buckets on which to condition
  agm_list <- split(dat, dat$gender)
  agm_list[[1]] <- split(agm_list[[1]], agm_list[[1]]$age)
  agm_list[[2]] <- split(agm_list[[2]], agm_list[[2]]$age)
  
  edu_m <- edu_vec[which(substr(names(edu_vec), 1,1) == "m")]
  edu_f <- edu_vec[which(substr(names(edu_vec), 1,1) == "f")]
  
  # 3. Apply education levels
  # assumptions: under 15 == less than high school education
  # assumptions: 15-17 == some high school education
  edu_levels <- c("lt_hs", "some_hs", "hs_grad", "some_col", "assoc_dec", "ba_deg", "grad_deg")
  
  agm_list[[1]] <- do.call("rbind", 
       lapply(agm_list[[1]], edu_lapply, ht= age_ht, edu_v= edu_m, levels= edu_levels))
  agm_list[[2]] <- do.call("rbind", 
       lapply(agm_list[[2]], edu_lapply, ht= age_ht, edu_v= edu_f, levels= edu_levels))
  
  dat <- do.call("rbind", agm_list)
  dat <- factor_return(dat, prob_name= "p")
  return(list(dat, levels(dat$age)))
}


# helper function for synth_data_edu. 
edu_lapply <- function(l, ht, edu_v, levels) {
  if (is.na(l$age[1])) # error catch, break/next not allowed in lapply
    return(data.frame(age= "under15", gender= "Male", 
                      marital_status= "never_mar", edu_attain= "lt_hs", p= 0)) 
  else if (l$age[1] == "under15") {
    return(data.frame(age=l$age, gender= l$gender, 
                      marital_status= l$marital_status, edu_attain= "lt_hs", p= l$p))
  } else if (l$age[1] == "15_17") {
    return(data.frame(age=l$age, gender= l$gender, 
                      marital_status= l$marital_status, edu_attain= "some_hs", p= l$p))
  } else {
    l_age_comp <- ht[,2][which(l$age[1] == ht[,1])]
    edu_comp <- edu_v[which(grepl(l_age_comp, names(edu_v)))]
    if (sum(edu_comp) > 0) edu_comp <- (edu_comp / sum(edu_comp)) 
    
    st <- data.frame(pct= edu_comp, levels= factor(levels, levels= levels))
    st <- base::split(st, 1:nrow(st))
    
    dat <- replicate(length(levels), l, simplify = FALSE)
    return(do.call("rbind", mapply(add_synth_attr_level, dat= dat, prob_name= "p", 
                                   attr_name= "edu_attain", attr= st, 
                                   SIMPLIFY = FALSE)))
  }
}

