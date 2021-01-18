# @param age_gender_vec A vector containing counts of the total, male, and female population as well
# as counts breaking out the gender counts by age. 
# should be. \code{unlist(<synth_data>$estimates$age_by_sex[<row i>,])}
# @return A list containing two elements: (1) the synthetic data of age, gender, and associated
# probabilities, and (2) an integer length of the number of age labels.
synth_data_ag <- function(age_gender_vec) {
  # 00. error checking
  if (!sum(grepl("m_", names(age_gender_vec))) > 0) 
    stop("age_gender_vec must include entries with regex 'm_'.")
  if (!sum(grepl("f_", names(age_gender_vec))) > 0) 
    stop("age_gender_vec must include entries with regex 'f_'.")
  if (!sum(grepl("cnt_all", names(age_gender_vec))) > 0) 
    stop("age_gender_vec must include a total count entry: regex 'cnt_all'.")
  if(any(age_gender_vec %% 1 != 0)) stop("all elements of age_gender_vec must be integers.")
  
  # 01 & 02. gender and age
  #------------------------------------
  gender_labels <- c("Male", "Female")
  age_labels <- c("under15", "15_17", "18_24", "25_29", "30_34","35_39", "40_44", 
                  "45_49", "50_54", "55_59", "60_64", "65_69","70_74", "75_79", "80_84", "85up")
  n_ages <- length(age_labels)
  
  # pattern matching for counts, the convert to percentages
  ag_vec <- age_gender_vec[which(substr(names(age_gender_vec), 1,2) %in% c("m_", "f_"))]
  cnt_total <- age_gender_vec[grep("cnt_all", names(age_gender_vec))]
  p_gender_age <- ag_vec / cnt_total
  
  dat <- data.frame(expand.grid(age=age_labels, gender= gender_labels), p=p_gender_age)
  dat <- factor_return(dat, prob_name= "p")
  return(list(dat))
}

# @param ag_dat a list equivalent to the output of \code{synth_data_ag} (above)
# @param mar_status_vec A vector containing counts of total population in a given geography 
# (age 15+) with breakouts by age, gender, and marital status
# should be \code{unlist(<synth_data>$estimates$marital_status[<row i>,])}
synth_data_mar <- function(ag_dat, mar_status_vec) {
  # 00. error checking
  if (!sum(grepl("m_", names(mar_status_vec))) > 0) 
    stop("mar_status_vec must include entries with regex 'm_'.")
  if (!sum(grepl("f_", names(mar_status_vec))) > 0) 
    stop("mar_status_vec must include entries with regex 'f_'.")
  if(any(mar_status_vec %% 1 != 0)) stop("all elements of mar_status_vec must be integers.")
  
  
  dat <- ag_dat[[1]]
  # 1. create hash table of age/gender ages to employment status ages
  ht <- data.frame(age_old= c("under15", "15_17", "18_24", "25_29", "30_34","35_39", "40_44", 
                              "45_49", "50_54", "55_59", "60_64", "65_69","70_74", "75_79", "80_84", "85up"),
                   age_new= c(NA, "15_17", "18_24", "25_29", "30_34", "35_39", "40_44", "45_49",
                              "50_54", "55_59", "60_64", rep("65_74", 2), rep("75_84", 2), "85up"))
  
  # 2. create age/gender buckets on which to condition
  ag_list <- split(dat, dat$gender)
  ag_list[[1]] <- split(ag_list[[1]], ag_list[[1]]$age)
  ag_list[[2]] <- split(ag_list[[2]], ag_list[[2]]$age)
  
  m_mar_vec <- mar_status_vec[which(substr(names(mar_status_vec), 1,1) == "m")]
  f_mar_vec <- mar_status_vec[which(substr(names(mar_status_vec), 1,1) == "f")]
  
  # 3. apply marital status
  # next, loop through males, then females 
  mar_levels <- c("never_mar", "married", "mar_apart", "widowed", "divorced")
  
  ag_list[[1]] <- do.call("rbind", lapply(ag_list[[1]], mar_lapply, ht= ht, 
                                 v= m_mar_vec, levels= mar_levels))
  ag_list[[2]] <- do.call("rbind", lapply(ag_list[[2]], mar_lapply, ht= ht, 
                                 v= f_mar_vec, levels= mar_levels))
  
  dat <- do.call("rbind", ag_list)
  dat <- factor_return(dat, prob_name= "p")
  return(list(dat, levels(dat$age)))
}


# helper function for synth_data_emp. 
mar_lapply <- function(l, ht, v, levels) {
  if (is.na(l$age[1])) # error catch, break/next not allowed in lapply
    return(data.frame(age= "under15", gender= "Male", 
                      marital_status= "never_mar", p= 0)) 
  # first assume M/F under 15 are never married
  else if (l$age[1] == "under15") {
    return(data.frame(age= l$age,
                      gender= l$gender,
                      marital_status= "never_mar", 
                      p= l$p))
  } else {
    l_comp <- ht[,2][which(l$age[1] == ht[,1])]
    comp <- v[which(grepl(l_comp, names(v)))]
    if (sum(comp) > 0) comp <- (comp / sum(comp)) 
    
    st <- data.frame(pct= comp, levels= factor(levels, levels= levels))
    st <- base::split(st, 1:nrow(st))
    
    dat <- replicate(length(levels), l, simplify = FALSE)
    dat <- do.call("rbind", mapply(add_synth_attr_level, dat= dat, prob_name= "p", 
                                   attr_name= "marital_status", attr= st, 
                                   SIMPLIFY = FALSE))
    return(dat)
  }
}
