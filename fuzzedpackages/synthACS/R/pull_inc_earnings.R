
#' @title Pull ACS income and earnings data
#' @description Pull ACS data for a specified geography from base tables
#' B19083, B19301, B19326, B21001, B22001, B23020, B24011. Not yet 
#' implemented: B28004
#' Additional fields, mainly percentages and aggregations, are calculated.
#' @param endyear An integer, indicating the latest year of the data in the survey.
#' @param span An integer in \code{c(1,3,5)} indicating the span of the desired data.
#' @param geography a valid \code{geo.set} object specifying the census geography or 
#' geographies to be fetched.
#' @return A \code{list} containing the endyear, span, a \code{data.frame} of estimates,
#' a \code{data.frame} of standard errors, a character vector of the original column names,
#' and a \code{data.frame} of the geography metadata from \code{\link[acs]{acs.fetch}}.
#' @seealso \code{\link[acs]{acs.fetch}}, \code{\link[acs]{geo.make}}
#' @export
pull_inc_earnings <- function(endyear, span, geography) {
  # 00 -- error checking
  #----------------------------------------------
  check_geo_inputs(endyear= endyear, span= span, geography= geography)
  
  # 01 -- pull data and move to lists
  #----------------------------------------------
  oldw <- getOption("warn")
  options(warn= -1) # suppress warnings from library(acs) / ACS API
  on.exit(options(warn= oldw)) # turn warnings back on
  gini <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                    table.number = "B19083", col.names= "pretty")
  inc_pc <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                    table.number = "B19301", col.names= "pretty")
  med_inc <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                      table.number = "B19326", col.names= "pretty")
  vet_status <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                       table.number = "B21001", col.names= "pretty")
  snap <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                          table.number = "B22001", col.names= "pretty")
  mean_hrs_work <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                          table.number = "B23020", col.names= "pretty")
  med_earn <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                          table.number = "B24011", col.names= "pretty")
#   inc_by_int <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
#                           table.number = "B28004", col.names= "pretty")
  
  est <- list(gini= data.frame(gini@estimate),
              inc_pc= data.frame(inc_pc@estimate),
              med_inc= data.frame(med_inc@estimate),
              vet_status= data.frame(vet_status@estimate),
              snap= data.frame(snap@estimate),
              mean_hrs_work= data.frame(mean_hrs_work@estimate),
              med_earn= data.frame(med_earn@estimate[, c(1,4,5,7:9,11:14, 16:17,19,21:25,27,28,
                                                         30:32,34:36), drop= FALSE]))
  
  se <- list(gini= data.frame(gini@standard.error),
              inc_pc= data.frame(inc_pc@standard.error),
              med_inc= data.frame(med_inc@standard.error),
              vet_status= data.frame(vet_status@standard.error),
              snap= data.frame(snap@standard.error),
              mean_hrs_work= data.frame(mean_hrs_work@standard.error),
              med_earn= data.frame(med_earn@standard.error[, c(1,4,5,7:9,11:14, 16:17,19,21:25,27,28,
                                                         30:32,34:36), drop= FALSE]))
  
  geo <- gini@geography
  
  rm(gini, inc_pc, med_inc, vet_status, snap, mean_hrs_work, med_earn)
  
  ## 02 -- (A) combine columns and (B) calc percentages
  ### gini
  #----------------------------------------------
  names(est$gini) <- names(se$gini) <- "gini_coef"
  
  ### inc_pc
  names(est$inc_pc) <- names(se$inc_pc) <- "inc_per_capita"
  
  ### med_inc
  names(est$med_inc) <- names(se$med_inc) <- c("med_inc", "med_inc_m", "med_inc_m_FTwork",
    "med_inc_m_oth", "med_inc_f", "med_inc_f_FTwork", "med_inc_f_oth")
  
  ### vet_status
  s <- c(7:21, 25:39); s <- s[!(s %in% c(seq(7,19,3), seq(25,37,3)))]
  est$vet_status <- est$vet_status[, c(1:3,s)]
  se$vet_status <- se$vet_status[, c(1:3,s)]
  
  names(est$vet_status) <- names(se$vet_status) <- c("all", "vet", "non_vet", paste(
    rep(c("m", "f"), each=10),
    rep(paste(rep(c("18_34", "35_54", "55_64", "65_74", "75up"), each=2),
    rep(c("vet", "non_vet"), 5),sep= "_"), 2), sep= "_"))
  
  ### snap
  est$snap <- est$snap[, c(1,2,5)]
  se$snap <- se$snap[, c(1,2,5)]
  
  names(est$snap) <- names(se$snap) <- c("hh_total", "hh_w_snap", "hh_no_snap")
  
  ### mean_hrs
  names(est$mean_hrs_work) <- names(se$mean_hrs_work) <- c("mean_hrs_worked_all", "mean_hrs_worked_m",
                                                           "mean_hrs_worked_f")
  
  ### med_earn
  names(est$med_earn) <- names(se$med_earn) <- c("all_occ", "mgmt_occ", "bus_finance_occ", "comp_math_occ",
      "architect_eng_occ", "life_phys_soc_science_occ", "comm_social_serv_occ", "legal_occ", "edu_occ",
      "art_media_occ", "health_practictioners_occ", "health_tech_occ", "heath_support_occ", "fire_occ",
      "police_occ", "food_prep_occ", "cleaning_occ", "personal_care_occ", "sales_occ", "office_admin_occ", 
      "farm_fish_forest_occ", "construct_occ","maint_rep_occ", "production_occ", 
      "transport_occ", "material_move_occ")
  
  # 03 -- sort, combine and return
  #----------------------------------------------
  geo_sorted <- geo_alphabetize(geo= geo, est= est, se= se)
  geo <- geo_sorted[["geo"]]
  est <- geo_sorted[["est"]]
  se <- geo_sorted[["se"]]
  
  ret <- list(endyear= endyear, span= span,
              estimates= est,
              standard_error= se,
              geography= geo,
              geo_title= unlist(geography@geo.list))
  class(ret) <- "macroACS"
  names(ret$estimates) <- names(ret$standard_error) <- c("gini_index", "income_per_capita", 
    "med_inc_by_sex_work_exp", "vet_status_by_sex_age", "snap_receipt", "mean_hrs_worked",
    "occ_by_median_earn")
  
  return(ret)
}