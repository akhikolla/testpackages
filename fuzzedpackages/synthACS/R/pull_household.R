
#' @title Pull ACS data on households and housing units
#' @description Pull ACS data for a specified geography from base tables
#' B09019, B11011, B19081, B25002, B25003, B25004, B25010, B25024,
#' B25056, B25058, B25071, and B27001.
#' Additional fields, mainly percentages and aggregations, are calculated.
#' @param endyear An integer, indicating the latest year of the data in the survey.
#' @param span An integer in \code{c(1,3,5)} indicating the span of the desired data.
#' @param geography a valid \code{geo.set} object specifying the census geography or 
#' geographies to be fetched.
#' @return A \code{list} containing the endyear, span, a \code{data.frame} of estimates,
#' a \code{data.frame} of standard errors, a character vector of the original column names,
#' and a \code{data.frame} of the geography metadata from \code{\link[acs]{acs.fetch}}.
#' @export
#' @seealso \code{\link[acs]{acs.fetch}}, \code{\link[acs]{geo.make}}
pull_household <- function(endyear, span, geography) {
  # 00 -- error checking
  #----------------------------------------------
  check_geo_inputs(endyear= endyear, span= span, geography= geography)
  
  # 01 -- pull data and move to lists
  #----------------------------------------------
  oldw <- getOption("warn")
  options(warn= -1) # suppress warnings from library(acs) / ACS API
  on.exit(options(warn= oldw)) # turn warnings back on
  hh_type_r <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                         table.number = "B09019", col.names= "pretty")
  hh_type_units <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                         table.number = "B11011", col.names= "pretty")
  hh_inc <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                     table.number = "B19081", col.names= "pretty")
  hh_occ <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                      table.number = "B25002", col.names= "pretty")
  hh_tenure <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                      table.number = "B25003", col.names= "pretty")
  hh_vacancy <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                      table.number = "B25004", col.names= "pretty")
  hh_num_units <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                          table.number = "B25024", col.names= "pretty")
  hh_rent <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                            table.number = "B25056", col.names= "pretty")
  hh_med_rent <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                                table.number = "B25058", col.names= "pretty")
  hh_med_rent_v_inc <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                                      table.number = "B25071", col.names= "pretty")
  health_ins <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                          table.number = "B27001", col.names= "pretty")
  
  # --- these tables not available
  #' B28001 - TYPES OF COMPUTERS IN HOUSEHOLD
  #' B28002 - PRESENCE AND TYPES OF INTERNET SUBSCRIPTIONS IN HOUSEHOLD
   
  est <- list(hh_type_r= data.frame(hh_type_r@estimate),
              hh_type_units= data.frame(hh_type_units@estimate),
              hh_inc= data.frame(hh_inc@estimate),
              hh_occ= data.frame(hh_occ@estimate),
              hh_tenure= data.frame(hh_tenure@estimate),
              hh_vacancy= data.frame(hh_vacancy@estimate),
              hh_num_units= data.frame(hh_num_units@estimate),
              hh_rent= data.frame(hh_rent@estimate),
              health_ins= data.frame(health_ins@estimate))
  
  se <- list(hh_type_r= data.frame(hh_type_r@standard.error),
              hh_type_units= data.frame(hh_type_units@standard.error),
              hh_inc= data.frame(hh_inc@standard.error),
              hh_occ= data.frame(hh_occ@standard.error),
              hh_tenure= data.frame(hh_tenure@standard.error),
              hh_vacancy= data.frame(hh_vacancy@standard.error),
              hh_num_units= data.frame(hh_num_units@standard.error),
              hh_rent= data.frame(hh_rent@standard.error),
              health_ins= data.frame(health_ins@standard.error))
  
  geo <- hh_type_r@geography
  
  
  ## 02 (A) combine columns and (B) calc percentages
  #----------------------------------------------
  ### hh_type_r
  #----------------------------------------------
  est$hh_type_r <- est$hh_type_r[, -c(2,5,6,9:11,16:18,22:23,26:32,36:37)]
  se$hh_type_r  <- se$hh_type_r[, -c(2,5,6,9:11,16:18,22:23,26:32,36:37)]
  names(est$hh_type_r) <- names(se$hh_type_r) <- c("total", paste(
    c(rep("fam_hh", 11), rep("nonfam_hh", 5)),
    c("all", "headofhh", "spouse", "child", "grandchild", "sibling", "parent",
      "inlaw", "boarder", "roommate", "unmarried_partner", "all", "headofh", "boarder",
      "roommate", "unmarried_partner"), sep= "_"), "in_grp_qtrs")
  
  ### hh_type_units
  #----------------------------------------------
  est$hh_type_units$oth_hh_1parent_sing_fam_home <- apply(est$hh_type_units[, c(9,13)], 1 ,sum)
  est$hh_type_units$oth_hh_1parent_mult_fam_home <- apply(est$hh_type_units[, c(10,14)], 1 ,sum)
  est$hh_type_units$oth_hh_1parent_mobile_oth_home <- apply(est$hh_type_units[, c(11,15)], 1 ,sum)
  
  se$hh_type_units$oth_hh_1parent_sing_fam_home <- sqrt(apply(se$hh_type_units[, c(9,13)]^2, 1 ,sum))
  se$hh_type_units$oth_hh_1parent_mult_fam_home <- sqrt(apply(se$hh_type_units[, c(10,14)]^2, 1 ,sum))
  se$hh_type_units$oth_hh_1parent_mobile_oth_home <- sqrt(apply(se$hh_type_units[, c(11,15)]^2, 1 ,sum))
  
  est$hh_type_units <- est$hh_type_units[, c(1,2,4:6,7,20:22,16:19)]
  se$hh_type_units  <- se$hh_type_units[, c(1,2,4:6,7,20:22,16:19)]
  
  names(est$hh_type_units) <- names(se$hh_type_units) <- c("total",
    paste(rep(c("mar_couple_hh", "oth_fam_hh", "nonfam_hh"), each= 4),
          rep(c("all", "sing_fam_home", "mult_fam_home", "mobile_oth_home"), 3), sep= "_"))
  
  ## pcts
  est$hh_type_units$pct_sing_fam_home <- apply(est$hh_type_units[, c(3,7,11)], 1 ,sum)/ est$hh_type_units$total
  est$hh_type_units$pct_mult_fam_home <- apply(est$hh_type_units[, c(4,8,12)], 1 ,sum)/ est$hh_type_units$total
  
  se$hh_type_units$pct_sing_fam_home <- sqrt(apply(se$hh_type_units[, c(3,7,11)]^2, 1 ,sum) - (
    est$hh_type_units$pct_sing_fam_home^2 * se$hh_type_units$total^2)) / est$hh_type_units$total
  se$hh_type_units$pct_mult_fam_home <- sqrt(apply(se$hh_type_units[, c(4,8,12)]^2, 1 ,sum) - (
    est$hh_type_units$pct_mult_fam_home^2 * se$hh_type_units$total^2)) / est$hh_type_units$total
  
  ### hh_inc
  #----------------------------------------------
  names(est$hh_inc) <- names(se$hh_inc) <- c("mean_hh_inc_bottom_quintile",
                                             "mean_hh_inc_2nd_quintile",
                                             "mean_hh_inc_3rd_quintile",
                                             "mean_hh_inc_4th_quintile",
                                             "mean_hh_inc_top_quintile",
                                             "mean_hh_inc_top_5pct")
  gini_calc <- function(row) {
    v <- unlist(row)
    dim(v) <- c(1, length(v))
    v[5] <- 4/3 * (v[5] - .25*v[6]) # adj to mean income of 80th-95pctile
    val <- apply(v, 2, function(i) abs(i-v))
    p <- c(.2,.2,.2,.2,.15,.05) 
    dim(p) <- c(6,1) 
    return(sum(p %*% t(p) * val) / (2 * v[3])) # based on mean(inc) == mean inc of 3rd quintile
  }
  
  est$hh_inc$pseudo_gini_coef <- apply(est$hh_inc, 1, gini_calc)
  se$hh_inc$pseudo_gini_coef <- as.numeric(NA)
  
  ### hh_occ
  #----------------------------------------------
  names(est$hh_occ) <- names(se$hh_occ) <- c("total_housing_units", 
                                             "occupied_housing_units", 
                                             "vacant_housing_units")
  
  est$hh_occ$pct_occupied <- est$hh_occ[,2] / est$hh_occ[,1]
  est$hh_occ$pct_vacant   <- est$hh_occ[,3] / est$hh_occ[,1]
  
  se$hh_occ$pct_occupied <- sqrt(se$hh_occ[,2]^2 - (est$hh_occ$pct_occupied * 
                                 se$hh_occ[,1]^2)) / est$hh_occ[,1]
  se$hh_occ$pct_vacant   <- sqrt(se$hh_occ[,3]^2 - (est$hh_occ$pct_vacant * 
                                 se$hh_occ[,1]^2)) / est$hh_occ[,1]
  
  ### hh_tenure
  #----------------------------------------------
  names(est$hh_tenure) <- names(se$hh_tenure) <- c("total_tenure_units", 
                                             "owner_occupied_units", 
                                             "renter_occupied_units")
  
  est$hh_tenure$pct_owner_occupied <- est$hh_tenure[,2] / est$hh_tenure[,1]
  est$hh_tenure$pct_renter_occupied   <- est$hh_tenure[,3] / est$hh_tenure[,1]
  
  se$hh_tenure$pct_owner_occupied <- sqrt(se$hh_tenure[,2]^2 - (est$hh_tenure$pct_owner_occupied * 
                                                      se$hh_tenure[,1]^2)) / est$hh_tenure[,1]
  se$hh_tenure$pct_renter_occupied   <- sqrt(se$hh_tenure[,3]^2 - (est$hh_tenure$pct_renter_occupied * 
                                                      se$hh_tenure[,1]^2)) / est$hh_tenure[,1]
  
  ### hh_vacancy
  #----------------------------------------------
  names(est$hh_vacancy) <- names(se$hh_vacancy) <- paste("vacant_units",
    c("all", "for_rent", "rented_unoccupied", "for_sale", "sold_unoccupied",
      "seasonal_use", "migrant_workers", "other_vacant"), sep= "_")
  
  ### hh_num_units
  #----------------------------------------------
  est$hh_num_units$one <- apply(est$hh_num_units[,2:3],1,sum)
  est$hh_num_units <- est$hh_num_units[, c(1,12,4:11)]
  se$hh_num_units$one <- sqrt(apply(se$hh_num_units[,2:3]^2,1,sum))
  se$hh_num_units <- se$hh_num_units[, c(1,12,4:11)]
  
  names(est$hh_num_units) <- names(se$hh_num_units) <- paste("building_units", 
    c("all", "1", "2", "3_4", "5_9", "10_19", "20_49", "50up", "mobile_home", "boat_rv_etc"), sep= "_")
  
  
  ### hh_rent
  #----------------------------------------------
  rent_est <- data.frame(
    rent_units= est$hh_rent[,1],
    rent_cash_rent= est$hh_rent[,2],
    rent_nocash_rent= est$hh_rent[,24],
    med_rent= hh_med_rent@estimate[,1],
    med_rent_to_income= hh_med_rent_v_inc@estimate[,1],
    pct_rent_lt500= apply(est$hh_rent[, 3:11], 1, sum) / est$hh_rent[,2],
    pct_rent_500_lt750= apply(est$hh_rent[, 12:16], 1, sum) / est$hh_rent[,2],
    pct_rent_750_lt1000= apply(est$hh_rent[, 17:19], 1, sum) / est$hh_rent[,2],
    pct_rent_1000_lt1250= est$hh_rent[, 20] / est$hh_rent[,2],
    pct_rent_1250_lt1500= est$hh_rent[, 21] / est$hh_rent[,2],
    pct_rent_1500_lt2000= est$hh_rent[, 22] / est$hh_rent[,2],
    pct_rent_2000up= est$hh_rent[, 23] / est$hh_rent[,2]
  )
  
  rent_se <- data.frame(
    rent_units= est$hh_rent[,1],
    rent_cash_rent= est$hh_rent[,2],
    rent_nocash_rent= est$hh_rent[,24],
    med_rent= hh_med_rent@standard.error[,1],
    med_rent_to_income= hh_med_rent_v_inc@standard.error[,1],
    pct_rent_lt500= sqrt(apply(se$hh_rent[, 3:11]^2, 1, sum) - (
      rent_est$pct_rent_lt500^2 * se$hh_rent[,2]^2)) / est$hh_rent[,2],
    pct_rent_500_lt750= sqrt(apply(se$hh_rent[, 12:16]^2, 1, sum) - (
      rent_est$pct_rent_500_lt750^2 * se$hh_rent[,2]^2)) / est$hh_rent[,2],
    pct_rent_750_lt1000= sqrt(apply(se$hh_rent[, 17:19]^2, 1, sum) - (
      rent_est$pct_rent_750_lt1000^2 * se$hh_rent[,2]^2)) / est$hh_rent[,2],
    pct_rent_1000_lt1250= sqrt(se$hh_rent[, 20]^2 - (rent_est$pct_rent_1000_lt1250^2 * 
      se$hh_rent[,2]^2)) / est$hh_rent[,2],
    pct_rent_1250_lt1500= sqrt(se$hh_rent[, 21]^2 - (rent_est$pct_rent_1250_lt1500^2 * 
      se$hh_rent[,2]^2)) / est$hh_rent[,2],
    pct_rent_1500_lt2000= sqrt(se$hh_rent[, 22]^2 - (rent_est$pct_rent_1500_lt2000^2 * 
      se$hh_rent[,2]^2)) / est$hh_rent[,2],
    pct_rent_2000up=      sqrt(se$hh_rent[, 23]^2 - (rent_est$pct_rent_2000up^2 * 
      se$hh_rent[,2]^2)) / est$hh_rent[,2]
  )
  
  est$hh_rent <- rent_est
  se$hh_rent  <- rent_se
  rm(rent_est, rent_se)
  
  ### health_ins
  #----------------------------------------------
  est$health_ins <- est$health_ins[, -c(seq(3,27,3),seq(31,55,3))]
  se$health_ins  <- se$health_ins[, -c(seq(3,27,3),seq(31,55,3))]
  names(est$health_ins) <- names(se$health_ins) <- c("total", paste(
    rep(c("m", "f"), each= 19),
    rep(c("",rep(c("u6", "6_17", "18_24", "25_34", "35_44", "45_54", "55_64", "65_74", "75up"), each=2)), 2),
    rep(c("", rep(c("w_health_ins", "wo_health_ins"), 9)), 2), sep= "_"))
  
  ### pcts
  est$health_ins$pct_w_health_ins <- apply(est$health_ins[,c(seq(3,19,2), seq(22,38,2))], 1, sum) / 
    est$health_ins[,1]
  est$health_ins$pct_m_w_health_ins <- apply(est$health_ins[,seq(3,19,2)], 1, sum) / 
    est$health_ins[,2]
  est$health_ins$pct_f_w_health_ins <- apply(est$health_ins[,seq(22,38,2)], 1, sum) / 
    est$health_ins[,21]
  est$health_ins$pct_w_health_ins_u18 <- apply(est$health_ins[,c(3,5,22,24)], 1, sum) /
                                          apply(est$health_ins[,c(3:6,22:25)], 1, sum)
  est$health_ins$pct_w_health_ins_18_54 <- apply(est$health_ins[,c(seq(7,13,2), seq(26,32,2))], 1, sum) /
                                            apply(est$health_ins[,c(7:14,26:33)], 1, sum)
  est$health_ins$pct_w_health_ins_55_64 <- apply(est$health_ins[,c(15,34)], 1, sum) /
                                            apply(est$health_ins[,c(15:16,34:35)], 1, sum)
  
  se$health_ins$pct_w_health_ins <- sqrt(
    apply(se$health_ins[,c(seq(3,19,2), seq(22,38,2))]^2, 1, sum) - (
      est$health_ins$pct_w_health_ins^2 * se$health_ins[,1]^2)
    ) / est$health_ins[,1]
  
  se$health_ins$pct_m_w_health_ins <- sqrt(
    apply(se$health_ins[,seq(3,19,2)]^2, 1, sum) - (
      est$health_ins$pct_m_w_health_ins^2 * se$health_ins[,2]^2)
    ) / est$health_ins[,2]
  
  se$health_ins$pct_f_w_health_ins <- sqrt(
    apply(se$health_ins[,seq(22,38,2)]^2, 1, sum) - (
      est$health_ins$pct_f_w_health_ins^2 * se$health_ins[,2]^2)
    ) /  est$health_ins[,21]
  
  
  se$health_ins$pct_w_health_ins_u18 <- sqrt(
    apply(se$health_ins[,c(3,5,22,24)]^2, 1, sum) - (
      est$health_ins$pct_w_health_ins_u18^2 * apply(se$health_ins[,c(3:6,22:25)]^2, 1, sum))
    ) / apply(est$health_ins[,c(3:6,22:25)], 1, sum)
  
  se$health_ins$pct_w_health_ins_18_54 <- sqrt(
    apply(se$health_ins[,c(seq(7,13,2), seq(26,32,2))]^2, 1, sum) - (
      est$health_ins$pct_w_health_ins_18_54^2 * apply(se$health_ins[,c(7:14,26:33)]^2, 1, sum))
    ) / apply(est$health_ins[,c(7:14,26:33)], 1, sum)
  
  se$health_ins$pct_w_health_ins_55_64 <- sqrt(
    apply(se$health_ins[,c(15,34)]^2, 1, sum) - (
      est$health_ins$pct_w_health_ins_55_64^2 * apply(se$health_ins[,c(15:16,34:35)]^2, 1, sum))
    ) / apply(est$health_ins[,c(15:16,34:35)], 1, sum)
  
  
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
  names(ret$estimates) <- names(ret$standard_error) <- c("hh_type_by_relationship", "hh_type_by_units",
    "mean_hh_inc_quintiles", "hh_occ_status", "hh_tenure", "hh_vacancy_status", # "avg_hh_size",
    "units_in_structure", "contract_rent", "hh_ins_by_sex_age")
  
  return(ret)
}
