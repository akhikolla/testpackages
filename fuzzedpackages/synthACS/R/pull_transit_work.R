
#' @title Pull ACS transit and work data
#' @description Pull ACS data for a specified geography from base tables
#' B08012, B08101, B08121, B08103, B08124, B08016, B08017.
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
pull_transit_work <- function(endyear, span, geography) {
  # 00 -- error checking
  #----------------------------------------------
  check_geo_inputs(endyear= endyear, span= span, geography= geography)
  
  # 01 -- pull data and move to lists
  #----------------------------------------------
  oldw <- getOption("warn")
  options(warn= -1) # suppress warnings from library(acs) / ACS API
  on.exit(options(warn= oldw)) # turn warnings back on
  travel_to_work <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                              table.number = "B08012", col.names= "pretty")
  
  mode_travel_to_work <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                                   table.number = "B08101", col.names= "pretty")
  
  mode_transit_earn <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                                 table.number = "B08121", col.names= "pretty")
  
  med_age_transit <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                               table.number = "B08103", col.names= "pretty")
  
  transit_by_occ <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                              table.number = "B08124", col.names= "pretty")
  
  place_of_work_msa <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                                 table.number = "B08016", col.names= "pretty")
  
  place_of_work_microsa <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                                     table.number = "B08017", col.names= "pretty")
  
  est <- list(time_to_work_by_sex= data.frame(travel_to_work@estimate),
              mode_to_work_by_age= data.frame(mode_travel_to_work@estimate),
              transit_med_earn= data.frame(mode_transit_earn@estimate),
              mode_by_med_age= data.frame(med_age_transit@estimate),
              by_occ= data.frame(transit_by_occ@estimate),
              msa_place= data.frame(place_of_work_msa@estimate),
              microsa_place= data.frame(place_of_work_microsa@estimate))
  
  se <- list(time_to_work_by_sex= data.frame(travel_to_work@standard.error),
              mode_to_work_by_age= data.frame(mode_travel_to_work@standard.error),
              transit_med_earn= data.frame(mode_transit_earn@standard.error),
              mode_by_med_age= data.frame(med_age_transit@standard.error),
              by_occ= data.frame(transit_by_occ@standard.error),
              msa_place= data.frame(place_of_work_msa@standard.error),
              microsa_place= data.frame(place_of_work_microsa@standard.error))
  
  geo <- travel_to_work@geography
  
  rm(travel_to_work, mode_travel_to_work, mode_transit_earn, med_age_transit, 
     transit_by_occ, place_of_work_microsa, place_of_work_msa)
  
  
  ## 02 -- (A) combine columns and (B) calc percentages
  ### time_to_work_by_sex
  #----------------------------------------------
  est$time_to_work_by_sex$tr_lt10      <- apply(est$time_to_work_by_sex[, 2:3], 1, sum)
  est$time_to_work_by_sex$tr_10_lt30   <- apply(est$time_to_work_by_sex[, 4:7], 1, sum)
  est$time_to_work_by_sex$tr_30_lt45   <- apply(est$time_to_work_by_sex[, 8:10], 1, sum)
  est$time_to_work_by_sex$tr_lt10_m    <- apply(est$time_to_work_by_sex[, 15:16], 1, sum)
  est$time_to_work_by_sex$tr_10_lt30_m <- apply(est$time_to_work_by_sex[, 17:20], 1, sum)
  est$time_to_work_by_sex$tr_30_lt45_m <- apply(est$time_to_work_by_sex[, 21:23], 1, sum)
  est$time_to_work_by_sex$tr_lt10_f    <- apply(est$time_to_work_by_sex[, 28:29], 1, sum)
  est$time_to_work_by_sex$tr_10_lt30_f <- apply(est$time_to_work_by_sex[, 30:33], 1, sum)
  est$time_to_work_by_sex$tr_30_lt45_f <- apply(est$time_to_work_by_sex[, 34:36], 1, sum)
  est$time_to_work_by_sex <- est$time_to_work_by_sex[, c(1, 40:42, 11:14, 43:45, 24:27, 46:48,37:39)]
  
  se$time_to_work_by_sex$tr_lt10      <- sqrt(apply(se$time_to_work_by_sex[, 2:3]^2, 1, sum))
  se$time_to_work_by_sex$tr_10_lt30   <- sqrt(apply(se$time_to_work_by_sex[, 4:7]^2, 1, sum))
  se$time_to_work_by_sex$tr_30_lt45   <- sqrt(apply(se$time_to_work_by_sex[, 8:10]^2, 1, sum))
  se$time_to_work_by_sex$tr_lt10_m    <- sqrt(apply(se$time_to_work_by_sex[, 15:16]^2, 1, sum))
  se$time_to_work_by_sex$tr_10_lt30_m <- sqrt(apply(se$time_to_work_by_sex[, 17:20]^2, 1, sum))
  se$time_to_work_by_sex$tr_30_lt45_m <- sqrt(apply(se$time_to_work_by_sex[, 21:23]^2, 1, sum))
  se$time_to_work_by_sex$tr_lt10_f    <- sqrt(apply(se$time_to_work_by_sex[, 28:29]^2, 1, sum))
  se$time_to_work_by_sex$tr_10_lt30_f <- sqrt(apply(se$time_to_work_by_sex[, 30:33]^2, 1, sum))
  se$time_to_work_by_sex$tr_30_lt45_f <- sqrt(apply(se$time_to_work_by_sex[, 34:36]^2, 1, sum))
  se$time_to_work_by_sex <- se$time_to_work_by_sex[, c(1, 40:42, 11:14, 43:45, 24:27, 46:48,37:39)]
  
  names(est$time_to_work_by_sex) <- names(se$time_to_work_by_sex) <-
    paste(rep(c("all", "m", "f"), each= 7),
          rep(c("cnt", "lt10", "10_lt30", "30_lt45", "45_lt60", "60_lt90", "90up"), 3), sep= "_")
  
  ### pcts
  est$time_to_work_by_sex$pct_lt30 <- apply(est$time_to_work_by_sex[, 2:3], 1, sum) / est$time_to_work_by_sex$all_cnt
  est$time_to_work_by_sex$pct_lt60 <- apply(est$time_to_work_by_sex[, 2:5], 1, sum) / est$time_to_work_by_sex$all_cnt
  
  se$time_to_work_by_sex$pct_lt30 <- sqrt(apply(se$time_to_work_by_sex[, 2:3]^2, 1, sum) - 
      ((apply(est$time_to_work_by_sex[, 2:3], 1, sum) / est$time_to_work_by_sex$all_cnt)^2 * 
        se$time_to_work_by_sex$all_cnt^2)) / est$time_to_work_by_sex$all_cnt
  se$time_to_work_by_sex$pct_lt60 <- sqrt(apply(se$time_to_work_by_sex[, 2:5]^2, 1, sum) - 
      ((apply(est$time_to_work_by_sex[, 2:5], 1, sum) / est$time_to_work_by_sex$all_cnt)^2 * 
        se$time_to_work_by_sex$all_cnt^2)) / est$time_to_work_by_sex$all_cnt
  
  
  ## 03 -- (A) combine columns and (B) calc percentages
  ### mode_to_work_by_age
  #----------------------------------------------
  names(est$mode_to_work_by_age) <- names(se$mode_to_work_by_age) <- paste(
    rep(c("cnt", "drove_alone", "carpool", "transit", "walk", "other", "work_home"), each= 8),
    rep(c("all", "16_19", "20_24", "25_44", "45_54", "55_59", "60_64", "65up"), 7), sep= "_")
  
  ### pcts 
  est$mode_to_work_by_age$pct_drove_alone <- est$mode_to_work_by_age[,9] / est$mode_to_work_by_age[,1]
  est$mode_to_work_by_age$pct_carpool     <- est$mode_to_work_by_age[,17] / est$mode_to_work_by_age[,1]
  est$mode_to_work_by_age$pct_transit     <- est$mode_to_work_by_age[,25] / est$mode_to_work_by_age[,1]
  est$mode_to_work_by_age$pct_walk        <- est$mode_to_work_by_age[,33] / est$mode_to_work_by_age[,1]
  est$mode_to_work_by_age$pct_other       <- est$mode_to_work_by_age[,41] / est$mode_to_work_by_age[,1]
  est$mode_to_work_by_age$pct_work_home   <- est$mode_to_work_by_age[,49] / est$mode_to_work_by_age[,1]
  
  se$mode_to_work_by_age$pct_drove_alone <- sqrt(se$mode_to_work_by_age[,9]^2 - 
     (est$mode_to_work_by_age[,9] / est$mode_to_work_by_age[,1])^2 * se$mode_to_work_by_age[,1]^2) / 
        est$mode_to_work_by_age[,1]
  se$mode_to_work_by_age$pct_carpool     <- sqrt(se$mode_to_work_by_age[,17]^2 - 
    (est$mode_to_work_by_age[,17] / est$mode_to_work_by_age[,1])^2 * se$mode_to_work_by_age[,1]^2) / 
      est$mode_to_work_by_age[,1]
  se$mode_to_work_by_age$pct_transit     <- sqrt(se$mode_to_work_by_age[,25]^2 - 
    (est$mode_to_work_by_age[,25] / est$mode_to_work_by_age[,1])^2 * se$mode_to_work_by_age[,1]^2) / 
      est$mode_to_work_by_age[,1]
  se$mode_to_work_by_age$pct_walk        <- sqrt(se$mode_to_work_by_age[,33]^2 - 
    (est$mode_to_work_by_age[,33] / est$mode_to_work_by_age[,1])^2 * se$mode_to_work_by_age[,1]^2) / 
      est$mode_to_work_by_age[,1]
  se$mode_to_work_by_age$pct_other       <- sqrt(se$mode_to_work_by_age[,41]^2 - 
    (est$mode_to_work_by_age[,41] / est$mode_to_work_by_age[,1])^2 * se$mode_to_work_by_age[,1]^2) / 
      est$mode_to_work_by_age[,1]
  se$mode_to_work_by_age$pct_work_home   <- sqrt(se$mode_to_work_by_age[,49]^2 - 
    (est$mode_to_work_by_age[,49] / est$mode_to_work_by_age[,1])^2 * se$mode_to_work_by_age[,1]^2) / 
      est$mode_to_work_by_age[,1]
  
  est$mode_to_work_by_age <- est$mode_to_work_by_age[, -c(2:8)]
  se$mode_to_work_by_age  <- se$mode_to_work_by_age[, -c(2:8)]
  
  ## 04 -- (A) combine columns and (B) calc percentages
  ### transit_med_earn
  #----------------------------------------------
  names(est$transit_med_earn) <- names(se$transit_med_earn) <- c(
    "all_cnt", "drove_alone", "carpool", "transit", "walk", "other", "work_home")
  
  ### pcts
  est$transit_med_earn$drove_alone_to_median <- est$transit_med_earn$drove_alone / est$transit_med_earn$all_cnt
  est$transit_med_earn$carpool_to_median <- est$transit_med_earn$carpool / est$transit_med_earn$all_cnt
  est$transit_med_earn$transit_to_median <- est$transit_med_earn$transit / est$transit_med_earn$all_cnt
  est$transit_med_earn$walk_to_median <- est$transit_med_earn$walk / est$transit_med_earn$all_cnt
  est$transit_med_earn$other_to_median <- est$transit_med_earn$other / est$transit_med_earn$all_cnt
  est$transit_med_earn$work_home_to_median <- est$transit_med_earn$work_home / est$transit_med_earn$all_cnt
  
  se$transit_med_earn$drove_alone_to_median <- sqrt(se$transit_med_earn$drove_alone^2 - (
    est$transit_med_earn$drove_alone_to_median^2 * se$transit_med_earn$all_cnt^2)) / est$transit_med_earn$all_cnt
  se$transit_med_earn$carpool_to_median <- sqrt(se$transit_med_earn$carpool^2 - (
    est$transit_med_earn$carpool_to_median^2 * se$transit_med_earn$all_cnt^2)) / est$transit_med_earn$all_cnt
  se$transit_med_earn$transit_to_median <- sqrt(se$transit_med_earn$transit^2 - (
    est$transit_med_earn$transit_to_median^2 * se$transit_med_earn$all_cnt^2)) / est$transit_med_earn$all_cnt
  se$transit_med_earn$walk_to_median <- sqrt(se$transit_med_earn$walk^2 - (
    est$transit_med_earn$walk_to_median^2 * se$transit_med_earn$all_cnt^2)) / est$transit_med_earn$all_cnt
  se$transit_med_earn$other_to_median <- sqrt(se$transit_med_earn$other^2 - (
    est$transit_med_earn$other_to_median^2 * se$transit_med_earn$all_cnt^2)) / est$transit_med_earn$all_cnt
  se$transit_med_earn$work_home_to_median <- sqrt(se$transit_med_earn$work_home^2 - (
    est$transit_med_earn$work_home_to_median^2 * se$transit_med_earn$all_cnt^2)) / est$transit_med_earn$all_cnt
  
  
  ## 05 -- (A) combine columns and (B) calc percentages
  ### mode_by_med_age
  #----------------------------------------------
  names(est$mode_by_med_age) <- names(se$mode_by_med_age) <- c(
    "median_age", "med_age_drove_alone", "med_age_carpool", "med_age_transit", 
    "med_age_walk", "med_age_other", "med_age_work_home")
  
  ## 06 -- (A) combine columns and (B) calc percentages
  ### by_occ
  #----------------------------------------------
  names(est$by_occ) <- names(se$by_occ) <- paste(
    rep(c("all", "drove_alone", "carpool", "transit", "walk", "other", "work_home"), each= 7),
    rep(c("", "mgmt_bus_science_arts", "service", "sales_office", "resources_constr_maint",
          "production_transport", "mil"), 7), sep= "_")
  
  ### pcts
  est$by_occ$pct_mgmt_bus_science_arts <- est$by_occ$all_mgmt_bus_science_arts / 
    est$by_occ$all_
  est$by_occ$pct_service <- est$by_occ$all_service / 
    est$by_occ$all_
  est$by_occ$pct_sales_office <- est$by_occ$all_sales_office / 
    est$by_occ$all_
  est$by_occ$pct_resources_constr_maint <- est$by_occ$all_resources_constr_maint / 
    est$by_occ$all_
  est$by_occ$pct_production_transport <- est$by_occ$all_production_transport / 
    est$by_occ$all_
  est$by_occ$pct_mil <- est$by_occ$all_mil / 
    est$by_occ$all_
  
  se$by_occ$pct_mgmt_bus_science_arts <- sqrt(se$by_occ$all_mgmt_bus_science_arts^2 - (
    est$by_occ$pct_mgmt_bus_science_arts^2 * se$by_occ$all_^2)) / 
    est$by_occ$all_
  se$by_occ$pct_service <- sqrt(se$by_occ$all_service^2 - (
    est$by_occ$pct_service^2 * se$by_occ$all_^2)) / 
    est$by_occ$all_
  se$by_occ$pct_sales_office <- sqrt(se$by_occ$all_sales_office^2 - (
    est$by_occ$pct_sales_office^2 * se$by_occ$all_^2)) / 
    est$by_occ$all_
  se$by_occ$pct_resources_constr_maint <- sqrt(se$by_occ$all_resources_constr_maint^2 - (
    est$by_occ$pct_resources_constr_maint^2 * se$by_occ$all_^2)) / 
    est$by_occ$all_
  se$by_occ$pct_production_transport <- sqrt(se$by_occ$all_production_transport^2 - (
    est$by_occ$pct_production_transport^2 * se$by_occ$all_^2)) / 
    est$by_occ$all_
  se$by_occ$pct_mil <- sqrt(se$by_occ$all_mil^2 - (est$by_occ$pct_mil^2 * se$by_occ$all_^2)) / 
    est$by_occ$all_
  
  ## 07 -- (A) combine columns and (B) calc percentages
  ### msa_place
  #----------------------------------------------
  est$msa_place <- est$msa_place[, -c(4,5,7,8,10,11, 15,16,18,19,21,22)]
  se$msa_place  <- se$msa_place[, -c(4,5,7,8,10,11, 15,16,18,19,21,22)]
  names(est$msa_place) <- names(se$msa_place) <- paste("msa", c("cnt_all", paste(
    rep(c("live_city", "live_out_city"), each= 5),
    rep(c("all", "work_same_msa", "work_diff_msa", "work_microsa", "work_rural"), 2), sep= "_")), sep="_")
  
  
  ## 08 -- (A) combine columns and (B) calc percentages
  ### microsa_place
  #----------------------------------------------
  est$microsa_place <- est$microsa_place[, -c(4,5,7,8,10,11, 15,16,18,19,21,22)]
  se$microsa_place  <- se$microsa_place[, -c(4,5,7,8,10,11, 15,16,18,19,21,22)]
  names(est$microsa_place) <- names(se$microsa_place) <- paste("microsa", c("cnt_all", paste(
    rep(c("live_city", "live_out_city"), each= 5),
    rep(c("all", "work_same_microsa", "work_diff_microsa", "work_msa", "work_rural"), 2), sep= "_")), sep="_")
  
  # 09 -- sort, combine, and return
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
  names(ret$estimates) <- names(ret$standard_error) <- c("time_to_work_by_sex", "mode_transit_by_age",
    "median_earnings_by_mode_transit", "median_age_by_mode_transit", "mode_transit_by_occ", 
    "place_of_work_metroSA", "place_of_work_microSA")
  
  return(ret)
}
