
#' @title Pull ACS population data
#' @description Pull ACS data for a specified geography from base tables
#' B01001, B01002, B02001, B06007, B06008, B06009, B06010, B06011, AND B06012.
#' These tables reference population counts by a number of slices.
#' Multiple additional fields, mainly percentages and aggregations, are calculated.
#' @param endyear An integer, indicating the latest year of the data in the survey.
#' @param span An integer in \code{c(1,3,5)} indicating the span of the desired data.
#' @param geography a valid \code{geo.set} object specifying the census geography or 
#' geographies to be fetched.
#' @return A \code{list} containing the endyear, span, a \code{data.frame} of estimates,
#' a \code{data.frame} of standard errors, a character vector of the original column names,
#' and a \code{data.frame} of the geography metadata from \code{\link[acs]{acs.fetch}}.
#' @seealso \code{\link[acs]{acs.fetch}}, \code{\link[acs]{geo.make}}
#' @export
pull_population <- function(endyear, span, geography) {
  # 00 -- error checking
  #----------------------------------------------
  check_geo_inputs(endyear= endyear, span= span, geography= geography)
  
  # 01 -- pull data
  #----------------------------------------------
  oldw <- getOption("warn")
  options(warn= -1) # suppress warnings from library(acs) / ACS API
  on.exit(options(warn= oldw)) # turn warnings back on
  age_by_sex <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                          table.number = "B01001", col.names= "pretty")
  med_age <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                       table.number = "B01002", col.names = "pretty")
  pop_by_race <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                           table.number = "B02001", col.names = "pretty")
  birth_and_language <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                                  table.number = "B06007", col.names = "pretty")
  by_marital_status <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                                 table.number = "B06008", col.names = "pretty")
  by_edu <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                      table.number = "B06009", col.names = "pretty")
  by_inc_12mo <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                           table.number = "B06010", col.names = "pretty")
  med_inc_12mo <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                            table.number = "B06011", col.names = "pretty")
  by_pov_status <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                             table.number = "B06012", col.names = "pretty")
  
  # 02 -- create lists of EST and SE -- as data.frames
  #----------------------------------------------
  est <- list(age_by_sex= data.frame(age_by_sex@estimate),
              med_age= data.frame(med_age@estimate),
              pop_by_race= data.frame(pop_by_race@estimate[, c(1:7), drop= FALSE]),
              birth_and_language= data.frame(birth_and_language@estimate),
              by_marital_status= data.frame(by_marital_status@estimate),
              by_edu= data.frame(by_edu@estimate),
              by_inc_12mo= data.frame(by_inc_12mo@estimate),
              med_inc_12mo= data.frame(med_inc_12mo@estimate),
              by_pov_status= data.frame(by_pov_status@estimate))
  
  se <- list(age_by_sex= data.frame(age_by_sex@standard.error),
             med_age= data.frame(med_age@standard.error),
             pop_by_race= data.frame(pop_by_race@standard.error[, c(1:7), drop= FALSE]),
             birth_and_language= data.frame(birth_and_language@standard.error),
             by_marital_status= data.frame(by_marital_status@standard.error),
             by_edu= data.frame(by_edu@standard.error),
             by_inc_12mo= data.frame(by_inc_12mo@standard.error),
             med_inc_12mo= data.frame(med_inc_12mo@standard.error),
             by_pov_status= data.frame(by_pov_status@standard.error))
  
  geo <- age_by_sex@geography
  
  rm(age_by_sex, med_age, pop_by_race, birth_and_language, by_marital_status,
     by_edu, by_inc_12mo, by_pov_status, med_inc_12mo)
  
  # 03 -- combine columns
  #----------------------------------------------
  ### birth and language
  est$birth_and_language$born_us                     <- apply(est$birth_and_language[, c(9,17)], 1, sum)
  est$birth_and_language$born_us_only_eng            <- apply(est$birth_and_language[, c(10,18)], 1, sum)
  est$birth_and_language$born_us_spk_span            <- apply(est$birth_and_language[, c(11,19)], 1, sum)
  est$birth_and_language$born_us_spk_span_eng_vw     <- apply(est$birth_and_language[, c(12,20)], 1, sum)
  est$birth_and_language$born_us_spk_span_eng_lt_vw  <- apply(est$birth_and_language[, c(13,21)], 1, sum)
  est$birth_and_language$born_us_spk_other           <- apply(est$birth_and_language[, c(14,22)], 1, sum)
  est$birth_and_language$born_us_spk_other_eng_vw    <- apply(est$birth_and_language[, c(15,23)], 1, sum)
  est$birth_and_language$born_us_spk_other_eng_lt_vw <- apply(est$birth_and_language[, c(16,24)], 1, sum)
  est$birth_and_language <- est$birth_and_language[, c(1:8,25:48)]
  
  se$birth_and_language$born_us                     <- sqrt(apply(se$birth_and_language[, c(9,17)]^2, 1, sum))
  se$birth_and_language$born_us_only_eng            <- sqrt(apply(se$birth_and_language[, c(10,18)]^2, 1, sum))
  se$birth_and_language$born_us_spk_span            <- sqrt(apply(se$birth_and_language[, c(11,19)]^2, 1, sum))
  se$birth_and_language$born_us_spk_span_eng_vw     <- sqrt(apply(se$birth_and_language[, c(12,20)]^2, 1, sum))
  se$birth_and_language$born_us_spk_span_eng_lt_vw  <- sqrt(apply(se$birth_and_language[, c(13,21)]^2, 1, sum))
  se$birth_and_language$born_us_spk_other           <- sqrt(apply(se$birth_and_language[, c(14,22)]^2, 1, sum))
  se$birth_and_language$born_us_spk_other_eng_vw    <- sqrt(apply(se$birth_and_language[, c(15,23)]^2, 1, sum))
  se$birth_and_language$born_us_spk_other_eng_lt_vw <- sqrt(apply(se$birth_and_language[, c(16,24)]^2, 1, sum))
  se$birth_and_language <- se$birth_and_language[, c(1:8,25:48)]
  
  names(se$birth_and_language) <- names(est$birth_and_language) <- paste(
    rep(c("all", "citizen_born_out_us", "foreign_born", "us_born"), each= 8),
    rep(c("total", "only_eng", "spk_span", "spk_span_eng_vw", "spk_span_eng_lt_vw",
          "other", "other_eng_vw", "other_eng_lt_vw"), 4), sep= "_")
  
  ### by_marital_status
  est$by_marital_status$born_us               <- apply(est$by_marital_status[, c(7,13)], 1, sum)
  est$by_marital_status$born_us_never_married <- apply(est$by_marital_status[, c(8,14)], 1, sum)
  est$by_marital_status$born_us_divorced      <- apply(est$by_marital_status[, c(9,16)], 1, sum)
  est$by_marital_status$born_us_separated     <- apply(est$by_marital_status[, c(10,17)], 1, sum)
  est$by_marital_status$born_us_widowed       <- apply(est$by_marital_status[, c(11,18)], 1, sum)
  est$by_marital_status <- est$by_marital_status[, c(1,2,4:6,19,20,22:24,25,26,28:35)]
  
  se$by_marital_status$born_us               <- sqrt(apply(se$by_marital_status[, c(7,13)]^2, 1, sum))
  se$by_marital_status$born_us_never_married <- sqrt(apply(se$by_marital_status[, c(8,14)]^2, 1, sum))
  se$by_marital_status$born_us_divorced      <- sqrt(apply(se$by_marital_status[, c(9,16)]^2, 1, sum))
  se$by_marital_status$born_us_separated     <- sqrt(apply(se$by_marital_status[, c(10,17)]^2, 1, sum))
  se$by_marital_status$born_us_widowed       <- sqrt(apply(se$by_marital_status[, c(11,18)]^2, 1, sum))
  se$by_marital_status <- se$by_marital_status[, c(1,2,4:6,19,20,22:24,25,26,28:35)]
  
  names(se$by_marital_status) <- names(est$by_marital_status) <- paste(
    rep(c("all", "citizen_born_out_us", "foreign_born", "us_born"), each= 5),
    rep(c("total", "nvr_married", "divorced", "separated", "widowed"), 4), sep= "_")
  
  ### by edu
  est$by_edu$born_us           <- apply(est$by_edu[, c(7,13)], 1, sum)
  est$by_edu$born_us_lt_hs     <- apply(est$by_edu[, c(8,14)], 1, sum)
  est$by_edu$born_us_hs        <- apply(est$by_edu[, c(9,15)], 1, sum)
  est$by_edu$born_us_some_col  <- apply(est$by_edu[, c(10,16)], 1, sum)
  est$by_edu$born_us_bachelors <- apply(est$by_edu[, c(11,17)], 1, sum)
  est$by_edu$born_us_graduate  <- apply(est$by_edu[, c(12,18)], 1, sum)
  est$by_edu <- est$by_edu[, c(1:6,19:36)]
  
  se$by_edu$born_us           <- sqrt(apply(se$by_edu[, c(7,13)]^2, 1, sum))
  se$by_edu$born_us_lt_hs     <- sqrt(apply(se$by_edu[, c(8,14)]^2, 1, sum))
  se$by_edu$born_us_hs        <- sqrt(apply(se$by_edu[, c(9,15)]^2, 1, sum))
  se$by_edu$born_us_some_col  <- sqrt(apply(se$by_edu[, c(10,16)]^2, 1, sum))
  se$by_edu$born_us_bachelors <- sqrt(apply(se$by_edu[, c(11,17)]^2, 1, sum))
  se$by_edu$born_us_graduate  <- sqrt(apply(se$by_edu[, c(12,18)]^2, 1, sum))
  se$by_edu <- se$by_edu[, c(1:6,19:36)]
  
  names(est$by_edu) <- names(se$by_edu) <- paste(
    rep(c("all", "citizen_born_out_us", "foreign_born", "us_born"), each= 6),
    rep(c("total", "lt_hs", "hs", "some_col", "bachelors", "graduate"), 4), sep= "_")
  
  
  ### by inc 12 mo
  est$by_inc_12mo$born_us            <- apply(est$by_inc_12mo[, c(12,23)], 1, sum)
  est$by_inc_12mo$born_us_no_inc     <- apply(est$by_inc_12mo[, c(13,24)], 1, sum)
  est$by_inc_12mo$born_us_1_lt10k    <- apply(est$by_inc_12mo[, c(15,26)], 1, sum)
  est$by_inc_12mo$born_us_10k_lt15k  <- apply(est$by_inc_12mo[, c(16,27)], 1, sum)
  est$by_inc_12mo$born_us_15k_lt25k  <- apply(est$by_inc_12mo[, c(17,28)], 1, sum)
  est$by_inc_12mo$born_us_25k_lt35k  <- apply(est$by_inc_12mo[, c(18,29)], 1, sum)
  est$by_inc_12mo$born_us_35k_lt50k  <- apply(est$by_inc_12mo[, c(19,30)], 1, sum)
  est$by_inc_12mo$born_us_50k_lt_65k <- apply(est$by_inc_12mo[, c(20,31)], 1, sum)
  est$by_inc_12mo$born_us_65k_lt75k  <- apply(est$by_inc_12mo[, c(21,32)], 1, sum)
  est$by_inc_12mo$born_us_gt75k      <- apply(est$by_inc_12mo[, c(22,33)], 1, sum)
  est$by_inc_12mo <- est$by_inc_12mo[, c(1:2,4:11,34,35,37:46,48:65)]
  
  se$by_inc_12mo$born_us            <- sqrt(apply(se$by_inc_12mo[, c(12,23)]^2, 1, sum))
  se$by_inc_12mo$born_us_no_inc     <- sqrt(apply(se$by_inc_12mo[, c(13,24)]^2, 1, sum))
  se$by_inc_12mo$born_us_1_lt10k    <- sqrt(apply(se$by_inc_12mo[, c(15,26)]^2, 1, sum))
  se$by_inc_12mo$born_us_10k_lt15k  <- sqrt(apply(se$by_inc_12mo[, c(16,27)]^2, 1, sum))
  se$by_inc_12mo$born_us_15k_lt25k  <- sqrt(apply(se$by_inc_12mo[, c(17,28)]^2, 1, sum))
  se$by_inc_12mo$born_us_25k_lt35k  <- sqrt(apply(se$by_inc_12mo[, c(18,29)]^2, 1, sum))
  se$by_inc_12mo$born_us_35k_lt50k  <- sqrt(apply(se$by_inc_12mo[, c(19,30)]^2, 1, sum))
  se$by_inc_12mo$born_us_50k_lt_65k <- sqrt(apply(se$by_inc_12mo[, c(20,31)]^2, 1, sum))
  se$by_inc_12mo$born_us_65k_lt75k  <- sqrt(apply(se$by_inc_12mo[, c(21,32)]^2, 1, sum))
  se$by_inc_12mo$born_us_gt75k      <- sqrt(apply(se$by_inc_12mo[, c(22,33)]^2, 1, sum))
  se$by_inc_12mo <- se$by_inc_12mo[, c(1:2,4:11,34,35,37:46,48:65)]
  
  names(se$by_inc_12mo) <- names(est$by_inc_12mo) <- paste(
    rep(c("all", "citizen_born_out_us", "foreign_born", "us_born"), each= 10),
    rep(c("total", "no_inc", "1_lt10k", "10k_lt15k", "15k_lt25k", "25k_lt35k", "35k_lt50k", 
          "50k_lt_65k", "65k_lt75k", "gt75k"), 4), sep= "_")
  
  ### by pov status
  est$by_pov_status$bus1 <- apply(est$by_pov_status[, c(5,9)], 1, sum)
  est$by_pov_status$bus2 <- apply(est$by_pov_status[, c(6,10)], 1, sum)
  est$by_pov_status$bus3 <- apply(est$by_pov_status[, c(7,11)], 1, sum)
  est$by_pov_status$bus4 <- apply(est$by_pov_status[, c(8,12)], 1, sum)
  est$by_pov_status <- est$by_pov_status[, c(1:4,13:24)]
  
  se$by_pov_status$bus1 <- sqrt(apply(se$by_pov_status[, c(5,9)]^2, 1, sum))
  se$by_pov_status$bus2 <- sqrt(apply(se$by_pov_status[, c(6,10)]^2, 1, sum))
  se$by_pov_status$bus3 <- sqrt(apply(se$by_pov_status[, c(7,11)]^2, 1, sum))
  se$by_pov_status$bus4 <- sqrt(apply(se$by_pov_status[, c(8,12)]^2, 1, sum))
  se$by_pov_status <- se$by_pov_status[, c(1:4,13:24)]
  
  names(se$by_pov_status) <- names(est$by_pov_status) <- paste(
    rep(c("all", "citizen_born_out_us", "foreign_born", "us_born"), each= 4),
    rep(c("total", "0_lt100pct_pov_lvl", "100_lt150pct_pov_lvl", "gt150pct_pov_lvl"), 4), sep= "_")
  
  names(est$age_by_sex) <- names(se$age_by_sex) <- c("total", paste(
    rep(c("m", "f"), each= 24),
    rep(c("", "u5", "5_9", "10_14", "15_17", "18_19", "20", "21", "22_24", "25_29", "30_34",
          "35_39", "40_44", "45_49", "50_54", "55_59", "60_61", "62_64", "65_66", "67_69",
          "70_74", "75_79", "80_84", "85up"), 2), sep= "_"))
  
  names(est$pop_by_race) <- names(se$pop_by_race) <-  c("total", "white", "black_AA",
    "nat_amer", "asian", "pac_isl", "other")
  
  names(est$med_age) <- names(se$med_age) <- c("total", "male", "female")
  
  # 04 -- calc percentages
  #----------------------------------------------
  est$pop_by_race$pct_white <- est$pop_by_race$white / est$pop_by_race$total
  est$pop_by_race$pct_black <- est$pop_by_race$black_AA / est$pop_by_race$total
  est$pop_by_race$pct_nat_amer <- est$pop_by_race$nat_amer / est$pop_by_race$total
  est$pop_by_race$pct_asian <- est$pop_by_race$asian / est$pop_by_race$total
  est$pop_by_race$pct_pac_isl <- est$pop_by_race$pac_isl / est$pop_by_race$total
  est$pop_by_race$pct_other <- est$pop_by_race$other / est$pop_by_race$total
  
  se$pop_by_race$pct_white <- sqrt(se$pop_by_race$white^2 - (
    (est$pop_by_race$white / est$pop_by_race$total)^2 * se$pop_by_race$total^2)) / est$pop_by_race$total
  
  se$pop_by_race$pct_black <- sqrt(se$pop_by_race$black_AA^2 - (
    (est$pop_by_race$black_AA / est$pop_by_race$total)^2 * se$pop_by_race$total^2)) / est$pop_by_race$total
  
  se$pop_by_race$pct_nat_amer <- sqrt(se$pop_by_race$nat_amer^2 - (
    (est$pop_by_race$nat_amer / est$pop_by_race$total)^2 * se$pop_by_race$total^2)) / est$pop_by_race$total
  
  se$pop_by_race$pct_asian <- sqrt(se$pop_by_race$asian^2 - (
    (est$pop_by_race$asian / est$pop_by_race$total)^2 * se$pop_by_race$total^2)) / est$pop_by_race$total
  
  se$pop_by_race$pct_pac_isl <- sqrt(se$pop_by_race$pac_isl^2 - (
    (est$pop_by_race$pac_isl / est$pop_by_race$total)^2 * se$pop_by_race$total^2)) / est$pop_by_race$total
  
  se$pop_by_race$pct_other <- sqrt(se$pop_by_race$other^2 - (
    (est$pop_by_race$other / est$pop_by_race$total)^2 * se$pop_by_race$total^2)) / est$pop_by_race$total
  
  # 05 -- sort, combine, and return
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
  names(ret$estimates) <- names(ret$standard_error) <- c("sex_by_age", "median_age_by_sex", "pop_by_race",
    "place_birth_by_lang_at_home", "place_birth_by_mar_status", "place_birth_by_edu_attain",
    "place_birth_by_income", "median_income_by_place_birth", "place_birth_by_pov_status")
  
  return(ret)
}
