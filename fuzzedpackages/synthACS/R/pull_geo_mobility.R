

#' @title Pull ACS geographic mobility data
#' @description Pull ACS data for a specified geography from base tables
#' B07001, B07003, B07008, B07009, B07010, and B07012.
#' These tables provide data on geographic mobility in the past year by a number of slices.
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
pull_geo_mobility <- function(endyear, span, geography) {
  # 00 -- error checking
  #----------------------------------------------
  check_geo_inputs(endyear= endyear, span= span, geography= geography)
  
  # 01 -- pull data
  #----------------------------------------------
  oldw <- getOption("warn")
  options(warn= -1) # suppress warnings from library(acs) / ACS API
  on.exit(options(warn= oldw)) # turn warnings back on
  geo_by_age <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                          table.number = "B07001", col.names= "pretty")
  
  geo_by_sex <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                          table.number = "B07003", col.names= "pretty")
  
  geo_by_mar <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                          table.number = "B07008", col.names= "pretty")
  
  geo_by_edu <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                          table.number = "B07009", col.names= "pretty")
  
  geo_by_inc <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                          table.number = "B07010", col.names= "pretty")
  
  est <- list(by_age= data.frame(geo_by_age@estimate[, c(1,17:dim(geo_by_age@estimate)[[2]]), drop= FALSE]),
              by_sex= data.frame(geo_by_sex@estimate[, c(1,4:dim(geo_by_sex@estimate)[[2]]), drop= FALSE]),
              by_mar= data.frame(geo_by_mar@estimate[, c(1,7:dim(geo_by_mar@estimate)[[2]]), drop= FALSE]),
              by_edu= data.frame(geo_by_edu@estimate[, c(1,7:dim(geo_by_edu@estimate)[[2]]), drop= FALSE]),
              by_inc= data.frame(geo_by_inc@estimate[, c(1,12:dim(geo_by_inc@estimate)[[2]]), drop= FALSE]))
  
  se <- list(by_age= data.frame(geo_by_age@standard.error[, c(1,17:dim(geo_by_age@standard.error)[[2]]), 
                                                          drop= FALSE]),
             by_sex= data.frame(geo_by_sex@standard.error[, c(1,4:dim(geo_by_sex@standard.error)[[2]]), 
                                                          drop= FALSE]),
             by_mar= data.frame(geo_by_mar@standard.error[, c(1,7:dim(geo_by_mar@standard.error)[[2]]), 
                                                          drop= FALSE]),
             by_edu= data.frame(geo_by_edu@standard.error[, c(1,7:dim(geo_by_edu@standard.error)[[2]]), 
                                                          drop= FALSE]),
             by_inc= data.frame(geo_by_inc@standard.error[, c(1,12:dim(geo_by_inc@standard.error)[[2]]), 
                                                          drop= FALSE]))
  
  geo <- geo_by_age@geography
  
  rm(geo_by_age, geo_by_sex, geo_by_mar, geo_by_edu, geo_by_inc)
  
  ## 02 -- (A) combine columns and (B) calc percentages
  ### by age
  #----------------------------------------------
  est$by_age$same_house_18_34 <- apply(est$by_age[, 5:8], 1, sum)
  est$by_age$same_house_35_54 <- apply(est$by_age[, 9:12], 1, sum)
  est$by_age$same_house_55_64 <- apply(est$by_age[, 13:14], 1, sum)
  est$by_age$same_house_65up  <- apply(est$by_age[, 15:17], 1, sum)
  est$by_age$same_cnty_18_34  <- apply(est$by_age[, 21:24], 1, sum)
  est$by_age$same_cnty_35_54  <- apply(est$by_age[, 25:28], 1, sum)
  est$by_age$same_cnty_55_64  <- apply(est$by_age[, 29:30], 1, sum)
  est$by_age$same_cnty_65up   <- apply(est$by_age[, 31:33], 1, sum)
  est$by_age$same_st_18_34    <- apply(est$by_age[, 37:40], 1, sum)
  est$by_age$same_st_35_54    <- apply(est$by_age[, 41:44], 1, sum)
  est$by_age$same_st_55_64    <- apply(est$by_age[, 45:46], 1, sum)
  est$by_age$same_st_65up     <- apply(est$by_age[, 47:49], 1, sum)
  est$by_age$diff_st_18_34    <- apply(est$by_age[, 53:56], 1, sum)
  est$by_age$diff_st_35_54    <- apply(est$by_age[, 57:60], 1, sum)
  est$by_age$diff_st_55_64    <- apply(est$by_age[, 61:62], 1, sum)
  est$by_age$diff_st_65up     <- apply(est$by_age[, 63:65], 1, sum)
  est$by_age$abroad_18_34     <- apply(est$by_age[, 69:72], 1, sum)
  est$by_age$abroad_35_54     <- apply(est$by_age[, 73:76], 1, sum)
  est$by_age$abroad_55_64     <- apply(est$by_age[, 77:78], 1, sum)
  est$by_age$abroad_65up      <- apply(est$by_age[, 79:81], 1, sum)
  est$by_age <- est$by_age[, c(1:4, 82:85, 18:20, 86:89, 34:36, 90:93, 50:52, 94:97, 66:68, 98:101)]
  
  se$by_age$same_house_18_34 <- sqrt(apply(se$by_age[, 5:8]^2, 1, sum))
  se$by_age$same_house_35_54 <- sqrt(apply(se$by_age[, 9:12]^2, 1, sum))
  se$by_age$same_house_55_64 <- sqrt(apply(se$by_age[, 13:14]^2, 1, sum))
  se$by_age$same_house_65up  <- sqrt(apply(se$by_age[, 15:17]^2, 1, sum))
  se$by_age$same_cnty_18_34  <- sqrt(apply(se$by_age[, 21:24]^2, 1, sum))
  se$by_age$same_cnty_35_54  <- sqrt(apply(se$by_age[, 25:28]^2, 1, sum))
  se$by_age$same_cnty_55_64  <- sqrt(apply(se$by_age[, 29:30]^2, 1, sum))
  se$by_age$same_cnty_65up   <- sqrt(apply(se$by_age[, 31:33]^2, 1, sum))
  se$by_age$same_st_18_34    <- sqrt(apply(se$by_age[, 37:40]^2, 1, sum))
  se$by_age$same_st_35_54    <- sqrt(apply(se$by_age[, 41:44]^2, 1, sum))
  se$by_age$same_st_55_64    <- sqrt(apply(se$by_age[, 45:46]^2, 1, sum))
  se$by_age$same_st_65up     <- sqrt(apply(se$by_age[, 47:49]^2, 1, sum))
  se$by_age$diff_st_18_34    <- sqrt(apply(se$by_age[, 53:56]^2, 1, sum))
  se$by_age$diff_st_35_54    <- sqrt(apply(se$by_age[, 57:60]^2, 1, sum))
  se$by_age$diff_st_55_64    <- sqrt(apply(se$by_age[, 61:62]^2, 1, sum))
  se$by_age$diff_st_65up     <- sqrt(apply(se$by_age[, 63:65]^2, 1, sum))
  se$by_age$abroad_18_34     <- sqrt(apply(se$by_age[, 69:72]^2, 1, sum))
  se$by_age$abroad_35_54     <- sqrt(apply(se$by_age[, 73:76]^2, 1, sum))
  se$by_age$abroad_55_64     <- sqrt(apply(se$by_age[, 77:78]^2, 1, sum))
  se$by_age$abroad_65up      <- sqrt(apply(se$by_age[, 79:81]^2, 1, sum))
  se$by_age <- se$by_age[, c(1:4, 82:85, 18:20, 86:89, 34:36, 90:93, 50:52, 94:97, 66:68, 98:101)]
  
  names(est$by_age) <- names(se$by_age) <- c("total", paste(
    rep(c("same_house", "same_cnty", "same_st", "diff_st", "abroad"), each= 7),
    rep(c("all", "1_4", "5_17", "18_34", "35_54", "55_64", "65up"), 5), sep= "_"))
  
  ### pct's
  est$by_age$pct_same_house <- est$by_age$same_house_all / est$by_age$total
  est$by_age$pct_same_cnty <- est$by_age$same_cnty_all / est$by_age$total
  est$by_age$pct_same_st <- est$by_age$same_st_all / est$by_age$total
  est$by_age$pct_diff_st <- est$by_age$diff_st_all / est$by_age$total
  est$by_age$pct_abroad <- est$by_age$abroad_all / est$by_age$total
  
  se$by_age$pct_same_house <- sqrt(se$by_age$same_house_all^2 - (
    (est$by_age$same_house_all / est$by_age$total)^2 * se$by_age$total^2)) / est$by_age$total
  
  se$by_age$pct_same_cnty <- sqrt(se$by_age$same_cnty_all^2 - (
    (est$by_age$same_cnty_all / est$by_age$total)^2 * se$by_age$total^2)) / est$by_age$total
  
  se$by_age$pct_same_st <- sqrt(se$by_age$same_st_all^2 - (
    (est$by_age$same_st_all / est$by_age$total)^2 * se$by_age$total^2)) / est$by_age$total
  
  se$by_age$pct_diff_st <- sqrt(se$by_age$diff_st_all^2 - (
    (est$by_age$diff_st_all / est$by_age$total)^2 * se$by_age$total^2)) / est$by_age$total
  
  se$by_age$pct_abroad <- sqrt(se$by_age$abroad_all^2 - (
    (est$by_age$abroad_all / est$by_age$total)^2 * se$by_age$total^2)) / est$by_age$total
  
  ## 03 -- (A) combine columns and (B) calc percentages
  ### by sex
  #----------------------------------------------
  names(est$by_sex) <- names(se$by_sex) <- c("total", paste(
    rep(c("same_house", "same_cnty", "same_st", "diff_st", "abroad"), each= 3),
    rep(c("all", "male", "female"), 5), sep= "_"))
  
  ### pcts
  est$by_sex$pctof_same_house_f <- est$by_sex$same_house_female / est$by_sex$same_house_all
  est$by_sex$pctof_same_cnty_f  <- est$by_sex$same_cnty_female / est$by_sex$same_cnty_all
  est$by_sex$pctof_same_st_f    <- est$by_sex$same_st_female / est$by_sex$same_st_all
  est$by_sex$pctof_diff_st_f    <- est$by_sex$diff_st_female / est$by_sex$diff_st_all
  est$by_sex$pctof_abroad_f     <- est$by_sex$abroad_female / est$by_sex$abroad_all
  
  se$by_sex$pctof_same_house_f <- sqrt(se$by_sex$same_house_female^2 - (
    (est$by_sex$same_house_female / est$by_sex$same_house_all)^2 * se$by_sex$same_house_all^2)) /
    est$by_sex$same_house_all
  
  se$by_sex$pctof_same_cnty_f <- sqrt(se$by_sex$same_cnty_female^2 - (
    (est$by_sex$same_cnty_female / est$by_sex$same_cnty_all)^2 * se$by_sex$same_cnty_all^2)) /
    est$by_sex$same_cnty_all
  
  se$by_sex$pctof_same_st_f <- sqrt(se$by_sex$same_st_female^2 - (
    (est$by_sex$same_st_female / est$by_sex$same_st_all)^2 * se$by_sex$same_st_all^2)) /
    est$by_sex$same_st_all
  
  se$by_sex$pctof_diff_st_f <- sqrt(se$by_sex$diff_st_female^2 - (
    (est$by_sex$diff_st_female / est$by_sex$diff_st_all)^2 * se$by_sex$diff_st_all^2)) /
    est$by_sex$diff_st_all
  
  se$by_sex$pctof_abroad_f <- sqrt(se$by_sex$abroad_female^2 - (
    (est$by_sex$abroad_female / est$by_sex$abroad_all)^2 * se$by_sex$abroad_all^2)) /
    est$by_sex$abroad_all

  ## 04 -- (A) combine columns and (B) calc percentages
  ### by marital status
  #----------------------------------------------
  names(est$by_mar) <- names(se$by_mar) <- c("total", paste(
    rep(c("same_house", "same_cnty", "same_st", "diff_st", "abroad"), each= 6),
    rep(c("all", "nvr_married", "now_mar_exc_sep", "divorced", "separated",
          "widowed"), 5), sep= "_"))
  
  ## 05 -- (A) combine columns and (B) calc percentages
  ### by edu
  #----------------------------------------------
  names(est$by_edu) <- names(se$by_edu) <- c("total", paste(
    rep(c("same_house", "same_cnty", "same_st", "diff_st", "abroad"), each= 6),
    rep(c("all", "lt_hs", "hs", "some_col", "bachelors", "graduate"), 5), sep= "_"))
  
  
  ## 06 -- (A) combine columns and (B) calc percentages
  ### by inc
  #----------------------------------------------
  names(est$by_inc) <- names(se$by_inc) <- c("total", paste(
    rep(c("same_house", "same_cnty", "same_st", "diff_st", "abroad"), each= 11),
    rep(c("all", "no_inc", "inc", "1_lt10k", "10k_lt15k", "15k_lt_25k", "25k_lt35k",
          "35k_lt50k", "50k_lt_65k", "65k_lt_75k", "75k_up"), 5), sep= "_"))
  est$by_inc <- est$by_inc[, -c(4,15,26,37,48)]
  se$by_inc  <- se$by_inc[, -c(4,15,26,37,48)]
  
  # 07 -- sort, combine, and return
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
  names(ret$estimates) <- names(ret$standard_error) <- c("geo_mob_by_age", "geo_mob_by_sex",
    "geo_mob_by_mar_status", "geo_mob_by_edu_attain", "geo_mob_by_income") # , "geo_mob_by_pov_status")
  
  return(ret)
}

