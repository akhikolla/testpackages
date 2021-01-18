
#' @title Pull ACS income and earnings data
#' @description Pull ACS data for a specified geography from base tables
#' B17001, B17004, B18101, B19001, B19013, B19055, B19057. Not yet 
#' implemented: B17002
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
pull_pov_inc <- function(endyear, span, geography) {
  # 00 -- error checking
  #----------------------------------------------
  check_geo_inputs(endyear= endyear, span= span, geography= geography)
  
  # 01 -- pull data and move to lists
  #----------------------------------------------
  oldw <- getOption("warn")
  options(warn= -1) # suppress warnings from library(acs) / ACS API
  on.exit(options(warn= oldw)) # turn warnings back on
  pov_status <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                    table.number = "B17001", col.names= "pretty")
#   pov_to_inc <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
#                     table.number = "B17002", col.names= "pretty")
  pov_to_work <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                    table.number = "B17004", col.names= "pretty")
  dis_by_age <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                    table.number = "B18101", col.names= "pretty")
  hh_inc <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                    table.number = "B19001", col.names= "pretty")
  med_hh_inc <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                      table.number = "B19013", col.names= "pretty")
  hh_ss_inc <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                      table.number = "B19055", col.names= "pretty")
  hh_psa_inc <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                      table.number = "B19057", col.names= "pretty")
  
  est <- list(pov_status= data.frame(pov_status@estimate),
              #pov_to_inc= data.frame(pov_to_inc@estimate),
              pov_to_work= data.frame(t(pov_to_work@estimate[, -c(3,7,12,16)])),
              dis_by_age= data.frame(dis_by_age@estimate),
              hh_inc= data.frame(hh_inc@estimate),
              med_hh_inc= data.frame(med_hh_inc@estimate),
              hh_ss_inc= data.frame(hh_ss_inc@estimate),
              hh_psa_inc= data.frame(hh_psa_inc@estimate))
  
  se <- list(pov_status= data.frame(pov_status@standard.error),
              #pov_to_inc= data.frame(pov_to_inc@standard.error),
              pov_to_work= data.frame(t(pov_to_work@standard.error[, -c(3,7,12,16)])),
              dis_by_age= data.frame(dis_by_age@standard.error),
              hh_inc= data.frame(hh_inc@standard.error),
              med_hh_inc= data.frame(med_hh_inc@standard.error),
              hh_ss_inc= data.frame(hh_ss_inc@standard.error),
              hh_psa_inc= data.frame(hh_psa_inc@standard.error))
  
  geo <- pov_status@geography
  
  ## 02 -- (A) combine columns and (B) calc percentages
  ### pov_status
  #----------------------------------------------
  est$pov_status$m_below_pov_u18 <- apply(est$pov_status[, 4:9], 1, sum)
  est$pov_status$m_above_pov_u18 <- apply(est$pov_status[, 33:38], 1, sum)
  est$pov_status$f_below_pov_u18 <- apply(est$pov_status[, 18:23], 1, sum)
  est$pov_status$f_above_pov_u18 <- apply(est$pov_status[, 47:52], 1, sum)
  est$pov_status$pct_below_pov   <- est$pov_status[,2] / est$pov_status[,1]
  est$pov_status$pct_below_pov_m   <- est$pov_status[,3] / apply(est$pov_status[,c(3,32)], 1, sum)
  est$pov_status$pct_below_pov_f   <- est$pov_status[,17] / apply(est$pov_status[,c(17,46)],1,sum)
  
  se$pov_status$m_below_pov_u18 <- sqrt(apply(se$pov_status[, 4:9]^2, 1, sum))
  se$pov_status$m_above_pov_u18 <- sqrt(apply(se$pov_status[, 33:38]^2, 1, sum))
  se$pov_status$f_below_pov_u18 <- sqrt(apply(se$pov_status[, 18:23]^2, 1, sum))
  se$pov_status$f_above_pov_u18 <- sqrt(apply(se$pov_status[, 47:52]^2, 1, sum))
  se$pov_status$pct_below_pov   <- sqrt(se$pov_status[,2]^2 - (
    est$pov_status$pct_below_pov^2 * se$pov_status[,1]^2)) / est$pov_status[,1]
  se$pov_status$pct_below_pov_m   <- sqrt(se$pov_status[,3]^2 - (
    est$pov_status$pct_below_pov_m^2 * apply(se$pov_status[,c(3,32)]^2, 1, sum))) / 
    apply(est$pov_status[,c(3,32)], 1, sum)
  se$pov_status$pct_below_pov_f   <- sqrt(se$pov_status[,17]^2 - (
    est$pov_status$pct_below_pov_f^2 * apply(se$pov_status[,c(17,46)]^2,1,sum))) / 
    apply(est$pov_status[,c(17,46)],1,sum)
  
  est$pov_status <- est$pov_status[, c(1,3,32,17,46,64:66,60,10:16,61,24:30,61,39:45,62,53:59)]
  se$pov_status  <- se$pov_status[, c(1,3,32,17,46,64:66,60,10:16,61,24:30,61,39:45,62,53:59)]
  
  names(est$pov_status) <- names(se$pov_status) <- c("cnt_all", "m_below_pov", "m_above_pov",
    "f_below_pov", "f_above_pov", "pct_below_pov", "pct_below_pov_m", "pct_below_pov_f",
    paste(rep(c("m_below_pov", "f_below_pov", "m_above_pov", "f_above_pov"), each= 8),
          rep(c("u18", "18_24", "25_34", "35_44", "45_54", "55_64", "65_74", "75up"), 4), sep= "_"))
  
  ### pov_to_work
  est$pov_to_work$pct_below_pov_workFT   <- apply(est$pov_to_work[,c(3,6)],1,sum) / apply(est$pov_to_work[,c(3,6,10,13)],1,sum)
  est$pov_to_work$pct_below_pov_workFT_m <- est$pov_to_work[,3] / apply(est$pov_to_work[,c(3,10)],1,sum)
  est$pov_to_work$pct_below_pov_workFT_f <- est$pov_to_work[,6] / apply(est$pov_to_work[,c(6,13)],1,sum)
  est$pov_to_work$pct_below_pov_workPT   <- apply(est$pov_to_work[,c(4,7)],1,sum) / apply(est$pov_to_work[,c(4,7,11,14)],1,sum)
  est$pov_to_work$pct_below_pov_nowork   <- apply(est$pov_to_work[,c(5,8)],1,sum) / apply(est$pov_to_work[,c(5,8,12,15)],1,sum)
  
  se$pov_to_work$pct_below_pov_workFT   <- sqrt(apply(se$pov_to_work[,c(3,6)]^2,1,sum) - (
    est$pov_to_work$pct_below_pov_workFT^2 * apply(se$pov_to_work[,c(3,6,10,13)]^2,1,sum))) / 
    apply(est$pov_to_work[,c(3,6,10,13)],1,sum)
  se$pov_to_work$pct_below_pov_workFT_m <- sqrt(se$pov_to_work[,3]^2 - (
    est$pov_to_work$pct_below_pov_workFT_m^2 * apply(se$pov_to_work[,c(3,10)]^2,1,sum))) / 
    apply(est$pov_to_work[,c(3,10)],1,sum)
  se$pov_to_work$pct_below_pov_workFT_f <- sqrt(se$pov_to_work[,6]^2 - (
    est$pov_to_work$pct_below_pov_workFT_f^2 * apply(se$pov_to_work[,c(6,13)]^2,1,sum))) / 
    apply(est$pov_to_work[,c(6,13)],1,sum)
  se$pov_to_work$pct_below_pov_workPT   <- sqrt(apply(se$pov_to_work[,c(4,7)]^2,1,sum) - (
    est$pov_to_work$pct_below_pov_workPT^2 * apply(se$pov_to_work[,c(4,7,11,14)]^2,1,sum))) / 
    apply(est$pov_to_work[,c(4,7,11,14)],1,sum)
  se$pov_to_work$pct_below_pov_nowork   <- sqrt(apply(se$pov_to_work[,c(5,8)],1,sum) - (
    est$pov_to_work$pct_below_pov_nowork^2 * apply(se$pov_to_work[,c(5,8,12,15)]^2,1,sum))) / 
    apply(est$pov_to_work[,c(5,8,12,15)],1,sum)
  
  est$pov_to_work <- est$pov_to_work[, c(1,2,9,3:8,10:15,16:20)]
  se$pov_to_work  <- se$pov_to_work[, c(1,2,9,3:8,10:15,16:20)]
  
  names(est$pov_to_work) <- names(se$pov_to_work) <- c("cnt_all", "inc_below_pov", "inc_above_pov", 
    paste(rep(c("inc_below_pov", "inc_above_pov"), each= 6),
      rep(paste(rep(c("m","f"), each= 3),
      rep(c("workFT", "workPT", "nowork"), 2), sep= "_"), 2), sep= "_"), 
    "pct_below_pov_workFT", "pct_below_pov_workFT_m", "pct_below_pov_workFT_f", 
    "pct_below_pov_workPT", "pct_below_pov_nowork")
  
  ### dis_by_age
  est$dis_by_age <- est$dis_by_age[, -c(seq(3,18,3), seq(22,37,3))]
  est$dis_by_age <- est$dis_by_age[, c(1:2,15, 3:14, 16:27)]
  se$dis_by_age  <- se$dis_by_age[, -c(seq(3,18,3), seq(22,37,3))]
  se$dis_by_age <- se$dis_by_age[, c(1:2,15, 3:14, 16:27)]
  names(est$dis_by_age) <- names(se$dis_by_age) <- c("all", "male", "female", paste(
    rep(c("m", "f"), each= 12),
    paste(rep(c("u5", "5_17", "18_34", "35_64", "65_74", "75up"), each= 2), 
          rep(c("disability", "no_disability"), 6), sep= "_"), sep= "_"))
  
  est$dis_by_age$pct_disabled <- apply(est$dis_by_age[, seq(4,26,2)], 1, sum) / est$dis_by_age[,1]
  est$dis_by_age$pct_disabled_m <- apply(est$dis_by_age[, seq(4,14,2)], 1, sum) / est$dis_by_age[,2]
  est$dis_by_age$pct_disabled_f <- apply(est$dis_by_age[, seq(16,26,2)], 1, sum) / est$dis_by_age[,3]
  
  se$dis_by_age$pct_disabled <- sqrt(apply(se$dis_by_age[, seq(4,26,2)]^2, 1, sum) - (
    est$dis_by_age$pct_disabled^2 * se$dis_by_age[,1]^2)) / est$dis_by_age[,1]
  se$dis_by_age$pct_disabled_m <- sqrt(apply(se$dis_by_age[, seq(4,14,2)]^2, 1, sum) - (
    est$dis_by_age$pct_disabled_m^2 * se$dis_by_age[,2]^2)) / est$dis_by_age[,2]
  se$dis_by_age$pct_disabled_f <- sqrt(apply(se$dis_by_age[, seq(16,26,2)]^2, 1, sum) - (
    est$dis_by_age$pct_disabled_f^2 * se$dis_by_age[,3]^2)) / est$dis_by_age[,3]
  
  ### hh_inc
  names(est$hh_inc) <- names(se$hh_inc) <- c("cnt_all", paste(
    rep("hh_inc", 16),
    c("lt10k", "10_lt15k", "15_lt20k", "20_lt25k", "25_lt30k", "30_lt35k", "35_lt40k",
      "40_lt45k", "45_lt50k", "50_lt60k", "60_lt75k", "75_lt100k", "100_lt125k",
      "125_lt150k", "150_lt200k", "200k_up"), sep= "_"))
  
  ### med_hh_inc
  names(est$med_hh_inc) <- names(se$med_hh_inc) <- "med_hh_inc"
  
  ### hh_ss_inc
  names(est$hh_ss_inc) <- names(se$hh_ss_inc) <- c("ss_inc_total", "ss_inc_wSS", "ss_inc_woSS")
  
  ### hh_psa_inc
  names(est$hh_psa_inc) <- names(se$hh_psa_inc) <- c("psa_inc_total", "psa_inc_wPSA", "psa_inc_woPSA")
  
  # 03 -- sort, combine, and return
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
  names(ret$estimates) <- names(ret$standard_error) <- c("pov_status_by_age_sex", "pov_status_by_work_exp",
    "sex_age_disability_status", "hh_inc", "median_hh_inc", "soc_sec_inc_hh", "pub_assist_inc_hh")
  
  return(ret)
}
