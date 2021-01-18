
#' @title Pull ACS marital status data
#' @description Pull ACS data for a specified geography from base tables
#' B12001, B12006, B12007, 12501
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
pull_mar_status <- function(endyear, span, geography) {
  # 00 -- error checking
  #----------------------------------------------
  check_geo_inputs(endyear= endyear, span= span, geography= geography)
  
  # 01 -- pull data and move to lists
  #----------------------------------------------
  oldw <- getOption("warn")
  options(warn= -1) # suppress warnings from library(acs) / ACS API
  on.exit(options(warn= oldw)) # turn warnings back on
  by_sex <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                        table.number = "B12001", col.names= "pretty")
  by_labor <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                        table.number = "B12006", col.names= "pretty")
  med_age_mar <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                        table.number = "B12007", col.names= "pretty")
  num_mar <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                        table.number = "B12501", col.names= "pretty")
  
  est <- list(by_sex= data.frame(by_sex@estimate),
              by_labor= data.frame(by_labor@estimate),
              med_age_mar= data.frame(med_age_mar@estimate),
              num_mar= data.frame(num_mar@estimate))
  
  se <- list(by_sex= data.frame(by_sex@standard.error),
              by_labor= data.frame(by_labor@standard.error),
              med_age_mar= data.frame(med_age_mar@standard.error),
              num_mar= data.frame(num_mar@standard.error))
  
  geo <- by_sex@geography
  
  rm(by_sex, by_labor, med_age_mar, num_mar)
  
  ## 02 -- (A) combine columns and (B) calc percentages
  ### by_sex
  #----------------------------------------------
  est$by_sex <- est$by_sex[, c(1,3,5:6,9:10,12,14:15,18:19)]
  se$by_sex  <- se$by_sex[, c(1,3,5:6,9:10,12,14:15,18:19)]
  names(est$by_sex) <- names(se$by_sex) <- c("all", "nvr_mar_m", "now_mar_spouse_pres_m",
    "now_mar_spouse_abs_m", "widow_m", "divorced_m", "nvr_mar_f", "now_mar_spouse_pres_f",
    "now_mar_spouse_abs_f", "widow_f", "divorced_f")
  
  ### by_labor
  est$by_labor <- est$by_labor[, c(1,2,13,24,35,46,
                                   5:7,10:12,16:18,21:23,27:29,32:34,37:39,43:45,48:50,54:56)]
  se$by_labor  <- se$by_labor[, c(1,2,13,24,35,46,
                                   5:7,10:12,16:18,21:23,27:29,32:34,37:39,43:45,48:50,54:56)]
  names(est$by_labor) <- names(se$by_labor) <- c("all","nvr_mar", "now_mar", "sep", "widow", "div",
    paste(rep(c("nvr_mar", "now_mar", "sep", "widow", "div"), each= 6),
      rep(paste(rep(c("m", "f"), each= 3), 
            rep(c("employed", "unemployed", "not_labor_force"), 2), sep= "_"), 5), sep= "_"))
  
  est$by_labor$pct_emp <- apply(est$by_labor[, seq(7,34,3)], 1, sum) / 
    apply(est$by_labor[, c(seq(7,34,3), seq(8,35,3))], 1, sum)
  est$by_labor$pct_unemp <- apply(est$by_labor[, seq(8,35,3)], 1, sum) / 
    apply(est$by_labor[, c(seq(7,34,3), seq(8,35,3))], 1, sum)
  est$by_labor$pct_labor_particip <- apply(est$by_labor[,c(seq(7,34,3), seq(8,35,3))], 1, sum) / est$by_labor[,1]
  est$by_labor$pct_nvr_mar <- est$by_labor[,2] / est$by_labor[,1]
  est$by_labor$pct_now_mar <- est$by_labor[,3] / est$by_labor[,1]
  est$by_labor$pct_sep     <- est$by_labor[,4] / est$by_labor[,1]
  est$by_labor$pct_widow   <- est$by_labor[,5] / est$by_labor[,1]
  est$by_labor$pct_div     <- est$by_labor[,6] / est$by_labor[,1]
  
  se$by_labor$pct_emp <- apply(se$by_labor[, seq(7,34,3)], 1, sum) / 
    apply(se$by_labor[, c(seq(7,34,3), seq(8,35,3))], 1, sum)
  se$by_labor$pct_unemp <- apply(se$by_labor[, seq(8,35,3)], 1, sum) / 
    apply(se$by_labor[, c(seq(7,34,3), seq(8,35,3))], 1, sum)
  se$by_labor$pct_labor_particip <- apply(se$by_labor[,c(seq(7,34,3), seq(8,35,3))], 1, sum) / se$by_labor[,1]
  se$by_labor$pct_nvr_mar <- sqrt(se$by_labor[,2]^2 - (
    est$by_labor$pct_nvr_mar^2 * se$by_labor[,1]^2)) / est$by_labor[,1]
  se$by_labor$pct_now_mar <- sqrt(se$by_labor[,3]^2 - (
    est$by_labor$pct_now_mar^2 * se$by_labor[,1]^2)) / est$by_labor[,1]
  se$by_labor$pct_sep     <- sqrt(se$by_labor[,4]^2 - (
    est$by_labor$pct_sep^2 * se$by_labor[,1]^2)) / est$by_labor[,1]
  se$by_labor$pct_widow   <- sqrt(se$by_labor[,5]^2 - (
    est$by_labor$pct_widow^2 * se$by_labor[,1]^2)) / est$by_labor[,1]
  se$by_labor$pct_div     <- sqrt(se$by_labor[,6]^2 - (
    est$by_labor$pct_div^2 * se$by_labor[,1]^2)) / est$by_labor[,1]
  
  ### med_age_mar
  names(est$med_age_mar) <- names(se$med_age_mar) <- c("med_age_1st_mar_m", "med_age_1st_mar_f")
  
  ### num_mar
  names(est$num_mar) <- names(se$num_mar) <- c("all", paste(
    rep(c("m", "f"), each= 5),
    rep(c("all", "never_mar", "ever_mar", "ever_mar_last_yr", "ever_mar_not_last_yr"), 2), sep= "_"))
  
  
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
  names(ret$estimates) <- names(ret$standard_error) <- c("mar_status_by_sex", 
    "mar_status_by_labor_participation", "median_age_first_marriage", "marriages_last_year_by_sex_by_status")
  
  return(ret)
}
