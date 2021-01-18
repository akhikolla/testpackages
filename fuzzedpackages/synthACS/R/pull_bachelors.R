
#' @title Pull ACS data on field of bachelor's degree
#' @description Pull ACS data for a specified geography from base tables
#' B15011 and B15012. Note: only 2014 data is supplied by ACS
#' @param endyear An integer, indicating the latest year of the data in the survey.
#' @param span An integer in \code{c(1,3,5)} indicating the span of the desired data.
#' @param geography a valid \code{geo.set} object specifying the census geography or 
#' geographies to be fetched.
#' @return A \code{list} containing the endyear, span, a \code{data.frame} of estimates,
#' a \code{data.frame} of standard errors, and a \code{data.frame} of the geography 
#' metadata from \code{\link[acs]{acs.fetch}}.
#' @seealso \code{\link[acs]{acs.fetch}}, \code{\link[acs]{geo.make}}
#' @export
pull_bachelors <- function(endyear, span, geography) {
  # 00 -- error checking
  #----------------------------------------------
  check_geo_inputs(endyear= endyear, span= span, geography= geography)
  
  # 01 -- pull data and move to lists
  #----------------------------------------------
  oldw <- getOption("warn")
  options(warn= -1) # suppress warnings from library(acs) / ACS API
  on.exit(options(warn= oldw)) # turn warnings back on
  # ba_deg25up <- acs.fetch(endyear= 2014, span= 5, geography= la_tracts, table.number = "C15010",
  #                    col.names= "pretty")
  by_sex_age <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, table.number = "B15011",
                     col.names = "pretty")
  ba_total <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, table.number = "B15012",
                     col.names = "pretty")
  
  est <- list(by_sex_age= data.frame(by_sex_age@estimate[,c(1,2,21, 3:20, 22:39), drop= FALSE]),
              ba_total= data.frame(ba_total@estimate))
  se <- list(by_sex_age= data.frame(by_sex_age@standard.error[,c(1,2,21, 3:20, 22:39), drop= FALSE]),
             ba_total= data.frame(ba_total@standard.error))
  geo <- by_sex_age@geography
  
  rm(by_sex_age, ba_total)
  
  # 02 -- mung data
  #----------------------------------------------
  est$by_sex_age <- est$by_sex_age[,-c(seq(4,34,6))]
  se$by_sex_age  <- se$by_sex_age[,-c(seq(4,34,6))]
  names(est$by_sex_age) <- names(se$by_sex_age) <- c("ba_tot", "cnt_male", "cnt_female", paste(
    rep(c("m", "f"), each= 15), rep(paste(
      rep(c("25_39", "40_64", "65up"), each= 5),
      rep(c("STEM_SocSci", "sci_eng_related", "business", "education", "arts_humanities_oth"), 3), 
      sep= "_"), 2), sep= "_"))
  
  names(est$ba_total) <- names(se$ba_total) <- c("ba_tot", "comp_math_stat", "bio_ag_env_sci",
      "phys_sci", "psych", "soc_sci", "engineering", "mult_disp_studies", "sci_eng_related", "business", 
      "education", "lit_lang", "lib_arts_history", "arts", "communications", "other")
  
  
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
  names(ret$estimates) <- names(ret$standard_error) <- c("degrees_by_age_sex", "degrees_by_field")
  return(ret)
  
}