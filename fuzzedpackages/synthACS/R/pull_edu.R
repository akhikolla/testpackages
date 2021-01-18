
#' @title Pull ACS educational attainment and enrollment data
#' @description Pull ACS data for a specified geography from base tables
#' B14001, B14003, B15001, B15002. Not currently implemented: B15010, B28006
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
pull_edu <- function(endyear, span, geography) {
  # 00 -- error checking
  #----------------------------------------------
  check_geo_inputs(endyear= endyear, span= span, geography= geography)
  
  # 01 -- pull data and move to lists
  #----------------------------------------------
  oldw <- getOption("warn")
  options(warn= -1) # suppress warnings from library(acs) / ACS API
  on.exit(options(warn= oldw)) # turn warnings back on
  edu_enroll <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                              table.number = "B14001", col.names= "pretty")
  enroll_details <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                          table.number = "B14003", col.names= "pretty")
  edu_attain18 <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                          table.number = "B15001", col.names= "pretty")
  edu_attain25 <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                          table.number = "B15002", col.names= "pretty")
#   deg_major25 <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
#                           table.number = "B15010", col.names= "pretty"))
#   edu_internet <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
#                             table.number = "B28006", col.names= "pretty")
  
  est <- list(edu_enroll= data.frame(edu_enroll@estimate[,-2, drop= FALSE]),
              enroll_details= data.frame(enroll_details@estimate[, -c(2,30), drop= FALSE]),
              edu_attain18= data.frame(edu_attain18@estimate[,-c(2,43), drop= FALSE]),
              edu_attain25= data.frame(edu_attain25@estimate[, -c(2,19), drop= FALSE])) #,
              # deg_major25= data.frame(deg_major25@estimate)))
  
  se <- list(edu_enroll= data.frame(edu_enroll@standard.error[,-2, drop= FALSE]),
              enroll_details= data.frame(enroll_details@standard.error[, -c(2,30), drop= FALSE]),
              edu_attain18= data.frame(edu_attain18@standard.error[,-c(2,43), drop= FALSE]),
              edu_attain25= data.frame(edu_attain25@standard.error[, -c(2,19), drop= FALSE]))
  
  geo <- edu_enroll@geography
  
  rm(edu_enroll, enroll_details, edu_attain18, edu_attain25)
  
  ## 02 -- (A) combine columns and (B) calc percentages
  ### edu_enroll
  #----------------------------------------------
  names(est$edu_enroll) <- names(se$edu_enroll) <- paste("enroll_cnt",
    c("total", "preschool", "kindergarten", "grade1_4", "grade5_8", 
      "grade8_12", "undergrad", "gradschool", "not_enrolled"), sep= "_")
  
  ### enroll_details
  names(est$enroll_details) <- names(se$enroll_details) <- c("all", paste(
          rep(c("m_pubSch", "m_priSch", "m_notSch", "f_pubSch", "f_priSch", "f_notSch"), each= 9),
          rep(c("all", "y3_4", "y5_9", "y10_14", "y15_17", "y18_19", "y20_24", "y25_34", 
                "y35up"), 6), sep= "_"))
  
  ### edu_attain18
  names(est$edu_attain18) <- names(se$edu_attain18) <- c("all", paste(
    rep(c("m18_24", "m25_34", "m35_44", "m45_64", "m65up",
          "f18_24", "f25_34", "f35_44", "f45_64", "f65up"), each= 8),
    rep(c("cnt", "lt_9grade", "lt_hs", "hs", "lt_col", "ass_deg", "ba_deg", "grad_deg"), 10), 
    sep= "_"))
  
  ### edu_attain25
  est$edu_attain25$new1 <- apply(est$edu_attain25[, 3:5], 1, sum) # m lt hs
  est$edu_attain25$new2 <- apply(est$edu_attain25[, 6:9], 1, sum) # m some HS
  est$edu_attain25$new3 <- apply(est$edu_attain25[, 11:12], 1, sum) # m some col
  est$edu_attain25$new4 <- apply(est$edu_attain25[, 19:21], 1, sum) # f lt hs
  est$edu_attain25$new5 <- apply(est$edu_attain25[, 22:25], 1, sum) # f some HS
  est$edu_attain25$new6 <- apply(est$edu_attain25[, 27:28], 1, sum) # f some col
  # m pct lt col
  est$edu_attain25$new7 <- apply(est$edu_attain25[, 3:13], 1, sum) / apply(est$edu_attain25[, 2:17], 1, sum)
  # m pct col
  est$edu_attain25$new8 <- est$edu_attain25[,14] / apply(est$edu_attain25[, 2:17], 1, sum) 
  # m pct grad
  est$edu_attain25$new9 <- apply(est$edu_attain25[, 15:17], 1, sum) / apply(est$edu_attain25[, 2:17], 1, sum) 
  # f pct lt col
  est$edu_attain25$new10 <- apply(est$edu_attain25[, 19:29], 1, sum) / apply(est$edu_attain25[, 18:33], 1, sum)
  # f pct col
  est$edu_attain25$new11 <- est$edu_attain25[,30] / apply(est$edu_attain25[, 18:33], 1, sum)
  # f pct grad
  est$edu_attain25$new12 <- apply(est$edu_attain25[, 31:33], 1, sum) / apply(est$edu_attain25[, 18:33], 1, sum)
  # all pct no school
  est$edu_attain25$new13 <- apply(est$edu_attain25[, c(2,18)], 1, sum) / est$edu_attain25[, 1]
  # all pct lt col
  est$edu_attain25$new14 <- apply(est$edu_attain25[, c(3:13, 19:29)], 1, sum) / est$edu_attain25[, 1]
  # all pct col
  est$edu_attain25$new15 <- apply(est$edu_attain25[, c(14,30)], 1, sum) / est$edu_attain25[, 1]
  # all pct grad
  est$edu_attain25$new16 <- apply(est$edu_attain25[, c(15:17, 31:33)], 1, sum) / est$edu_attain25[, 1]
  
  se$edu_attain25$new1 <- sqrt(apply(se$edu_attain25[, 3:5]^2, 1, sum)) # m lt hs
  se$edu_attain25$new2 <- sqrt(apply(se$edu_attain25[, 6:9]^2, 1, sum)) # m some HS
  se$edu_attain25$new3 <- sqrt(apply(se$edu_attain25[, 11:12]^2, 1, sum)) # m some col
  se$edu_attain25$new4 <- sqrt(apply(se$edu_attain25[, 19:21]^2, 1, sum)) # f lt hs
  se$edu_attain25$new5 <- sqrt(apply(se$edu_attain25[, 22:25]^2, 1, sum)) # f some HS
  se$edu_attain25$new6 <- sqrt(apply(se$edu_attain25[, 27:28]^2, 1, sum)) # f some col
  
  se$edu_attain25$new7 <- sqrt(apply(se$edu_attain25[, 3:13]^2, 1, sum) - (
    est$edu_attain25$new7^2 * apply(se$edu_attain25[, 2:17]^2, 1, sum)
    )) / apply(est$edu_attain25[, 2:17], 1, sum)
  se$edu_attain25$new8 <- sqrt(se$edu_attain25[,14]^2 - (est$edu_attain25$new8^2 * 
     apply(se$edu_attain25[, 2:17]^2, 1, sum))) / apply(est$edu_attain25[, 2:17], 1, sum)
  se$edu_attain25$new9 <- sqrt(apply(se$edu_attain25[, 15:17]^2, 1, sum) - (
    est$edu_attain25$new9^2 * apply(se$edu_attain25[, 2:17]^2, 1, sum)
  )) / apply(est$edu_attain25[, 2:17], 1, sum)
  se$edu_attain25$new10 <- sqrt(apply(se$edu_attain25[, 19:29]^2, 1, sum) - (
    est$edu_attain25$new10^2 * apply(se$edu_attain25[, 18:33]^2, 1, sum)
  )) / apply(est$edu_attain25[, 18:33], 1, sum)
  se$edu_attain25$new11 <- sqrt(se$edu_attain25[,30]^2 - (est$edu_attain25$new11^2 * 
     apply(se$edu_attain25[, 18:33]^2, 1, sum))) / apply(est$edu_attain25[, 18:33], 1, sum)
  se$edu_attain25$new12 <- sqrt(apply(se$edu_attain25[, 31:33]^2, 1, sum) - (
    est$edu_attain25$new12^2 * apply(se$edu_attain25[, 18:33]^2, 1, sum)
  )) / apply(est$edu_attain25[, 18:33], 1, sum)
  se$edu_attain25$new13 <- sqrt(apply(se$edu_attain25[, c(2,18)]^2, 1, sum) - (
    est$edu_attain25$new13^2 * se$edu_attain25[, 1]^2
  )) / est$edu_attain25[,1]
  se$edu_attain25$new14 <- sqrt(apply(se$edu_attain25[, c(3:13, 19:29)]^2, 1, sum) - (
    est$edu_attain25$new14^2 * se$edu_attain25[, 1]
  )) / est$edu_attain25[,1]
  se$edu_attain25$new15 <- sqrt(apply(se$edu_attain25[, c(14,30)]^2, 1, sum) - (
    est$edu_attain25$new15^2 * se$edu_attain25[, 1]
  )) / est$edu_attain25[,1]
  se$edu_attain25$new16 <- sqrt(apply(se$edu_attain25[, c(15:17, 31:33)]^2, 1, sum) - (
    est$edu_attain25$new16^2 * se$edu_attain25[, 1]
  )) / est$edu_attain25[,1]
  
  est$edu_attain25 <- est$edu_attain25[, c(1:2,34,35,10,36,13:17,40:42, 18,37:38,26,39,29:33,43:49)]
  se$edu_attain25 <- se$edu_attain25[, c(1:2,34,35,10,36,13:17,40:42, 18,37:38,26,39,29:33,43:49)]
  
  names(est$edu_attain25) <- names(se$edu_attain25) <- c("all", paste(
    rep(c("m25up", "f25up"), each= 13),
    rep(c("no_sch", "ltHS", "someHS", "hs_deg", "some_col", "ass_deg", "ba_deg", "ma_deg",
          "prof_deg", "PhD", "pct_ltcol", "pct_col", "pct_grad"), 2), sep= "_"),
    "pct_no_school", "pct_ltcol", "pct_col", "pct_grad")
  
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
  names(ret$estimates) <- names(ret$standard_error) <- c("school_enrollment_by_type", 
      "school_enrollment_by_age", "edu_attain_by_age_sex18", "edu_attain_age25up")
  return(ret)
  
}