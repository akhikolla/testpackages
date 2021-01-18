

#' @title Pull ACS race data
#' @description Pull ACS data for a specified geography from base tables
#' B01001B-I and B02001. ' These tables reference population counts by race.
#' @param endyear An integer, indicating the latest year of the data in the survey.
#' @param span An integer in \code{c(1,3,5)} indicating the span of the desired data.
#' @param geography a valid \code{geo.set} object specifying the census geography or 
#' geographies to be fetched.
#' @return A \code{list} containing the endyear, span, a \code{data.frame} of estimates,
#' a \code{data.frame} of standard errors, and a \code{data.frame} of the geography 
#' metadata from \code{\link[acs]{acs.fetch}}.
#' @seealso \code{\link[acs]{acs.fetch}}, \code{\link[acs]{geo.make}}
#' @export
pull_race_data <- function(endyear, span, geography) {
  # 00 -- error checking
  #----------------------------------------------
  check_geo_inputs(endyear= endyear, span= span, geography= geography)
  
  # 01 -- pull data
  #----------------------------------------------
  oldw <- getOption("warn")
  options(warn= -1) # suppress warnings from library(acs) / ACS API
  on.exit(options(warn= oldw)) # turn warnings back on
  race_all <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                        table.number = "B02001", col.names = "pretty")
  race_aa <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                       table.number = "B01001B", col.names = "pretty")
  race_nat <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                        table.number = "B01001C", col.names = "pretty")
  race_asian <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                       table.number = "B01001D", col.names = "pretty")
  race_isl <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                       table.number = "B01001E", col.names = "pretty")
  race_oth <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                       table.number = "B01001F", col.names = "pretty")
  race_2p <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                        table.number = "B01001G", col.names = "pretty")
  race_white <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                       table.number = "B01001H", col.names = "pretty")
  race_hisp <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                       table.number = "B01001I", col.names = "pretty")
  
  # 02 -- create lists of EST and SE -- as data.frames
  #----------------------------------------------
  est <- list(tot_pop= data.frame(race_all@estimate[, c(1:2), drop= FALSE]),
              aa_pop = data.frame(race_aa@estimate[,c(1,2,17), drop= FALSE]),
              nat_pop= data.frame(race_nat@estimate[,c(1,2,17), drop= FALSE]),
              asn_pop= data.frame(race_asian@estimate[,c(1,2,17), drop= FALSE]),
              isl_pop= data.frame(race_isl@estimate[,c(1,2,17), drop= FALSE]),
              oth_pop= data.frame(race_oth@estimate[,c(1,2,17), drop= FALSE]),
              r2_pop = data.frame(race_2p@estimate[,c(1,2,17), drop= FALSE]),
              whi_pop= data.frame(race_white@estimate[,c(1,2,17), drop= FALSE]),
              his_pop= data.frame(race_hisp@estimate[,c(1,2,17), drop= FALSE]))
  
  se <- list(tot_pop= data.frame(race_all@standard.error[, c(1:2), drop= FALSE]),
              aa_pop = data.frame(race_aa@standard.error[,c(1,2,17), drop= FALSE]),
              nat_pop= data.frame(race_nat@standard.error[,c(1,2,17), drop= FALSE]),
              asn_pop= data.frame(race_asian@standard.error[,c(1,2,17), drop= FALSE]),
              isl_pop= data.frame(race_isl@standard.error[,c(1,2,17), drop= FALSE]),
              oth_pop= data.frame(race_oth@standard.error[,c(1,2,17), drop= FALSE]),
              r2_pop = data.frame(race_2p@standard.error[,c(1,2,17), drop= FALSE]),
              whi_pop= data.frame(race_white@standard.error[,c(1,2,17), drop= FALSE]),
              his_pop= data.frame(race_hisp@standard.error[,c(1,2,17), drop= FALSE]))
  
  geo <- race_all@geography
  
  rm(race_all, race_aa, race_nat, race_asian, race_isl, race_oth, race_2p, race_white, race_hisp)
  
  # 03 -- combine columns
  #----------------------------------------------
  est <- data.frame(do.call("cbind", est)); se <- data.frame(do.call("cbind", se))
  est <- est[, c(1:2, seq(3,26,3), seq(4,26,3), seq(5,26,3))]
  se <-  se[, c(1:2, seq(3,26,3), seq(4,26,3), seq(5,26,3))]
  
  names(est) <- names(se) <- c(paste("agg", c("count", "white"), sep= "_"), 
    paste(rep(c("total", "m", "f"), each= 8),
      rep(c("black", "nat_amer", "asian", "pac_isl", "other", "2p_races", "white_alone", "hisp_lat"), 3), sep="_"))
  
  
  # 04 -- sort, combine and return
  #----------------------------------------------
  geo_sorted <- geo_alphabetize(geo= geo, est= est, se= se)
  geo <- geo_sorted[["geo"]]
  est <- geo_sorted[["est"]]
  se <- geo_sorted[["se"]]
  
  ret <- list(endyear= endyear, span= span,
              estimates= list(pop_by_race= est),
              standard_error= list(pop_by_race= se),
              geography= geo,
              geo_title= unlist(geography@geo.list))
  class(ret) <- "macroACS"
  return(ret)
  
}

