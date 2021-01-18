#' @title Pull ACS base tables
#' @description A wrapper function to pull multiple base tables from ACS API via 
#' \code{\link[acs]{acs.fetch}}. 
#' @param endyear An integer, indicating the latest year of the data in the survey.
#' @param span An integer in \code{c(1,3,5)} indicating the span of the desired data.
#' @param geography a valid \code{geo.set} object specifying the census geography or 
#' geographies to be fetched.
#' @param table_vec A \code{character} vector specifying ACS base tables.
#' @return A \code{'macroACS'} class object
#' @export
#' @references \url{https://www.census.gov/programs-surveys/acs/technical-documentation/summary-file-documentation.html}
#' @examples \dontrun{
#' # make geography
#' la_geo <- acs::geo.make(state= "CA", county= "Los Angeles")
#' # pull data 
#' la_dat <- pull_acs_basetables(endyear= 2015, span= 1, geography= la_geo, 
#'   table_vec= c("B01001", "B01002", "B01003"))
#' }
pull_acs_basetables <- function(endyear, span, geography, table_vec) {
  check_geo_inputs(endyear, span, geography)
  if (!is.character(table_vec) | length(table_vec) < 1L) 
    stop("table_vec must be a character vector of ACS base table numbers")
  
  nr <- length(table_vec) 
  temp_dat <- vector("list", length= nr)
  oldw <- getOption("warn")
  options(warn= -1) # suppress warnings from library(acs) / ACS API
  on.exit(options(warn= oldw)) # turn warnings back on
  for (i in 1:nr) {
    temp_dat[[i]] <- acs::acs.fetch(endyear, span, geography, table.number= table_vec[i],
                                    col.names= "pretty")
  }
  
  ret <- list(endyear= endyear, span= span,
              estimates=lapply(temp_dat, function(l) {return(as.data.frame(l@estimate))}),
              standard_error= lapply(temp_dat, function(l) {return(as.data.frame(l@standard.error))}),
              geography= temp_dat[[1]]@geography,
              geo_title= unlist(geography@geo.list)
  )
  class(ret) <- "macroACS"
  return(ret)
}

#' @title Pull ACS data for synthetic data creation.
#' @description Pull ACS data for a specified geography from base tables
#' B01001, B02001, B12002, B15001, B06001, B06010, B23001, B17005, and B17005.
#' These tables reference population counts by a number of slices.
#' Multiple additional fields, mainly percentages and aggregations, are calculated.
#' @param endyear An integer, indicating the latest year of the data in the survey.
#' @param span An integer in \code{c(1,3,5)} indicating the span of the desired data.
#' @param geography a valid \code{geo.set} object specifying the census geography or 
#' geographies to be fetched.
#' @return A \code{list} containing the endyear, span, a list of \code{data.frame}s of estimates,
#' a list of \code{data.frame}s of standard errors, 
#' and the geography metadata from \code{\link[acs]{acs.fetch}}.
#' @seealso \code{\link[acs]{acs.fetch}}, \code{\link[acs]{geo.make}}
#' @export
#' 
#' @examples \dontrun{
#' # make geography
#' la_geo <- acs::geo.make(state= "CA", county= "Los Angeles", tract= "*")
#' # pull data elements for creating synthetic data
#' la_dat <- pull_synth_data(2014, 5, la_geo)
#' }
pull_synth_data <- function(endyear, span, geography) {
  # 00 -- error checking
  #----------------------------------------------
  check_geo_inputs(endyear= endyear, span= span, geography= geography)
  
  # 01 -- pull data
  #----------------------------------------------
  oldw <- getOption("warn")
  options(warn= -1) # suppress warnings from library(acs) / ACS API
  on.exit(options(warn= oldw)) # turn warnings back on
  age_by_sex <- acs::acs.fetch(endyear= endyear, span= span, geography= geography, 
                          table.number = "B01001", col.names= "pretty") # total pop
  pop_by_race <- pull_race_data(endyear= endyear, span= span, geography= geography)
  marital_status <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                                 table.number = "B12002", col.names = "pretty") # (age 15+, by age and gender)
  edu <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                      table.number = "B15001", col.names = "pretty") # (age 18+, by age and gender)
  nativity <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                        table.number = "B06001", col.names = "pretty") # by age; total US population
  by_inc_12mo <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                           table.number = "B06010", col.names = "pretty") #  15+ years by nativity status
  # geo_mob_mar_stat <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
  #                          table.number = "B07008", col.names = "pretty") # 15+
  geo_mob_inc <- acs::acs.fetch(endyear = endyear, span= span, geography = geography,
                           table.number = "B07010", col.names = "pretty") # 15+
  geo_mob_edu <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                           table.number = "B07009", col.names = "pretty") # 25+
  emp_status <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                          table.number = "B23001", col.names = "pretty") # 16+ by age/sex
  pov_status1 <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
                             table.number = "B17005", col.names = "pretty") # (age 16+, by sex & employment status)
  # pov_status2 <- acs::acs.fetch(endyear = endyear, span= span, geography = geography, 
  #                          table.number = "B17007", col.names = "pretty") # B17007 - age 15+, by age and gender ... 
  
  # 02 -- create lists of EST and SE -- as data.frames
  #----------------------------------------------
  est <- list(age_by_sex= data.frame(age_by_sex@estimate),
              pop_by_race= pop_by_race$estimate$pop_by_race,
              marital_status= data.frame(marital_status@estimate),
              edu= data.frame(edu@estimate),
              nativity= data.frame(nativity@estimate),
              by_inc_12mo= data.frame(by_inc_12mo@estimate),
              # geo_mob_mar_stat= data.frame(geo_mob_mar_stat@estimate[, c(1,7:dim(geo_mob_mar_stat@estimate)[[2]])]),
              geo_mob_edu= data.frame(
                geo_mob_edu@estimate[, c(1,7:dim(geo_mob_edu@estimate)[[2]]), drop= FALSE]
              ),
              ind_inc= data.frame(geo_mob_inc@estimate[, 1:11, drop= FALSE]),
              emp_status= data.frame(emp_status@estimate),
              pov_status1= data.frame(pov_status1@estimate)) #,
              # pov_status2= data.frame(pov_status2@estimate))
  
  se <- list(age_by_sex= data.frame(age_by_sex@standard.error),
             pop_by_race= pop_by_race$standard_error$pop_by_race,
             marital_status= data.frame(marital_status@standard.error),
             edu= data.frame(edu@standard.error),
             nativity= data.frame(nativity@standard.error),
             by_inc_12mo= data.frame(by_inc_12mo@standard.error),
             # geo_mob_mar_stat= data.frame(geo_mob_mar_stat@standard.error[, c(1,7:dim(geo_mob_mar_stat@standard.error)[[2]])]),
             geo_mob_edu= data.frame(
               geo_mob_edu@standard.error[, c(1,7:dim(geo_mob_edu@standard.error)[[2]]), drop= FALSE]
             ),
             ind_inc= data.frame(geo_mob_inc@standard.error[,1:11, drop= FALSE]),
             emp_status= data.frame(emp_status@standard.error),
             pov_status1= data.frame(pov_status1@standard.error)) #,
             # pov_status2= data.frame(pov_status2@standard.error))
  
  geo <- age_by_sex@geography
  
  # 03 -- combine columns
  #----------------------------------------------
  ### age by sex
  #-----------
  est$age_by_sex$m_u15 <- apply(est$age_by_sex[, c(3:5)],1, sum)
  est$age_by_sex$m_18_24 <- apply(est$age_by_sex[, c(7:10)],1, sum)
  est$age_by_sex$f_u15 <- apply(est$age_by_sex[, c(27:29)],1, sum)
  est$age_by_sex$f_18_24 <- apply(est$age_by_sex[, c(31:34)],1, sum)
  est$age_by_sex$m_60_64 <- apply(est$age_by_sex[, c(18:19)],1, sum)
  est$age_by_sex$m_65_69 <- apply(est$age_by_sex[, c(20:21)],1, sum)
  est$age_by_sex$f_60_64 <- apply(est$age_by_sex[, c(42:43)],1, sum)
  est$age_by_sex$f_65_69 <- apply(est$age_by_sex[, c(44:45)],1, sum)
  
  se$age_by_sex$m_u15 <- sqrt(apply(se$age_by_sex[, c(3:5)]^2, 1, sum))
  se$age_by_sex$m_18_24 <- sqrt(apply(se$age_by_sex[, c(7:10)]^2, 1, sum))
  se$age_by_sex$f_u15 <- sqrt(apply(se$age_by_sex[, c(27:29)]^2, 1, sum))
  se$age_by_sex$f_18_24 <- sqrt(apply(se$age_by_sex[, c(31:34)]^2, 1, sum))
  se$age_by_sex$m_60_64 <- sqrt(apply(se$age_by_sex[, c(18:19)],1, sum))
  se$age_by_sex$m_65_69 <- sqrt(apply(se$age_by_sex[, c(20:21)],1, sum))
  se$age_by_sex$f_60_64 <- sqrt(apply(se$age_by_sex[, c(42:43)],1, sum))
  se$age_by_sex$f_65_69 <- sqrt(apply(se$age_by_sex[, c(44:45)],1, sum))
  
  est$age_by_sex <- est$age_by_sex[, c(1,2,26,50,6,51,11:17,54:55,22:25,52,30,53,35:41,56:57,46:49)]
  se$age_by_sex  <- se$age_by_sex[, c(1,2,26,50,6,51,11:17,54:55,22:25,52,30,53,35:41,56:57,46:49)]
  
  names(est$age_by_sex) <- names(se$age_by_sex) <- c("pop_age_cnt_all", "male_pop", "female_pop", 
    paste(rep(c("m", "f"), each= 16),
          rep(c("u15", "15_17", "18_24", "25_29", "30_34","35_39", "40_44", "45_49", "50_54", 
                "55_59", "60_64", "65_69","70_74", "75_79", "80_84", "85up"), 2), sep= "_"))
  
  ### marital_status
  #-----------
  # est
  est$marital_status$m_18_24nvr_married <- apply(est$marital_status[,5:6], 1, sum)
  est$marital_status$m_18_24married <- apply(est$marital_status[,21:22], 1, sum)
  est$marital_status$m_18_24widowed <- apply(est$marital_status[,67:68], 1, sum)
  est$marital_status$m_18_24divorced <- apply(est$marital_status[,82:83], 1, sum)
  est$marital_status$m_15_17mar_apart <- apply(est$marital_status[,c(36,51)], 1, sum)
  est$marital_status$m_18_24mar_apart <- apply(est$marital_status[,c(37,38,52,53)], 1, sum)
  est$marital_status$m_25_29mar_apart <- apply(est$marital_status[,c(39,54)], 1, sum)
  est$marital_status$m_30_34mar_apart <- apply(est$marital_status[,c(40,55)], 1, sum)
  est$marital_status$m_35_39mar_apart <- apply(est$marital_status[,c(41,56)], 1, sum)
  est$marital_status$m_40_44mar_apart <- apply(est$marital_status[,c(42,57)], 1, sum)
  est$marital_status$m_45_49mar_apart <- apply(est$marital_status[,c(43,58)], 1, sum)
  est$marital_status$m_50_54mar_apart <- apply(est$marital_status[,c(44,59)], 1, sum)
  est$marital_status$m_55_59mar_apart <- apply(est$marital_status[,c(45,60)], 1, sum)
  est$marital_status$m_60_64mar_apart <- apply(est$marital_status[,c(46,61)], 1, sum)
  est$marital_status$m_65_74mar_apart <- apply(est$marital_status[,c(47,62)], 1, sum)
  est$marital_status$m_75_84mar_apart <- apply(est$marital_status[,c(48,63)], 1, sum)
  est$marital_status$m_85up_mar_apart <- apply(est$marital_status[,c(49,64)], 1, sum)
  
  est$marital_status$f_18_24nvr_married <- apply(est$marital_status[,98:99], 1, sum)
  est$marital_status$f_18_24married <- apply(est$marital_status[,114:115], 1, sum)
  est$marital_status$f_18_24widowed <- apply(est$marital_status[,67:68], 1, sum)
  est$marital_status$f_18_24divorced <- apply(est$marital_status[,82:83], 1, sum)
  est$marital_status$f_15_17mar_apart <- apply(est$marital_status[,c(129,144)], 1, sum)
  est$marital_status$f_18_24mar_apart <- apply(est$marital_status[,c(130:131,145:146)], 1, sum)
  est$marital_status$f_25_29mar_apart <- apply(est$marital_status[,c(132,147)], 1, sum)
  est$marital_status$f_30_34mar_apart <- apply(est$marital_status[,c(133,148)], 1, sum)
  est$marital_status$f_35_39mar_apart <- apply(est$marital_status[,c(134,149)], 1, sum)
  est$marital_status$f_40_44mar_apart <- apply(est$marital_status[,c(135,150)], 1, sum)
  est$marital_status$f_45_49mar_apart <- apply(est$marital_status[,c(136,151)], 1, sum)
  est$marital_status$f_50_54mar_apart <- apply(est$marital_status[,c(137,152)], 1, sum)
  est$marital_status$f_55_59mar_apart <- apply(est$marital_status[,c(138,153)], 1, sum)
  est$marital_status$f_60_64mar_apart <- apply(est$marital_status[,c(139,154)], 1, sum)
  est$marital_status$f_65_74mar_apart <- apply(est$marital_status[,c(140,155)], 1, sum)
  est$marital_status$f_75_84mar_apart <- apply(est$marital_status[,c(141,156)], 1, sum)
  est$marital_status$f_85up_mar_apart <- apply(est$marital_status[,c(142,157)], 1, sum)
  
  # se
  se$marital_status$m_18_24nvr_married <- sqrt(apply(se$marital_status[,5:6], 1, sum))
  se$marital_status$m_18_24married <- sqrt(apply(se$marital_status[,21:22], 1, sum))
  se$marital_status$m_18_24widowed <- sqrt(apply(se$marital_status[,67:68], 1, sum))
  se$marital_status$m_18_24divorced <- sqrt(apply(se$marital_status[,82:83], 1, sum))
  se$marital_status$m_15_17mar_apart <- sqrt(apply(se$marital_status[,c(36,51)], 1, sum))
  se$marital_status$m_18_24mar_apart <- sqrt(apply(se$marital_status[,c(37,38,52,53)], 1, sum))
  se$marital_status$m_25_29mar_apart <- sqrt(apply(se$marital_status[,c(39,54)], 1, sum))
  se$marital_status$m_30_34mar_apart <- sqrt(apply(se$marital_status[,c(40,55)], 1, sum))
  se$marital_status$m_35_39mar_apart <- sqrt(apply(se$marital_status[,c(41,56)], 1, sum))
  se$marital_status$m_40_44mar_apart <- sqrt(apply(se$marital_status[,c(42,57)], 1, sum))
  se$marital_status$m_45_49mar_apart <- sqrt(apply(se$marital_status[,c(43,58)], 1, sum))
  se$marital_status$m_50_54mar_apart <- sqrt(apply(se$marital_status[,c(44,59)], 1, sum))
  se$marital_status$m_55_59mar_apart <- sqrt(apply(se$marital_status[,c(45,60)], 1, sum))
  se$marital_status$m_60_64mar_apart <- sqrt(apply(se$marital_status[,c(46,61)], 1, sum))
  se$marital_status$m_65_74mar_apart <- sqrt(apply(se$marital_status[,c(47,62)], 1, sum))
  se$marital_status$m_75_84mar_apart <- sqrt(apply(se$marital_status[,c(48,63)], 1, sum))
  se$marital_status$m_85up_mar_apart <- sqrt(apply(se$marital_status[,c(49,64)], 1, sum))
  
  se$marital_status$f_18_24nvr_married <- sqrt(apply(se$marital_status[,98:99], 1, sum))
  se$marital_status$f_18_24married <- sqrt(apply(se$marital_status[,114:115], 1, sum))
  se$marital_status$f_18_24widowed <- sqrt(apply(se$marital_status[,160:161], 1, sum))
  se$marital_status$f_18_24divorced <- sqrt(apply(se$marital_status[,175:176], 1, sum))
  se$marital_status$f_15_17mar_apart <- sqrt(apply(se$marital_status[,c(129,144)], 1, sum))
  se$marital_status$f_18_24mar_apart <- sqrt(apply(se$marital_status[,c(130:131,145:146)], 1, sum))
  se$marital_status$f_25_29mar_apart <- sqrt(apply(se$marital_status[,c(132,147)], 1, sum))
  se$marital_status$f_30_34mar_apart <- sqrt(apply(se$marital_status[,c(133,148)], 1, sum))
  se$marital_status$f_35_39mar_apart <- sqrt(apply(se$marital_status[,c(134,149)], 1, sum))
  se$marital_status$f_40_44mar_apart <- sqrt(apply(se$marital_status[,c(135,150)], 1, sum))
  se$marital_status$f_45_49mar_apart <- sqrt(apply(se$marital_status[,c(136,151)], 1, sum))
  se$marital_status$f_50_54mar_apart <- sqrt(apply(se$marital_status[,c(137,152)], 1, sum))
  se$marital_status$f_55_59mar_apart <- sqrt(apply(se$marital_status[,c(138,153)], 1, sum))
  se$marital_status$f_60_64mar_apart <- sqrt(apply(se$marital_status[,c(139,154)], 1, sum))
  se$marital_status$f_65_74mar_apart <- sqrt(apply(se$marital_status[,c(140,155)], 1, sum))
  se$marital_status$f_75_84mar_apart <- sqrt(apply(se$marital_status[,c(141,156)], 1, sum))
  se$marital_status$f_85up_mar_apart <- sqrt(apply(se$marital_status[,c(142,157)], 1, sum))
  
  est$marital_status <- 
    est$marital_status[, c(1,4,188,7:17,20,189,23:33,192:204,66,190,69:79,81,191,84:94,
                           97,205,100:110,113,206,116:126,209:221,159,207,162:172,174,208,177:187)]
  se$marital_status <- 
    se$marital_status[, c(1,4,188,7:17,20,189,23:33,192:204,66,190,69:79,81,191,84:94,
                           97,205,100:110,113,206,116:126,209:221,159,207,162:172,174,208,177:187)]
  
  names(se$marital_status) <- names(est$marital_status) <- c("cnt_mar_status_all", paste(
    rep(c("m", "f"), each= 65),
    rep(paste(rep(c("nvr_married", "married", "mar_apart", "widowed", "divorced"), each= 13),
    rep(c("15_17", "18_24", "25_29", "30_34", "35_39", "40_44", "45_49", "50_54", "55_59", "60_64",
          "65_74", "75_84", "85up"), 5), sep= "_"), 2), sep= "_"))
  
  ### edu
  #-----------
  est$edu <- est$edu[, -c(seq(3,42,8), seq(44,83,8))]
  est$edu <- est$edu[, c(1:2,38, 3:37, 39:73)]
  se$edu  <- se$edu[, -c(seq(3,42,8), seq(44,83,8))]
  se$edu  <- se$edu[, c(1:2,38, 3:37, 39:73)]
  
  names(est$edu) <- names(se$edu) <- c("cnt_edu_all", "edu_male", "edu_female", paste(
    rep(c("m", "f"), each= 35),
    paste(rep(c("18_24", "25_34", "35_44", "45_64", "65up"), each= 7),
          rep(c("lt_hs", "some_hs", "hs_grad", "some_col", "assoc_dec", "ba_deg", "grad_deg"), 5), 
    sep= "_"), sep= "_"))
  
  
  ### nativity
  #-----------
  est$nativity$cnt_u18 <- apply(est$nativity[,2:3],1,sum)
  est$nativity$born_st_u18 <- apply(est$nativity[,14:15],1,sum)
  est$nativity$born_out_st_u18 <- apply(est$nativity[,26:27],1,sum)
  est$nativity$born_out_us_u18 <- apply(est$nativity[,38:39],1,sum)
  est$nativity$foreigner_u18 <- apply(est$nativity[,50:51],1,sum)
  est$nativity$cnt_60_64 <- apply(est$nativity[,9:10],1,sum)
  est$nativity$born_st_60_64 <- apply(est$nativity[,21:22],1,sum)
  est$nativity$born_out_st_60_64 <- apply(est$nativity[,33:34],1,sum)
  est$nativity$born_out_us_60_64 <- apply(est$nativity[,45:46],1,sum)
  est$nativity$foreigner_60_64 <- apply(est$nativity[,57:58],1,sum)
  
  se$nativity$cnt_u18 <- sqrt(apply(se$nativity[,2:3],1,sum))
  se$nativity$born_st_u18 <- sqrt(apply(se$nativity[,14:15],1,sum))
  se$nativity$born_out_st_u18 <- sqrt(apply(se$nativity[,26:27],1,sum))
  se$nativity$born_out_us_u18 <- sqrt(apply(se$nativity[,38:39],1,sum))
  se$nativity$foreigner_u18 <- sqrt(apply(se$nativity[,50:51],1,sum))
  se$nativity$cnt_60_64 <- apply(se$nativity[,9:10],1,sum)
  se$nativity$born_st_60_64 <- apply(se$nativity[,21:22],1,sum)
  se$nativity$born_out_st_60_64 <- apply(se$nativity[,33:34],1,sum)
  se$nativity$born_out_us_60_64 <- apply(se$nativity[,45:46],1,sum)
  se$nativity$foreigner_60_64 <- apply(se$nativity[,57:58],1,sum)
  
  est$nativity <- est$nativity[, c(1,61,4:8,66,11:12,62,16:20,67,23:24,63,28:32,68,35:36,
                                   64,40:44,69,47:48,65,52:56,70,59:60)]
  se$nativity  <- se$nativity[, c(1,61,4:8,66,11:12,62,16:20,67,23:24,63,28:32,68,35:36,
                                  64,40:44,69,47:48,65,52:56,70,59:60)]
  
  names(est$nativity) <- names(se$nativity) <- c("cnt_nativity_all", paste(
    rep(c("cnt", "born_st_res", "born_out_state", "born_out_us", "foreigner"), each= 9),
    rep(c("u18", "18_24", "25_34", "35_44", "45_54", "55_59", "60_64", "65_74", "75up"), 5)
    , sep= "_"))
  
  ### by inc 12 mo
  #-----------
  est$by_inc_12mo <- est$by_inc_12mo[, c(1:2,4:11,13,15:22,24,26:33,35,37:44,46,48:55)]
  se$by_inc_12mo <-  se$by_inc_12mo[, c(1:2,4:11,13,15:22,24,26:33,35,37:44,46,48:55)]
  
  names(se$by_inc_12mo) <- names(est$by_inc_12mo) <- c("inc_nativity_total",paste(
    rep(c("cnts", "cit_born_st_res", "cit_born_other_st", "cit_born_out_us","foreign_born"), each= 9),
    rep(c("no_inc", "1_lt10k", "10k_lt15k", "15k_lt25k", "25k_lt35k", "35k_lt50k", 
          "50k_lt65k", "65k_lt75k", "gt75k"), 5), sep= "_"))
  
  ### by marital status
  #-----------
  # names(est$geo_mob_mar_stat) <- names(se$geo_mob_mar_stat) <- c("total", paste(
  #   rep(c("same_house", "same_cnty", "same_st", "diff_st", "abroad"), each= 6),
  #   rep(c("all", "nvr_married", "now_mar_exc_sep", "divorced", "separated",
  #         "widowed"), 5), sep= "_"))
  
  ### by edu
  #-----------
  names(est$geo_mob_edu) <- names(se$geo_mob_edu) <- c("total", paste(
    rep(c("same_house", "same_cnty", "same_st", "diff_st", "abroad"), each= 6),
    rep(c("all", "lt_hs", "high_sch", "some_col", "bachelors", "graduate"), 5), sep= "_"))
  
  ### by inc
  #-----------
  names(est$ind_inc) <- names(se$ind_inc) <- c("ind_inc_cnt_total", "no_inc", "has_inc", 
    "inc_lt10k", "inc10_lt15k", "inc15_lt25k", "inc25_lt35k", "inc35_lt50k", "inc50_lt65k",
    "inc65_lt75k", "inc75k_up")
  
  ### employment status
  #-----------
  est$emp_status$m16_19_empl <- apply(est$emp_status[, c(5,7)], 1, sum)
  est$emp_status$m20_24_empl <- apply(est$emp_status[, c(12,14,19,21)], 1, sum)
  est$emp_status$m20_24_unemp <- apply(est$emp_status[, c(15,22)], 1, sum)
  est$emp_status$m20_24_nolabor <- apply(est$emp_status[, c(16,23)], 1, sum)
  est$emp_status$m25_29_empl <- apply(est$emp_status[, c(26,28)], 1, sum)
  est$emp_status$m30_34_empl <- apply(est$emp_status[, c(33,35)], 1, sum)
  est$emp_status$m35_44_empl <- apply(est$emp_status[, c(40,42)], 1, sum)
  est$emp_status$m45_54_empl <- apply(est$emp_status[, c(47,49)], 1, sum)
  est$emp_status$m55_59_empl <- apply(est$emp_status[, c(54,56)], 1, sum)
  est$emp_status$m60_64_empl <- apply(est$emp_status[, c(61,63,68,70)], 1, sum)
  est$emp_status$m60_64_unemp <- apply(est$emp_status[, c(64,71)], 1, sum)
  est$emp_status$m60_64_nolabor <- apply(est$emp_status[, c(65,72)], 1, sum)
  est$emp_status$f16_19_empl <- apply(est$emp_status[, c(91,93)], 1, sum)
  est$emp_status$f20_24_empl <- apply(est$emp_status[, c(98,100,105,107)], 1, sum)
  est$emp_status$f20_24_unemp <- apply(est$emp_status[, c(101,108)], 1, sum)
  est$emp_status$f20_24_nolabor <- apply(est$emp_status[, c(102,109)], 1, sum)
  est$emp_status$f25_29_empl <- apply(est$emp_status[, c(112,114)], 1, sum)
  est$emp_status$f30_34_empl <- apply(est$emp_status[, c(119,121)], 1, sum)
  est$emp_status$f35_44_empl <- apply(est$emp_status[, c(126,128)], 1, sum)
  est$emp_status$f45_54_empl <- apply(est$emp_status[, c(133,135)], 1, sum)
  est$emp_status$f55_59_empl <- apply(est$emp_status[, c(140,142)], 1, sum)
  est$emp_status$f60_64_empl <- apply(est$emp_status[, c(147,149,154,156)], 1, sum)
  est$emp_status$f60_64_unemp <- apply(est$emp_status[, c(150,157)], 1, sum)
  est$emp_status$f60_64_nolabor <- apply(est$emp_status[, c(151,158)], 1, sum)
  
  se$emp_status$m16_19_empl <- sqrt(apply(se$emp_status[, c(5,7)], 1, sum))
  se$emp_status$m20_24_empl <- sqrt(apply(se$emp_status[, c(12,14,19,21)], 1, sum))
  se$emp_status$m20_24_unemp <- sqrt(apply(se$emp_status[, c(15,22)], 1, sum))
  se$emp_status$m20_24_nolabor <- sqrt(apply(se$emp_status[, c(16,23)], 1, sum))
  se$emp_status$m25_29_empl <- sqrt(apply(se$emp_status[, c(26,28)], 1, sum))
  se$emp_status$m30_34_empl <- sqrt(apply(se$emp_status[, c(33,35)], 1, sum))
  se$emp_status$m35_44_empl <- sqrt(apply(se$emp_status[, c(40,42)], 1, sum))
  se$emp_status$m45_54_empl <- sqrt(apply(se$emp_status[, c(47,49)], 1, sum))
  se$emp_status$m55_59_empl <- sqrt(apply(se$emp_status[, c(54,56)], 1, sum))
  se$emp_status$m60_64_empl <- sqrt(apply(se$emp_status[, c(61,63,68,70)], 1, sum))
  se$emp_status$m60_64_unemp <- sqrt(apply(se$emp_status[, c(64,71)], 1, sum))
  se$emp_status$m60_64_nolabor <- sqrt(apply(se$emp_status[, c(65,72)], 1, sum))
  se$emp_status$f16_19_empl <- sqrt(apply(se$emp_status[, c(91,93)], 1, sum))
  se$emp_status$f20_24_empl <- sqrt(apply(se$emp_status[, c(98,100,105,107)], 1, sum))
  se$emp_status$f20_24_unemp <- sqrt(apply(se$emp_status[, c(101,108)], 1, sum))
  se$emp_status$f20_24_nolabor <- sqrt(apply(se$emp_status[, c(102,109)], 1, sum))
  se$emp_status$f25_29_empl <- sqrt(apply(se$emp_status[, c(112,114)], 1, sum))
  se$emp_status$f30_34_empl <- sqrt(apply(se$emp_status[, c(119,121)], 1, sum))
  se$emp_status$f35_44_empl <- sqrt(apply(se$emp_status[, c(126,128)], 1, sum))
  se$emp_status$f45_54_empl <- sqrt(apply(se$emp_status[, c(133,135)], 1, sum))
  se$emp_status$f55_59_empl <- sqrt(apply(se$emp_status[, c(140,142)], 1, sum))
  se$emp_status$f60_64_empl <- sqrt(apply(se$emp_status[, c(147,149,154,156)], 1, sum))
  se$emp_status$f60_64_unemp <- sqrt(apply(se$emp_status[, c(150,157)], 1, sum))
  se$emp_status$f60_64_nolabor <- sqrt(apply(se$emp_status[, c(151,158)], 1, sum))
  
  est$emp_status <- est$emp_status[, c(1,174,8:9,175:178,29:30,179,36:37,180,43:44,181,50:51,
       182,57:58,183:185,75:77,80:82,85:87,
       186,94:95,187:190,115:116,191,122:123,192,129:130,193,136:137,194,143:144,195:197,
       161:163,166:168,171:173)]
  se$emp_status <- se$emp_status[, c(1,174,8:9,175:178,29:30,179,36:37,180,43:44,181,50:51,
       182,57:58,183:185,75:77,80:82,85:87,
       186,94:95,187:190,115:116,191,122:123,192,129:130,193,136:137,194,143:144,195:197,
       161:163,166:168,171:173)]
  
  names(est$emp_status) <- names(se$emp_status) <- c("emp_status_cnt_all", paste(
    rep(c("m", "f"), each= 33),
    paste(rep(c("16_19", "20_24", "25_29", "30_34", "35_44", "45_54", "55_59", "60_64",
                "65_69", "70_74", "75up"), each= 3),
          rep(c("employed", "unemp", "no_labor"), 11), sep= "_"), sep= "_"))
  
  ### pov status 1
  #-----------
  est$pov_status1 <- est$pov_status1[, c(1,5:7,10:12,16:18,21:23)]
  se$pov_status1  <- se$pov_status1[, c(1,5:7,10:12,16:18,21:23)]
  names(est$pov_status1) <- names(se$pov_status1) <- c("pov_status1_cnt_total", paste(
    rep(c("m_lt_pov", "f_lt_pov", "m_gt_eq_pov", "f_gt_eq_pov"), each=3),
    rep(c("employed", "unemployed", "not_in_labor_force"), 4), sep= "_"))
  
  ### pov status2
  #-----------
  # est$pov_status2$m_lt_15_17 <- apply(est$pov_status2[, 4:5], 1, sum)
  # est$pov_status2$f_lt_15_17 <- apply(est$pov_status2[, 14:15], 1, sum)
  # est$pov_status2$m_gt_15_17 <- apply(est$pov_status2[, 25:26], 1, sum)
  # est$pov_status2$f_gt_15_17 <- apply(est$pov_status2[, 35:36], 1, sum)
  # 
  # se$pov_status2$m_lt_15_17 <- sqrt(apply(se$pov_status2[, 4:5], 1, sum))
  # se$pov_status2$f_lt_15_17 <- sqrt(apply(se$pov_status2[, 14:15], 1, sum))
  # se$pov_status2$m_gt_15_17 <- sqrt(apply(se$pov_status2[, 25:26], 1, sum))
  # se$pov_status2$f_gt_15_17 <- sqrt(apply(se$pov_status2[, 35:36], 1, sum))
  # 
  # est$pov_status2 <- est$pov_status2[, c(1,44,6:12,45,16:22,46,27:33,47,37:43)]
  # se$pov_status2  <- se$pov_status2[, c(1,44,6:12,45,16:22,46,27:33,47,37:43)]
  # 
  # names(est$pov_status2) <- names(se$pov_status2) <-  c("pov_status1_cnt_total", paste(
  #   rep(c("m", "f", "m", "f"), each= 8),
  #   rep(c("lt_pov", "gt_eq_pov"), each= 16),
  #   rep(c("15_17", "18_24", "25_34", "35_44", "45_54", "55_64", "65_74", "75up")), sep= "_"))
  
  # 05 -- sort, combine and return
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
  
  return(ret)
  
}

