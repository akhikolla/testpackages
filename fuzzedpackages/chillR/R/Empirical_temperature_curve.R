#' Empirical daily temperature curve
#' 
#' This function derives an empirical daily temperature curve from observed hourly
#' temperature data. The mean temperature during each hour of the day is expressed
#' as a function of the daily minimum and maximum temperature. This is done
#' separately for each month of the year. The output is a data.frame that can then
#' be used with the
#' \code{\link[=Empirical_hourly_temperatures]{Empirical_hourly_temperatures}}
#' function to generate hourly temperatures from data on daily minimum (\code{Tmin})
#' and maximum (\code{Tmax}) temperatures.
#' 
#' 
#' @param Thourly data.frame containing hourly temperatures. Must contain columns
#' \code{Year} (year of observation), \code{Month} (month of observation), \code{Day} (day of 
#' observation), \code{Hour} (hour of observation) and \code{Temp} (Observed temperature).
#' If multiple observations within an hour are available, these are averaged.
#' 
#' @import dplyr
#' @importFrom rlang .data
#' 
#' @return data.frame containing three columns: \code{Month} (month for which coefficient
#' applies), \code{Hour} (hour for which coefficient applies) and \code{Prediction_coefficient}
#' (the coefficient used for empirical temperature prediction). Coefficients
#' indicate, by what fraction of the daily temperature range the temperature during
#' the specified hour is above the daily minimum temperature.
#' 
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#'  
#' Empirical_daily_temperature_curve(Winters_hours_gaps)
#' 
#' @export Empirical_daily_temperature_curve
Empirical_daily_temperature_curve<-function(Thourly)
{
  ### checks to ensure that the dataset is suitable for this analysis.
  assertthat::assert_that(is.data.frame(Thourly),msg="Thourly not a data.frame")
  assertthat::assert_that("Year" %in% colnames(Thourly),msg="Column 'Year' missing")
  assertthat::assert_that("Month" %in% colnames(Thourly),msg="Column 'Month' missing")
  assertthat::assert_that("Day" %in% colnames(Thourly),msg="Column 'Day' missing")
  assertthat::assert_that("Hour" %in% colnames(Thourly),msg="Column 'Hour' missing")
  assertthat::assert_that("Temp" %in% colnames(Thourly),msg="Column 'Temp' missing")
  assertthat::assert_that(is.numeric(Thourly$Year),msg="'Year' column not numeric")
  assertthat::assert_that(is.numeric(Thourly$Month),msg="'Month' column not numeric")
  assertthat::assert_that(is.numeric(Thourly$Day),msg="'Day' column not numeric")
  assertthat::assert_that(is.numeric(Thourly$Temp),msg="'Temp' column not numeric")
  
 
  ## derive daily extremes from (sub)hourly data
  Thourly[,"YEARMODA"]<-Thourly$Year*10000+Thourly$Month*100+Thourly$Day
  Tday <- Thourly %>% group_by(.data$Year, .data$Month, .data$Day) %>% summarise(Tmin = min(.data$Temp),
                                                               Tmax = max(.data$Temp))
  
  # summarize (sub)hourly data into hourly data (using the mean of subhourly values)
  Thours <- Thourly %>% group_by(.data$Year, .data$Month, .data$Day,.data$Hour) %>% summarise(Temp = mean(.data$Temp))
  
  merged<-merge(Thours,Tday)
  # temperatures are scaled by the daily temperature range, i.e. we compute how much of the overall
  # daytime warming (or cooling) has occurred by a particular hour of the day. A possible improvement
  # would be considering the temperatures of the previous day, but that seems to complicated right now...
  scaled<-merged
  scaled[,"Tscaled"]<-(scaled$Temp-scaled$Tmin)/(scaled$Tmax-scaled$Tmin)
  
  scaled_summ<-scaled %>% group_by(.data$Month,.data$Hour) %>% summarise(Tscaled = mean(.data$Tscaled))
  
  for (mm in unique(scaled_summ$Month))
  {
    sel<-which(scaled_summ$Month==mm)
    scaled_summ[sel,"Tscale_adj"]<-scaled_summ[sel,"Tscaled"]-min(scaled_summ[sel,"Tscaled"])
    scaled_summ[sel,"Tscale_adj"]<-scaled_summ[sel,"Tscale_adj"]/max(scaled_summ[sel,"Tscale_adj"])
  }
  
  empi_coeffs<-as.data.frame(scaled_summ[,c(1,2,4)])
  colnames(empi_coeffs)[3]<-"Prediction_coefficient"
  
  return(empi_coeffs)
}



#' Empirical daily temperature prediction
#' 
#' This function generates hourly temperatures from daily minimum and maximum
#' temperatures, based on an empirical relationship of these two daily 
#' temperature extremes with the hourly temperature. Usually, this relationship
#' will have been determined with the
#' \code{\link[=Empirical_daily_temperature_curve]{Empirical_daily_temperature_curve}}
#' function.
#' 
#' @param Tdaily data.frame containing daily minimum and maximum temperatures.
#' Must contain columns \code{Year} (year of observation), \code{Month} (month of observation),
#' \code{Day} (day of observation), \code{Tmin} (Minimum daily temperature) and \code{Tmax}
#' (Maximum daily temperature).
#' @param empi_coeffs data.frame containing coefficients for the hourly temperature
#' prediction, e.g. generated with the function
#' \code{\link[=Empirical_daily_temperature_curve]{Empirical_daily_temperature_curve}}.
#' Needs to contain the following columns: \code{Month} (month for which coefficient
#' applies), \code{Hour} (hour for which coefficient applies) and \code{Prediction_coefficient}
#' (the coefficient used for empirical temperature prediction). Coefficients
#' indicate, by what fraction of the daily temperature range the temperature during
#' the specified hour is above the daily minimum temperature.
#' 
#' @import dplyr
#' 
#' @return data.frame containing all columns of the \code{Tdaily} dataset, but also the
#' columns \code{Hour} and \code{Temp}, for the hour of the day and the predicted temperature,
#' respectively.
#' 
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#'  
#' coeffs<-Empirical_daily_temperature_curve(Winters_hours_gaps)
#' Winters_daily<-make_all_day_table(Winters_hours_gaps, input_timestep="hour")
#' Empirical_hourly_temperatures(Winters_daily,coeffs)
#' 
#' @export Empirical_hourly_temperatures
Empirical_hourly_temperatures<-function(Tdaily,empi_coeffs)
{
  ##checks to ensure that the inputs are suitable
  assertthat::assert_that(is.data.frame(Tdaily),msg="Tdaily not a data.frame")
  assertthat::assert_that("Year" %in% colnames(Tdaily),msg="Tdaily: Column 'Year' missing")
  assertthat::assert_that("Month" %in% colnames(Tdaily),msg="Tdaily: Column 'Month' missing")
  assertthat::assert_that("Day" %in% colnames(Tdaily),msg="Tdaily: Column 'Day' missing")
  assertthat::assert_that("Tmin" %in% colnames(Tdaily),msg="Tdaily: Column 'Tmin' missing")
  assertthat::assert_that("Tmax" %in% colnames(Tdaily),msg="Tdaily: Column 'Tmax' missing")
  assertthat::assert_that(is.numeric(Tdaily$Year),msg="Tdaily: 'Year' column not numeric")
  assertthat::assert_that(is.numeric(Tdaily$Month),msg="Tdaily: 'Month' column not numeric")
  assertthat::assert_that(is.numeric(Tdaily$Day),msg="Tdaily: 'Day' column not numeric")
  assertthat::assert_that(is.numeric(Tdaily$Tmin),msg="Tdaily: 'Tmin' column not numeric")
  assertthat::assert_that(is.numeric(Tdaily$Tmax),msg="Tdaily: 'Tmax' column not numeric")  
  assertthat::assert_that(is.data.frame(empi_coeffs),msg="empi_coeffs not a data.frame")
  assertthat::assert_that("Month" %in% colnames(empi_coeffs),msg="empi_coeffs: Column 'Month' missing")
  assertthat::assert_that("Hour" %in% colnames(empi_coeffs),msg="empi_coeffs: Column 'Hour' missing")
  assertthat::assert_that("Prediction_coefficient" %in% colnames(empi_coeffs),msg="Column 'Prediction_coefficient' missing")
  assertthat::assert_that(is.numeric(empi_coeffs$Month),msg="empi_coeffs: 'Month' column not numeric")
  assertthat::assert_that(is.numeric(empi_coeffs$Hour),msg="empi_coeffs: 'Hour' column not numeric")
  assertthat::assert_that(is.numeric(empi_coeffs$Prediction_coefficient),msg="empi_coeffs: 'Prediction_coefficient' column not numeric")
  
  frame<-stack_hourly_temps(as.data.frame(Tdaily),latitude=0)$hourtemps
  
  frame[,"YEARMODAHO"]<-frame$Year*1000000+frame$Month*10000+frame$Day*100+frame$Hour
  frame[,"MonthHour"]<-frame$Month*100+frame$Hour
  empi_coeffs[,"MonthHour"]<-empi_coeffs$Month*100+empi_coeffs$Hour
  merged<-merge(frame,empi_coeffs[,c("MonthHour","Prediction_coefficient")],by="MonthHour")
  merged<-merged[order(merged$YEARMODAHO),]
  merged[,"Temp_empirical"]<-merged$Tmin+(merged$Tmax-merged$Tmin)*merged$Prediction_coefficient
  
  colnames(merged)[which(colnames(merged)=="Temp")]<-"Temp_idealized"
  merged[,"Temp"]<-merged$Temp_empirical
  
  merged<-merged[,which(!colnames(merged) %in% c("MonthHour","YEARMODAHO","Temp_idealized",
                                                 "Temp_empirical","Prediction_coefficient"))]
  return(merged)
}


