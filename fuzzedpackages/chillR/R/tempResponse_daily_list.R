#' Calculation of climatic metrics from lists of daily temperature records
#' 
#' Wrapper for the tempResponse function, to facilitate its use on lists
#' of daily temperature records, e.g. those produced by the 
#' \code{\link[=temperature_generation]{temperature_generation}}
#' function. Daily temperature records are converted
#' into hourly records using either the 
#' \code{\link[=stack_hourly_temps]{stack_hourly_temps}}
#' function or an empirical relationship between observed hourly temperatures
#' and daily temperature extremes (see
#' \code{\link[=Empirical_hourly_temperatures]{Empirical_hourly_temperatures}}
#' for details). These hourly
#' records are then used as input into the
#' \code{\link[=tempResponse]{tempResponse}} function, to which
#' most parameters are passed. See the documentation of
#' \code{\link[=tempResponse]{tempResponse}} for
#' more details.
#' 
#' @param temperature_list list of daily temperature records, as produced
#' by \code{\link[=temperature_generation]{temperature_generation}}.
#' @param latitude latitude of the location of interest (used for generating
#' hourly records).
#' @param Start_JDay the start date (in Julian date, or day of the year) of the
#' period, for which chill and heat should be quantified.
#' @param End_JDay the end date (in Julian date, or day of the year) of the
#' period, for which chill and heat should be quantified.
#' @param models named list of models that should be applied to the hourly
#' temperature data. These should be functions that take as input a vector of
#' hourly temperatures. This defaults to the set of models provided by the
#' chilling function.
#' @param misstolerance maximum percentage of values for a given season that
#' can be missing without the record being removed from the output. Defaults to
#' 50.
#' @param whole_record boolean parameter indicating whether the metrics should
#' be summed over the entire temperature record. If set to \code{TRUE} (default is
#' \code{FALSE}), then the function ignores the specified start and end dates and
#' simply returns the totals of each metric that accumulated over the entire
#' temperature record.
#' @param empirical indicates whether hourly temperatures should be generated
#' based on an idealized temperature curve (set to \code{NULL}, the default) or an
#' empirically derived relationship between hourly temperatures and daily
#' temperature extremes (see
#' \code{\link[=Empirical_hourly_temperatures]{Empirical_hourly_temperatures}} and
#' \code{\link[=Empirical_daily_temperature_curve]{Empirical_daily_temperature_curve}},
#' also for the format of the empirical prediction coefficient data.frame). If the
#' latter, this parameter needs to be a data.frame including columns \code{Month},
#' \code{Hour} and \code{Prediction_coefficients}. See
#' \code{\link[=Empirical_daily_temperature_curve]{Empirical_daily_temperature_curve}}
#' for further details on the format.
#' @param mean_out boolean parameter indicating whether the mean of the input
#' metric (e.g. temperature) should be returned in a column named "Input_mean".
#' @return data frame showing totals for all specified models for the
#' respective periods for all seasons included in the temperature records.
#' Columns are \code{Season}, \code{End_year} (the year when the period ended)
#' and \code{Days} (the duration of the period), as well as one column per model,
#' which receives the same name as the function in the models list.
#' If the weather input consisted of a list with elements \code{stack} and \code{QC},
#' the output also contains columns from
#' \code{QC} that indicate the completeness of the weather record that the
#' calculations are based on.
#' @author Eike Luedeling
#' @references The chillR package:
#' 
#' Luedeling E, Kunz A and Blanke M, 2013. Identification of chilling and heat
#' requirements of cherry trees - a statistical approach. International Journal
#' of Biometeorology 57,679-689.
#' @keywords chill and heat calculation
#' @examples
#' 
#' 
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2006),])
#' temperature_list<-list(weather,weather,weather)
#' 
#' tempResponse_daily_list(temperature_list,latitude=50.4)
#' 
#' 
#' @export tempResponse_daily_list
tempResponse_daily_list <-
function (temperature_list,latitude,Start_JDay=1,End_JDay=366,
          models=list(Chilling_Hours=Chilling_Hours,Utah_Chill_Units=Utah_Model,Chill_Portions=Dynamic_Model,GDH=GDH),
          misstolerance=50,whole_record=FALSE,empirical=NULL,mean_out=FALSE)             
     {
  
  if(is.data.frame(temperature_list))
    temperature_list<-list(temperature_list)
  
  output<-list()
  
  for(i in 1:length(temperature_list))
    {
    if(is.null(empirical))  
      hourtemps<-stack_hourly_temps(temperature_list[[i]],latitude=latitude)
    # more checks needed here
    if(is.data.frame(empirical))
      hourtemps<-Empirical_hourly_temperatures(temperature_list[[i]],empirical)
     output[[i]]<-tempResponse(hourtemps,Start_JDay=Start_JDay,End_JDay=End_JDay,
                            models=models,misstolerance=misstolerance,
                            whole_record=whole_record,mean_out=mean_out)}
  names(output)<-names(temperature_list)
  return(output)
}


