#' Check a daily or hourly temperature record for compliance with chillR's standards
#'
#' This function performs basic tests to determine whether a temperature record
#' complies with chillR's formatting rules. If desired, the function also checks
#' whether the record is complete (has rows for all time units in the interval)
#' and how many values are missing.
#'
#' @param weather object to be tested for whether it contains chillR-compatible
#' temperature data. 
#' @param hourly boolean parameter indicating whether temp_record contains hourly
#' data. If not, it is assumed to consist of daily records (the default).
#' @param completeness_check boolean parameter indicating whether the records should be
#' checked for completeness.
#' @param no_variable_check boolean parameter to indicate whether the function should
#' check if the dataset contains the usual chillR temperature variables. Defaults to
#' TRUE, but should be set to FALSE for different data formats.
#' @return list containing the following elements: 'data_frequency' ("daily or "hourly),
#' 'weather_object' (boolean, indicates whether records are in a sub-object called weather),
#' 'chillR_compliant' (boolean, indicates whether the object was found to conform to
#' chillR format standards) and 'error' (contains error messages generated during the
#' checking procedure).
#' @note This function doesn't check whether there are faulty data. It only tests
#' whether the data is compatible with the requirements of chillR's major functions.
#' 
#' @author Eike Luedeling
#' @references The chillR package:
#'
#' Luedeling E, Kunz A and Blanke M, 2013. Identification of chilling and heat
#' requirements of cherry trees - a statistical approach. International Journal
#' of Biometeorology 57,679-689.
#' @keywords utilities
#' @examples
#'
#' check_temperature_record(KA_weather)
#'
#' @export check_temperature_record
check_temperature_record<-function(weather,hourly=FALSE,completeness_check=TRUE,no_variable_check=FALSE)
{
  if(hourly) out<-list(data_frequency="hourly",weather_object=FALSE,
                       chillR_compliant=FALSE,error="none",dates_as_YEARMODA=FALSE)
  if(!hourly) out<-list(data_frequency="daily",weather_object=FALSE,
                        chillR_compliant=FALSE,error="none",dates_as_YEARMODA=FALSE)
  
  if(length(weather)==0) {out$error<-"no data provided (weather is null)"
  warning(paste("Error -",out$error))
  return(out)}
  
  if(!is.data.frame(weather))
  {if(is.null(names(weather)))
  {out$error<-"no weather data.frame found"
  warning(paste("Error -",out$error))
  return(out)} else
  {if(hourly)
  {if("hourtemps" %in% names(weather))
  {weather<-weather$hourtemps
  out$weather_object<-TRUE} else
  {out$error<-"weather is no data.frame and no argument 'hourfile' was found"
  warning(paste("Error -",out$error))
  return(out)}
  } else
  {if("weather" %in% names(weather))
  {weather<-weather$weather
  out$weather_object<-TRUE} else
  {out$error<-"weather is no data.frame and no argument 'weather' was found"
  warning(paste("Error -",out$error))
  return(out)}
  }
  }
  }
  
  if(!hourly)
    if(!no_variable_check) required_columns<-c("Year","Month","Day","Tmin","Tmax") else
      required_columns<-c("Year","Month","Day")
    if(hourly)
      if(!no_variable_check) required_columns<-c("Year","Month","Day","Hour","Temp") else
        required_columns<-c("Year","Month","Day","Hour")
      
      
      miss_cols<-c()
      for(rc in required_columns) if(!rc %in% colnames(weather)) miss_cols<-c(miss_cols,rc)
      
      if(!hourly)
        if("YEARMODA" %in% colnames(weather)&
           "Year" %in% miss_cols&"Month" %in% miss_cols&"Day" %in% miss_cols)
        {out$dates_as_YEARMODA<-TRUE
        weather[,"Year"]<-floor(weather$YEARMODA/10000)
        weather[,"Month"]<-floor(weather$YEARMODA%%10000/100)
        weather[,"Day"]<-weather$YEARMODA%%100
        miss_cols<-miss_cols[which(!(miss_cols %in% c("Year","Month","Day")))]}
      
      if(hourly)
        if("YEARMODAHO" %in% colnames(weather)&
           "Year" %in% miss_cols&"Month" %in% miss_cols&"Day" %in% miss_cols&"Hour" %in% miss_cols)
        {out$dates_as_YEARMODA<-TRUE
        weather[,"Year"]<-floor(weather$YEARMODAHO/1000000)
        weather[,"Month"]<-floor(weather$YEARMODAHO%%1000000/10000)
        weather[,"Day"]<-floor(weather$YEARMODAHO%%10000/100)
        weather[,"Hour"]<-weather$YEARMODAHO%%100
        miss_cols<-miss_cols[which(!(miss_cols %in% c("Year","Month","Day","Hour")))]} 
      
      
      if(length(miss_cols)>0)
      {if(length(miss_cols)==1)
        out$error<-paste("Column missing:",toString(miss_cols)) else
          out$error<-paste("Columns missing:",toString(miss_cols))
        warning(paste("Error -",out$error))
        return(out)
      }
      
      nnum_cols<-c()
      for(rc in required_columns) if(!is.numeric(weather[[rc]])) nnum_cols<-c(nnum_cols,rc)
      if(length(nnum_cols)>0)
      {if(length(nnum_cols)==1)
        out$error<-paste("One column is not numeric:",toString(nnum_cols)) else
          out$error<-paste("The following columns are not numeric:",toString(nnum_cols))
        warning(paste("Error -",out$error))
        return(out)
      }
      
      if(completeness_check)
      {
        #check if sufficient records
        
        if(hourly) alltimes<-ISOdatetime(weather$Year,weather$Month,weather$Day,weather$Hour,
                                         0,0,tz="GMT") else
                                           alltimes<-ISOdate(weather$Year,weather$Month,weather$Day)                                      
                                         
                                         
        if(!alltimes[1]<alltimes[length(alltimes)])
           {out$error<-paste("The first date of the record must be before the last date")
            warning(paste("Error -",out$error))
            return(out)
           }
                                         
        if(hourly) datevec<-seq(alltimes[1],alltimes[length(alltimes)], "hour") else
           datevec<-seq(alltimes[1],alltimes[length(alltimes)],"DSTday")
                                         
        missing_records<-length(which(is.na(match(datevec,alltimes))))
        repeated_records<-length(which(table(alltimes)>1))
        
        if(hourly)
           {if(no_variable_check)
               total_missing_temps<-missing_records else
                 total_missing_temps<-missing_records+length(which(is.na(weather$Temp)))
            out[["completeness_check"]]<-
               list(missing_records=missing_records,
                    repeated_records=repeated_records,
                    total_missing_Temps=total_missing_temps)
            }  else
            {if(no_variable_check)
               {total_missing_Tmin<-missing_records
                total_missing_Tmax<-missing_records} else
                   {total_missing_Tmin<-missing_records+length(which(is.na(weather$Tmin)))
                    total_missing_Tmax<-missing_records+length(which(is.na(weather$Tmax)))}
             out[["completeness_check"]]<-
               list(missing_records=missing_records,
                    repeated_records=repeated_records,
                    total_missing_Tmin=total_missing_Tmin,
                    total_missing_Tmax=total_missing_Tmax)
            }
      }
      out$chillR_compliant<-TRUE
      return(out)
}

