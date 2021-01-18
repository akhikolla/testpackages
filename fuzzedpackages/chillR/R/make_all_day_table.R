#' Fill in missing days in incomplete time series
#' 
#' Time series often have gaps, and these are often not marked by 'no data'
#' values but simply missing from the dataset. This function completes the time
#' series by adding lines for all these missing records. For these lines, all
#' values are set to 'NA'. By setting timestep<-"hour", this function can also
#' process hourly data. Where data are provided at a time resolution that is finer
#' than timestep, values are aggregated (by calculating the mean) to timestep
#' resolution (e.g. when data are at 15-minute resolution, they will be aggregated
#' to hourly average values - at timestep=="hour" - or daily average values - at
#' timestep=="day").
#' 
#' 
#' @param tab a data.frame containing a time series dataset. It should have
#' columns c("Year", "Month", "Day") or c("YEAR", "MONTH","DAY") or "YEARMODA".
#' @param timestep time step for the table. This defaults to 'day' but can also be 'hour'
#' @param input_timestep can also be 'day' or 'hour' and defaults to the value assigned
#' to timestep. If timestep is 'day' and input_timestep is 'hour', hourly records are
#' aggregated to daily Tmin, Tmean and Tmax.
#' @param tz timezone. Defaults to GMT. While it isn't important in what time zone
#' the temperatures were recorded, the onset of daylight savings time can cause problems.
#' 'GMT' is the correct setting in cases were the recorded times weren't adjusted
#' according to daylight savings time (i.e. no hours omitted or double-counted because
#' of such adjustment).
#' @param add.DATE boolean parameter indicating whether a column called DATE which
#' contains the IOSdate should be added to the output data.frame.
#' @param no_variable_check boolean parameter to indicate whether the function should
#' check if the dataset contains the usual chillR temperature variables. Defaults to
#' TRUE, but should be set to FALSE for different data formats.
#' @param aggregation_hours vector or list consisting of three integers that specify how
#' the function should search for daily minimum and maximum temperatures in hourly
#' datasets, when not all hourly temperatures have been observed. This is only relevant
#' during conversion from hourly to daily data. Tmin and Tmax can only be derived when
#' temperatures have been recorded during the coldest and warmest parts of the day,
#' respectively. The function should therefore check if records are available for these
#' times. The elements of `aggregation_hours` describe window sizes for the times (as
#' number of hours), during which the coldest and warmest temperature typically occurs.
#' The first two elements (which can be named `min_hours` and `max_hours`)
#' specify the number of hours contained in these windows for the cold and warm parts
#' of the day, respectively. These hours are determined by computing mean hourly
#' temperatures over the entire weather record, disaggregated by month to account for
#' the impact of daylength. The third element, `hours_needed` specifies how many
#' records during these windows have to have been recorded. `aggregation_hours`
#' defaults to NULL, in which case the parameter is ignored.
#' @return data frame containing all the columns of the input data frame, but
#' one row for each day between the start and end of the dataset. Data values
#' for the missing rows are filled in as 'NA'. Dates are expressed as
#' c("YEARMODA","DATE","Year","Month","Day"). In this, 'DATE' is the date in
#' ISOdate format.
#' @author Eike Luedeling
#' @references Luedeling E, Kunz A and Blanke M, 2013. Identification of
#' chilling and heat requirements of cherry trees - a statistical approach.
#' International Journal of Biometeorology 57,679-689.
#' @keywords utility
#' @examples
#' 
#' #fill in missing lines in a weather dataset (modified from KA_weather)
#' day_to_day<-make_all_day_table(KA_weather[c(1:10,20:30),],timestep="day")
#' 
#' #fill in missing hours in the Winters_hours_gaps dataset
#' Winters_hours<-subset(Winters_hours_gaps, select = -c(Temp_gaps))[1:2000,]
#' hour_to_hour<-make_all_day_table(Winters_hours,timestep="hour",input_timestep="hour")
#' 
#' #convert Winters_hours_gaps dataset into daily temperature data (min, max, mean)
#' hour_to_day<-make_all_day_table(Winters_hours,timestep="day",input_timestep="hour")
#' hour_to_day<-make_all_day_table(Winters_hours,timestep="day",input_timestep="hour",
#'                                aggregation_hours=c(3,3,2))
#' 
#' @export make_all_day_table
make_all_day_table<-function(tab,timestep="day",input_timestep=timestep,tz="GMT",
                             add.DATE=TRUE,no_variable_check=FALSE,aggregation_hours=NULL) #tab should have columns named Year, Month and Day (or YEAR, MONTH, DAY; or YEARMODA)
{
  only_numeric_mean<-function(x) if(is.numeric(x)) mean(x) else x[1]
  
  if(!timestep %in% c("day","hour"))
    return(warning("timestep needs to be 'hour' or 'day'"))
  if(!input_timestep %in% c("day","hour"))
    return(warning("input_timestep needs to be 'hour' or 'day'"))
  
  if(input_timestep=="hour")
    check<-check_temperature_record(tab,hourly=TRUE,no_variable_check=no_variable_check)
  if(input_timestep=="day")
    check<-check_temperature_record(tab,hourly=FALSE,no_variable_check=no_variable_check)
  
  if(check$dates_as_YEARMODA)
  {if(timestep=="hour")
  {tab[,"Year"]<-floor(tab$YEARMODAHO/1000000)
  tab[,"Month"]<-floor(tab$YEARMODAHO%%1000000/10000)
  tab[,"Day"]<-floor(tab$YEARMODAHO%%10000/100)
  tab[,"Hour"]<-tab$YEARMODAHO%%100} 
    if(timestep=="day")
    {tab[,"Year"]<-floor(tab$YEARMODA/10000)
    tab[,"Month"]<-floor(tab$YEARMODA%%10000/100)
    tab[,"Day"]<-tab$YEARMODA%%100}}
  
  
  if(timestep=="day"&input_timestep=="hour")
    convert_hour_to_day<-TRUE else convert_hour_to_day<-FALSE
    
    tab<-tab[which(!is.na(tab$Year)&!is.na(tab$Month)&!is.na(tab$Day)),]
    if(input_timestep=="hour")
      tab<-tab[which(!is.na(tab$Hour)),]
    
    if(timestep=="hour"&input_timestep=="day")
      return(warning("Converting daily to hourly data isn't possible with this function. ",
                     "Check out the interpolate_gaps_hourly function."))
    
    if(!check$chillR_compliant) return(warning("tab doesn't comply with chillR's data formatting standard"))
    
    if(!convert_hour_to_day)
      if(check$chillR_compliant&
         check$completeness_check$missing_records==0&
         check$completeness_check$repeated_records==0)
      {if((!"DATE" %in% colnames(tab))&add.DATE) 
        if(timestep=="hour")
          tab[,"DATE"]<-ISOdatetime(tab$Year,tab$Month,tab$Day,tab$Hour,0,0,tz=tz) else
            tab[,"DATE"]<-ISOdate(tab$Year,tab$Month,tab$Day)
          tab<-tab[,c("DATE",colnames(tab)[which(!colnames(tab)=="DATE")])]
          return(tab) #maybe some format changes?
      }
    
    if(input_timestep=="hour")
    {alltimes<-ISOdatetime(tab$Year,tab$Month,tab$Day,tab$Hour,0,0,tz=tz)
    datevec<-seq(min(alltimes,na.rm=TRUE),max(alltimes,na.rm=TRUE), "hour")}
    if(input_timestep=="day")
    {alltimes<-ISOdate(tab$Year,tab$Month,tab$Day)
    datevec<-seq(min(alltimes,na.rm=TRUE),max(alltimes,na.rm=TRUE),"DSTday")}
    
    outcols<-colnames(tab)
    
    tab[,"Alltime"]<-alltimes
    
    if(check$completeness_check$repeated_records>0)
    {reps<-which(table(as.character(alltimes))>1)
    oks<-subset(tab,is.na(match(as.character(alltimes),names(reps))))
    notoks<-subset(tab,!is.na(match(as.character(alltimes),names(reps))))
    notokfixed<-aggregate(notoks,by=list(notoks$Alltime),FUN=only_numeric_mean)
    notokfixed<-notokfixed[,colnames(notoks)]
    tab<-rbind(oks,notokfixed)
    tab<-tab[order(tab$Alltime),]
    alltimes<-tab$Alltime
    }
    
    if(convert_hour_to_day)
    {
      alltimes<-ISOdate(tab$Year,tab$Month,tab$Day)
    tab$Alltime<-alltimes
    datevec<-seq(min(alltimes),max(alltimes),"DSTday")      
    
    if(!is.null(aggregation_hours))
      {
        if(length(aggregation_hours)==3)
          {if("min_hours" %in% names(aggregation_hours)&"max_hours" %in% names(aggregation_hours)&"hours_needed" %in% names(aggregation_hours))
          {min_hours<-as.numeric(aggregation_hours["min_hours"])
           max_hours<-as.numeric(aggregation_hours["max_hours"])
           hours_needed<-as.numeric(aggregation_hours["hours_needed"])
          } else
          {min_hours<-aggregation_hours[1]
          max_hours<-aggregation_hours[2]
          hours_needed<-aggregation_hours[3]
          }} else return(warning("not clear how many values need to be present for deriving temperature extremes from incomplete hourly records"))

      mean_hour_temps<-aggregate(tab,by=list(tab$Hour,tab$Month),function(x) mean(x,na.rm=TRUE))[,c("Month","Hour","Temp")]
      month_min_hours<-list()
      month_max_hours<-list()
      for(m in 1:12)
        {month_min_hours[[m]]<-mean_hour_temps$Hour[which(mean_hour_temps$Month==m)][order(mean_hour_temps$Temp[which(mean_hour_temps$Month==m)])][1:min_hours]
         month_max_hours[[m]]<-mean_hour_temps$Hour[which(mean_hour_temps$Month==m)][order(mean_hour_temps$Temp[which(mean_hour_temps$Month==m)],decreasing = TRUE)][1:max_hours]}
        
      means<-aggregate(tab,by=list(tab$Alltime),FUN=mean)
      colnames(means)[which(colnames(means)=="Temp")]<-"Tmean"
      means[,"Tmin"]<-do.call(rbind,lapply(split(tab,tab$Alltime),function(chunk) if(length(which(chunk$Hour[which(!is.na(chunk$Temp))] %in% month_min_hours[[mean(chunk$Month)]]))>=hours_needed)
        min(chunk$Temp,na.rm=TRUE) else (NA)))
      means[,"Tmax"]<-do.call(rbind,lapply(split(tab,tab$Alltime),function(chunk) if(length(which(chunk$Hour[which(!is.na(chunk$Temp))] %in% month_max_hours[[mean(chunk$Month)]]))>=hours_needed)
        max(chunk$Temp,na.rm=TRUE) else (NA)))
    }
 
    if(is.null(aggregation_hours))
      {means<-aggregate(tab,by=list(tab$Alltime),FUN=only_numeric_mean)
      colnames(means)[which(colnames(means)=="Temp")]<-"Tmean"
      means[,"Tmin"]<-aggregate(tab[,c("Alltime","Temp")],by=list(tab$Alltime),FUN=min)[,"Temp"]
      means[,"Tmax"]<-aggregate(tab[,c("Alltime","Temp")],by=list(tab$Alltime),FUN=max)[,"Temp"]}

    outcols<-c(outcols[which(!outcols %in% c("Tmin","Tmean","Tmax","Temp","Hour"))],"Tmin","Tmean","Tmax")
    tab<-means
    alltimes<-tab$Alltime
    }
    
    
    #   matches<-match(datevec,alltimes)
    
    dv<-data.frame(Date=datevec)
    dv$Year<-as.numeric(format(datevec, "%Y"))
    dv$Month<-as.numeric(format(datevec, "%m"))
    dv$Day<-as.numeric(format(datevec, "%d"))
    if(timestep=="hour") dv$Hour<-as.numeric(format(datevec, "%H"))
    
    for(cc in colnames(tab)[which(is.na(match(colnames(tab),c("Year","Month","Day","Hour","Group.1","Date"))))])
      dv[which(datevec %in% alltimes),cc]<-tab[which(alltimes %in% datevec),cc]
    
    if(add.DATE)
    {dv[,"DATE"]<-dv[,"Date"]
    outcols<-c("DATE",outcols)}
    
    return(dv[,outcols])
}
