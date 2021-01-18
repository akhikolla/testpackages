#' Make hourly temperature record from daily data
#' 
#' This function generates hourly temperature records for a particular location
#' from daily minimum and maximum temperatures and latitude.
#' 
#' Temperature estimates are based on an idealized daily temperature curve that
#' uses a sine curve for daytime warming and a logarithmic decay function for
#' nighttime cooling. The input data frame can have more columns, which are
#' preserved, but ignored in the processing. References to papers outlining the
#' procedures are given below.
#' 
#' Note that this function should be able to generate hourly temperatures for
#' all latitudes, but it uses an algorithm designed for locations with regular
#' day/night behavior. It may therefore be that the curves aren't very realistic
#' for very short or very long days, or especially for polar days and nights.
#' 
#' @param latitude the geographic latitude (in decimal degrees) of the location
#' of interest
#' @param year_file a data frame containing data on daily minimum temperature
#' (called Tmin), daily maximum temperature (called Tmax), and date
#' information. Dates can either be specified by two columns called Year and
#' JDay, which contain the Year and Julian date (day of the year), or as three
#' columns called Year, Month and Day. year_file cannot have any missing
#' values, so it may be a good idea to process the relevant columns with
#' make_all_day_table and interpolate_gaps before.
#' @param keep_sunrise_sunset boolean variable indicating whether information
#' on sunrise, sunset and daylength, which is calculated for producing hourly
#' temperature records, should be preserved in the output. Defaults to FALSE.
#' @return data frame containing all the columns of year_file, plus 24 columns
#' for hourly temperatures (called Hour_1 ... Hour_24).
#' @author Eike Luedeling
#' @references Luedeling E, Kunz A and Blanke M, 2013. Identification of
#' chilling and heat requirements of cherry trees - a statistical approach.
#' International Journal of Biometeorology 57,679-689.
#' 
#' Luedeling E, Girvetz EH, Semenov MA and Brown PH, 2011. Climate change
#' affects winter chill for temperate fruit and nut trees. PLoS ONE 6(5),
#' e20155.
#' 
#' The temperature interpolation is described in
#' 
#' Linvill DE, 1990. Calculating chilling hours and chill units from daily
#' maximum and minimum temperature observations. HortScience 25(1), 14-16.
#' 
#' Calculation of sunrise, sunset and daylength was done according to
#' 
#' Spencer JW, 1971. Fourier series representation of the position of the Sun.
#' Search 2(5), 172.
#' 
#' Almorox J, Hontoria C and Benito M, 2005. Statistical validation of
#' daylength definitions for estimation of global solar radiation in Toledo,
#' Spain. Energy Conversion and Management 46(9-10), 1465-1471)
#' @keywords utility
#' @examples
#' 
#' weather<-fix_weather(KA_weather)
#' 
#' THourly<-make_hourly_temps(50.4,weather$weather)
#' 
#' #in most cases, you're probably better served by stack_hour_temperatures
#' 
#' @export make_hourly_temps
make_hourly_temps <-
  function (latitude,year_file,keep_sunrise_sunset=FALSE)
    
  {
    
    if(missing(latitude)) stop("'latitude' not specified")
    if(length(latitude)>1) stop("'latitude' has more than one element")
    if(!is.numeric(latitude)) stop("'latitude' is not numeric")
    if(latitude>90|latitude<(-90)) warning("'latitude' is usually between -90 and 90")
    
    
    year_file<-year_file[which(!is.na(year_file$Tmin)&!is.na(year_file$Tmax)),]
    
    if(!"JDay" %in% colnames(year_file))
      year_file[,"JDay"]<-strptime(paste(year_file$Month,"/",year_file$Day,"/",year_file$Year,sep=""),"%m/%d/%Y")$yday+1
    
    preserve_columns<-colnames(year_file)
    
    Day_times<-daylength(latitude=latitude,
                         JDay=c(year_file$JDay[1]-1,
                                year_file$JDay,
                                year_file$JDay[nrow(year_file)]+1))
    Day_times$Sunrise[which(Day_times$Sunrise==99)]<-0
    Day_times$Sunrise[which(Day_times$Sunrise==-99)]<-12
    Day_times$Sunset[which(Day_times$Sunset==99)]<-24
    Day_times$Sunset[which(Day_times$Sunset==-99)]<-12
    
    year_file$Sunrise<-Day_times$Sunrise[2:(length(Day_times$Sunrise)-1)]
    year_file$Sunset<-Day_times$Sunset[2:(length(Day_times$Sunset)-1)]
    year_file$Daylength<-Day_times$Daylength[2:(length(Day_times$Daylength)-1)]
    year_file$prev_Sunset<-Day_times$Sunset[1:(length(Day_times$Sunset)-2)]
    year_file$next_Sunrise<-Day_times$Sunrise[3:length(Day_times$Sunrise)]
    year_file$prev_max<-year_file$Tmax[c(NA,1:(nrow(year_file)-1))]
    year_file$next_min<-year_file$Tmin[c(2:nrow(year_file),NA)]
    year_file$prev_min<-year_file$Tmin[c(NA,1:(nrow(year_file)-1))]
    year_file$Tsunset<-year_file$Tmin+(year_file$Tmax-year_file$Tmin)*
      sin((pi*(year_file$Sunset-year_file$Sunrise)/(year_file$Daylength+4)))
    year_file$prev_Tsunset<-year_file$prev_min+(year_file$prev_max-year_file$prev_min)*
      sin((pi*(year_file$Daylength)/(year_file$Daylength+4)))
    colnum<-ncol(year_file)+1
    
    hourcol<-c(colnum:(colnum+23))
    
    
    
    
    for (hour in 0:23)
    {
      hourcount<-hour+1
      
      #if(length(which(year_file$Daylength==-99))>0)
      #{
      no_riseset<-which(year_file$Daylength %in% c(0,24,-99))
      year_file[no_riseset,colnum+hour]<-((year_file$Tmax+year_file$Tmin)/2)[no_riseset]
      #}
      
      c_morn<-which(hour<=year_file$Sunrise)
      if(1 %in% c_morn)
        if(!length(c_morn)==1) 
          c_morn<-c_morn[2:length(c_morn)]
      else c_morn<-c() #can't compute temperatures before sunrise for day 1
      c_day<-which(hour>year_file$Sunrise&hour<=year_file$Sunset)
      c_eve<-which(hour>=year_file$Sunset)
      if(nrow(year_file) %in% c_eve) c_eve<-c_eve[1:(length(c_eve)-1)] #can't compute temperatures after sunset for last day
      
      
      year_file[c_morn,colnum+hour]<-
        year_file$prev_Tsunset[c_morn]-  #prev temp at sunset
        ((year_file$prev_Tsunset[c_morn]-year_file$Tmin[c_morn])/
           log(max(1,24-(year_file$prev_Sunset[c_morn]-year_file$Sunrise[c_morn])))*
           log(hour+24-year_file$prev_Sunset[c_morn]+1))
      
      year_file[c_day,colnum+hour]<-
        year_file$Tmin[c_day]+
        (year_file$Tmax[c_day]-year_file$Tmin[c_day])*
        sin((pi*(hour-year_file$Sunrise[c_day])/
               (year_file$Daylength[c_day]+4)))
      
      year_file[c_eve,colnum+hour]<-
        year_file$Tsunset[c_eve]- #temp at sunset
        ((year_file$Tsunset[c_eve]-year_file$next_min[c_eve])/
           log(24-(year_file$Sunset[c_eve]-year_file$next_Sunrise[c_eve])+1)*
           log(hour-year_file$Sunset[c_eve]+1))
      
    }
    colnames(year_file)[(ncol(year_file)-23):(ncol(year_file))]<-c(paste("Hour_",0:23,sep=""))
    if (!keep_sunrise_sunset)
      year_file<-year_file[,c(preserve_columns,paste("Hour_",0:23,sep=""))]
    if (keep_sunrise_sunset)
      year_file<-year_file[,c(preserve_columns,"Sunrise","Sunset","Daylength",paste("Hour_",0:23,sep=""))]
    year_file[1,(ncol(year_file)-23):(ncol(year_file))][which(is.na(year_file[1,(ncol(year_file)-23):(ncol(year_file))]))]<-year_file[1,"Tmin"]
    year_file[nrow(year_file),(ncol(year_file)-23):(ncol(year_file))][which(is.na(year_file[nrow(year_file),(ncol(year_file)-23):(ncol(year_file))]))]<-year_file[nrow(year_file),"Tmin"]
    
    return(year_file)
  }
