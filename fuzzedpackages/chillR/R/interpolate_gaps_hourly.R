#' Interpolate gaps in hourly temperature records
#' 
#' Using idealized temperature curves for guidance, this function interpolated hourly
#' temperature data.
#' 
#' Many agroclimatic metrics are calculated from hourly temperature data. chillR
#' provides functions for generating hourly data from daily records, which are
#' often available. Small gaps in such daily records can easily be closed
#' through linear interpolation, with relatively small errors, so that complete
#' hourly records can be generated. However, many sites have recorded actual
#' hourly temperatures, which allow much more accurate site-specific assessments.
#' Such records quite often have gaps, which need to be closed before
#' calculating most agroclimatic metrics (such as Chill Portions). Linear
#' interpolation is not a good option for this, because daily temperature curves
#' are not linear. Moreover, when gaps exceed a certain number of hours,
#' important featured would be missed (e.g. interpolating between temperatures
#' at 8 pm and 8 am may miss all the cool hours of the day, which would greatly
#' distort chill estimates).
#' 
#' This function solves this problem by using an idealized daily temperature
#' curve as guide to the interpolation of hourly temperature data.
#' 
#' These are the steps:
#' 1) produce an idealized temperature curve for the site (which requires
#' site latitude as an input), assuming minimum and maximum temperatures of
#' 0 and 1 degrees C, respectively. The calculations are based on equations published
#' by Spencer (1971), Almorox et al. (2005) and Linvill (1990, though I
#' modified these slightly to produce a smooth curve). This curve describes
#' the expected relationship of the temperature for the respective hour with
#' minimum and maximum temperatures of the same, previous or next day
#' (depending on the time of day), according to idealized temperature curve.
#' At this point, however, these daily minimum or maximum temperatures aren't
#' known yet.
#' 
#' 2) determine minimum and maximum temperatures for each day. For each minimum
#' and maximum daily temperature, the expected relationships between hourly
#' temperatures and daily extremes determined in step 1, combined with the
#' hourly temperatures that were observed can be interpreted as an
#' overdetermined set of equations that define these temperatures. Since few
#' days will follow the ideal curve precisely, and there are usually more than
#' two equations that define the same daily temperature extreme value, these
#' equations can only be solved numerically. This is implemented with the
#' qr.solve function, which can provide estimates of the minimum and maximum
#' temperatures for all days from the available hourly records.
#' 
#' 3) interpolate gaps in the record of estimated daily temperature extremes.
#' There can be days, when the number of recorded hourly temperatures isn't
#' sufficient for inferring daily minimum or maximum temperatures. The resulting
#' gaps are closed by linear interpolation (this may produce poor results if
#' gaps are really large, but this isn't currently addressed).
#' 
#' 4) compute an idealized daily temperature
#' curve for all days, based on estimated daily temperature extremes (using
#' the make_hourly_temperatures function).
#' 
#' 5) calculate deviation of recorded temperatures from idealized curve.
#' 
#' 6) linearly interpolate deviation values using the interpolate_gaps
#' function.
#' 
#' 7) add interpolated deviation values to idealized temperature curve.
#' 
#' 
#' @param hourtemps data.frame containing hourly temperatures. This has to
#' contain columns c("Year","Month","Day","Hour","Temp").
#' @param latitude the geographic latitude (in decimal degrees) of the location
#' of interest
#' @return data frame containing interpolated temperatures for all hours within
#' the interval defined by the first and last day of the hourtemps input.
#' @param daily_temps list of (chillR compliant) daily temperature data sets
#' for patching gaps in the record.
#' @param interpolate_remaining boolean parameter indicating whether gaps remaining
#' after the daily record has been patched (or after solving temperature equations,
#' if (daily_temps==NULL)) should be linearly interpolated.
#' @param return_extremes boolean parameters indicating whether daily minimum
#' and maximum temperatures used for the interpolation should be part of the
#' output table. Defaults to FALSE.
#' @param minimum_values_for_solving integer specifying the minimum number of hourly
#' temperature values that must be available for the solving function to be
#' applied. Must be greater than 1 (otherwise you get an error). Since according to
#' the idealized temperature curves used here, a given daily extreme temperature
#' is related to hourly temperatures of about a 12-hour period, values above 12
#' are not useful. Note that relatively large numbers for this parameter raise the
#' reliability of the interpolated values, but they restrict the number of missing
#' values in a day, for which the procedure produces results.
#' @param runn_mean_test_length integer specifying the length of the period, for
#' which a running mean test for is applied to daily records after the solving
#' procedure. This aims to remove spurious values that can sometimes arise during
#' solving. This test checks for all daily minimum and maximum temperature values,
#' if they differ from the mean of the surrounding values by more than
#' runn_mean_test_diff. If this is the case, they are set to NA, and have to be
#' filled by other means (from proxy data or by interpolation). Defaults to 5,
#' which means each value is compared to the mean of the 2 previous and 2
#' following days.
#' @param runn_mean_test_diff integer specifying the maximum tolerable difference
#' between solved daily extreme temperature values and the mean for the
#' surrounding days. See description of runn_mean_test_length for more details.
#' Defaults to 5.
#' @param daily_patch_max_mean_bias maximum acceptable mean difference between
#' the daily extreme temperatures of daily temperature records used as proxy and
#' daily extreme temperatures in the dataset that is to be interpolated. If the
#' bias between stations is greater than this, the station is not considered
#' a useful proxy and not used for filling gaps.
#' @param daily_patch_max_stdev_bias maximum acceptable standard deviation of
#' the difference between the daily extreme temperatures of daily temperature
#' records used as proxy and daily extreme temperatures in the dataset that is
#' to be interpolated. If the bias between stations is greater than this,
#' the station is not considered a useful proxy and not used for filling gaps.  
#' @author Eike Luedeling
#' @references 
#' Linvill DE, 1990. Calculating chilling hours and chill units from daily
#' maximum and minimum temperature observations. HortScience 25(1), 14-16.
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
#'
#' Winters_gaps<-make_JDay(Winters_hours_gaps[1:2000,])
#' colnames(Winters_gaps)[5:6]<-c("Temp","original_Temp")
#' interp<-interpolate_gaps_hourly(hourtemps=Winters_gaps,latitude=38.5)
#' 
#' #plot results: interpolated temperatures are shown in red, measured temperatures in black.
#' plot(interp$weather$Temp[1:120]~c(interp$weather$JDay[1:120]+
#'    interp$weather$Hour[1:120]/24),type="l",
#'    col="RED",lwd=2,xlab="JDay",ylab="Temperature")
#' lines(interp$weather$Temp_measured[1:120]~c(interp$weather$JDay[1:120]+
#'    interp$weather$Hour[1:120]/24),lwd=2)
#' 
#' @export interpolate_gaps_hourly
interpolate_gaps_hourly<-function(hourtemps,latitude=50,daily_temps=NULL,
                                  interpolate_remaining=TRUE,return_extremes=FALSE,
                                  minimum_values_for_solving=5,runn_mean_test_length=5,
                                  runn_mean_test_diff=5,
                                  daily_patch_max_mean_bias=NA,
                                  daily_patch_max_stdev_bias=NA)
  
{
  if((length(names(hourtemps))==2) & ("hourtemps" %in% names(hourtemps))) {QC<-hourtemps$QC; hourtemps<-hourtemps$hourtemps} else QC<-NULL
  
  hs<-hourtemps
  
  if(!("Year" %in% names(hs)&
       "Month" %in% names(hs)&
       "Day" %in% names(hs)&
       "Hour" %in% names(hs)&
       "Temp" %in% names(hs))) stop(paste("File is missing one of the following",
                                          "columns: Year, Month, Day, Hour, Temp"))
  
  dailyextremes<-make_all_day_table(rbind(hs[1,],hs[nrow(hs),])[,c("Year","Month","Day")],
                                    no_variable_check = TRUE)
  #hs,"day",input_timestep="hour")[,c("Year","Month","Day")]
  dailyextremes[,"Tmin"]<-0
  dailyextremes[,"Tmax"]<-1
  ideal_temps<-stack_hourly_temps(dailyextremes,latitude=latitude,
                                  keep_sunrise_sunset = TRUE)$hourtemps
  ideal_temps[,"YEARMODAHO"]<-ideal_temps[,"Year"]*1000000+
    ideal_temps[,"Month"]*10000+ideal_temps[,"Day"]*100+ideal_temps[,"Hour"]
  ideal_temps<-ideal_temps[,c("YEARMODAHO","Sunrise","Sunset","Daylength","Temp")]
  colnames(ideal_temps)[which(colnames(ideal_temps)=="Temp")]<-"ideal_temp"
  
  hs<-make_all_day_table(hs,"hour")
  hs[,"YEARMODA"]<-hs[,"Year"]*10000+hs[,"Month"]*100+hs[,"Day"]
  hs[,"YEARMODAHO"]<-hs[,"YEARMODA"]*100+hs[,"Hour"]
  
  hs<-merge(hs,ideal_temps,by.x="YEARMODAHO",by.y="YEARMODAHO")
  # hs[,"TminDay"]<-NA
  #  hs[,"TmaxDay"]<-NA
  
  hs[,"BeforeDay"]<-Date2YEARMODA(YEARMODA2Date(hs[,"YEARMODA"])-86400)
  hs[,"AfterDay"]<-Date2YEARMODA(YEARMODA2Date(hs[,"YEARMODA"])+86400)
  
  hs[which(hs[,"Hour"]<hs[,"Sunrise"]),"TminDay"]<-hs[which(hs[,"Hour"]<hs[,"Sunrise"]),"YEARMODA"]
  hs[which(hs[,"Hour"]<hs[,"Sunrise"]),"TmaxDay"]<-hs[which(hs[,"Hour"]<hs[,"Sunrise"]),"BeforeDay"]
  hs[which(hs[,"Hour"]>hs[,"Sunset"]),"TminDay"]<-hs[which(hs[,"Hour"]>hs[,"Sunset"]),"AfterDay"]
  hs[which(hs[,"Hour"]>hs[,"Sunset"]),"TmaxDay"]<-hs[which(hs[,"Hour"]>hs[,"Sunset"]),"YEARMODA"]
  hs[which(is.na(hs[,"TmaxDay"])),"TmaxDay"]<-hs[which(is.na(hs[,"TmaxDay"])),"YEARMODA"]
  hs[which(is.na(hs[,"TminDay"])),"TminDay"]<-hs[which(is.na(hs[,"TminDay"])),"YEARMODA"]  
  
  
  for(Tmd in hs[,"TminDay"])
  {d<-which(hs[,"TminDay"]==Tmd)
  Tminsolve<-hs[d,]
  A<-matrix(c(1-Tminsolve[,"ideal_temp"],
              (Tminsolve$TmaxDay<Tmd)*Tminsolve[,"ideal_temp"],
              (Tminsolve$TmaxDay>=Tmd)*Tminsolve[,"ideal_temp"]),nrow(Tminsolve))
  b<-hs[d,"Temp"]
  if(!length(which(!is.na(b)))<minimum_values_for_solving)
  {A1<-A[which(!is.na(b)),]
  A1<-A1[,which(colSums(A1)>0)]
  b1<-b[which(!is.na(b))]
  sol<-qr.solve(A1, b1)
  hs[d,"Tmin_solved"]<-sol[1]}
  }
  
  for(Tmd in hs[,"TmaxDay"])
  {d<-which(hs[,"TmaxDay"]==Tmd)
  Tmaxsolve<-hs[d,]
  A<-matrix(c(Tmaxsolve[,"ideal_temp"],
              (Tmaxsolve$TminDay==Tmd)*(1-Tmaxsolve[,"ideal_temp"]),
              (Tmaxsolve$TminDay>Tmd)*(1-Tmaxsolve[,"ideal_temp"])),nrow(Tmaxsolve))
  b<-hs[d,"Temp"]
  if(!length(which(!is.na(b)))<minimum_values_for_solving)
  {A1<-A[which(!is.na(b)),]
  A1<-A1[,which(colSums(A1)>0)]
  b1<-b[which(!is.na(b))]
  sol<-qr.solve(A1, b1)
  hs[d,"Tmax_solved"]<-sol[1]}
  }
  
  solved<-hs
  hs<-solved
  
  
  # hs[which(hs[,"Tmin_solved"]>hs[,"Tmax_solved"]),c("Tmin_solved","Tmax_solved")]<-NA
  
  Tmin_hs<-hs[,c("TminDay","Tmin_solved")]
  colnames(Tmin_hs)<-c("YEARMODA","Tmin")
  allmins<-make_all_day_table(Tmin_hs,no_variable_check = TRUE)
  runn_diff<-abs(allmins[,"Tmin"]-runn_mean(allmins[,"Tmin"],runn_mean_test_length,na.rm=TRUE,exclude_central_value = TRUE))
  allmins[which(runn_diff>runn_mean_test_diff),"Tmin"]<-NA
  #allmins[,"Tmin_source"]<-NA
  
  Tmax_hs<-hs[,c("TmaxDay","Tmax_solved")]
  colnames(Tmax_hs)<-c("YEARMODA","Tmax")
  allmaxs<-make_all_day_table(Tmax_hs,no_variable_check = TRUE)
  runn_diff<-abs(allmaxs[,"Tmax"]-runn_mean(allmaxs[,"Tmax"],runn_mean_test_length,na.rm=TRUE,exclude_central_value = TRUE))
  allmaxs[which(runn_diff>runn_mean_test_diff),"Tmax"]<-NA
  allmaxs[,"Tmax_source"]<-NA
  
  patch_report<-data.frame(Var=c("Tmin","Tmax"),Proxy="solved",mean_bias=NA,stdev_bias=NA,filled=c(length(which(!is.na(allmins$Tmin))),
                                                                                                   length(which(!is.na(allmaxs$Tmax)))),
                           gaps_remain=c(length(which(is.na(allmins$Tmin))),
                                         length(which(is.na(allmaxs$Tmax)))),
                           stringsAsFactors=FALSE)
  
  if(!is.null(daily_temps))
  {
    min_patched<-patch_daily_temperatures(weather=allmins,patch_weather=daily_temps,vars="Tmin",
                                          max_mean_bias = daily_patch_max_mean_bias, max_stdev_bias=daily_patch_max_stdev_bias)
    patch_report_min<-min_patched$statistics
    min_patched<-min_patched$weather
    colnames(min_patched)[which(colnames(min_patched)=="YEARMODA")]<-"TminDay"
    if(!"Tmin_source" %in% names(min_patched))
      min_patched[,"Tmin_source"]<-NA
    min_patched$Tmin_source[which(!is.na(min_patched$Tmin)&is.na(min_patched$Tmin_source))]<-"solved"
    
    max_patched<-patch_daily_temperatures(weather=allmaxs,patch_weather=daily_temps,vars="Tmax",
                                          max_mean_bias = daily_patch_max_mean_bias, max_stdev_bias=daily_patch_max_stdev_bias)
    patch_report_max<-max_patched$statistics
    max_patched<-max_patched$weather
    colnames(max_patched)[which(colnames(max_patched)=="YEARMODA")]<-"TmaxDay"
    if(!"Tmax_source" %in% names(max_patched))
      min_patched[,"Tmax_source"]<-NA
    max_patched$Tmax_source[which(!is.na(max_patched$Tmax)&is.na(max_patched$Tmax_source))]<-"solved"
    
    allmins<-min_patched
    allmaxs<-max_patched
    
    rep_mins<-as.data.frame(t(sapply(patch_report_min,FUN=function(x) x[1,])))
    rep_maxs<-as.data.frame(t(sapply(patch_report_max,FUN=function(x) x[1,])))
    
    if(is.null(names(patch_report_min)))
      names(patch_report_min)<-paste("station_",1:length(patch_report_min),sep="")
    if(is.null(names(patch_report_max)))
      names(patch_report_max)<-paste("station_",1:length(patch_report_max),sep="")
    
    patch_report<-rbind(patch_report[1,],
                        data.frame(Var=rep("Tmin",length(daily_temps)),Proxy=names(patch_report_min),
                                   mean_bias=unlist(rep_mins$mean_bias),
                                   stdev_bias=unlist(rep_mins$stdev_bias),
                                   filled=unlist(rep_mins$filled),
                                   gaps_remain=unlist(rep_mins$gaps_remain),
                                   stringsAsFactors=FALSE),
                        patch_report[2,],
                        data.frame(Var=rep("Tmax",length(daily_temps)),Proxy=names(patch_report_max),
                                   mean_bias=unlist(rep_maxs$mean_bias),
                                   stdev_bias=unlist(rep_maxs$stdev_bias),
                                   filled=unlist(rep_maxs$filled),
                                   gaps_remain=unlist(rep_maxs$gaps_remain),
                                   stringsAsFactors=FALSE)
    )
  }
  
  if(is.null(daily_temps))
  {if(!"Tmin_source" %in% colnames(allmins)) allmins[,"Tmin_source"]<-NA
  allmins$Tmin_source[which(!is.na(allmins$Tmin))]<-"solved"
  colnames(allmins)[which(colnames(allmins)=="YEARMODA")]<-"TminDay"
  if(!"Tmax_source" %in% colnames(allmaxs)) allmaxs[,"Tmax_source"]<-NA
  allmaxs$Tmax_source[which(!is.na(allmaxs$Tmax))]<-"solved"
  colnames(allmaxs)[which(colnames(allmaxs)=="YEARMODA")]<-"TmaxDay"}  
  
  
  
  if(interpolate_remaining)
  {allmins$Tmin<-interpolate_gaps(allmins$Tmin)$interp
  allmins$Tmin_source[which(!is.na(allmins$Tmin)&
                              is.na(allmins$Tmin_source))]<-"interpolated"
  patch_report<-rbind(patch_report[which(patch_report$Var=="Tmin"),],c(Var="Tmin",Proxy="interpolated",mean_bias=NA,stdev_bias=NA,filled=length(which(allmins$Tmin_source=="interpolated")),gaps_remain=0),
                      patch_report[which(patch_report$Var=="Tmax"),])
  allmaxs$Tmax<-interpolate_gaps(allmaxs$Tmax)$interp
  allmaxs$Tmax_source[which(!is.na(allmaxs$Tmax)&
                              is.na(allmaxs$Tmax_source))]<-"interpolated"
  patch_report<-rbind(patch_report,c(Var="Tmax",Proxy="interpolated",mean_bias=NA,stdev_bias=NA,filled=length(which(allmaxs$Tmax_source=="interpolated")),gaps_remain=0))
  }
  
  patch_report[which(is.na(patch_report$filled)),"filled"]<-0
  for(ll in 2:nrow(patch_report))
    if(patch_report[ll-1,"gaps_remain"]==0&is.na(patch_report[ll,"gaps_remain"])) patch_report[ll,"gaps_remain"]<-0
  
  row.names(patch_report)<-NULL
  
  hs_all<-merge(hs[,which(!(colnames(hs) %in% c("Tmin","Tmax")))],allmins[,c("TminDay","Tmin","Tmin_source")],by="TminDay")  
  hs_all<-merge(hs_all,allmaxs[,c("TmaxDay","Tmax","Tmax_source")],by="TmaxDay")
  
  hs_all<-hs_all[,which(!colnames(hs_all) %in% c("Tmin_solved","Tmax_solved",
                                                 "BeforeDay","AfterDay"))]
  
  
  hs_all[,"Temp_idealized"]<-hs_all[,"Tmin"]+
    hs_all[,"ideal_temp"]*(hs_all[,"Tmax"]-hs_all[,"Tmin"])
  
  hs_all[,"deviation"]<-hs_all[,"Temp"]-hs_all[,"Temp_idealized"]
  hs_all[,"deviation"]<-interpolate_gaps(hs_all[,"deviation"])$interp
  hs_all[,"Temp_interp"]<-hs_all[,"Temp_idealized"]+hs_all[,"deviation"]
  hs_all<-hs_all[,which(!colnames(hs_all) %in% c("deviation","Temp_idealized","Tmax_solved",
                                                 "Tmin_solved","TmaxDay","TminDay","ideal_temp",
                                                 "Sunrise","Sunset","Daylength","BeforeDay","AfterDay"))]
  colnames(hs_all)[which(colnames(hs_all)=="Temp")]<-"Temp_measured"
  colnames(hs_all)[which(colnames(hs_all)=="Temp_interp")]<-"Temp"
  hs_all$Tmin_source[which(!is.na(hs_all$Temp_measured))]<-NA
  hs_all$Tmax_source[which(!is.na(hs_all$Temp_measured))]<-NA
  if(!return_extremes) hs_all<-hs_all[,which(!colnames(hs_all) %in% c("Tmin","Tmax"))]
  hs_all<-hs_all[,which(!colnames(hs_all) %in% c("YEARMODA","YEARMODAHO","DATE"))]
  
  return(list(weather=hs_all,daily_patch_report=patch_report))
}

