#' Quality filter for temperature records
#' 
#' This function attempts to remove erroneous temperature readings. This is
#' tricky because of the wide range of errors that can occur, so this isn't
#' necessarily sufficient for problems of particular records.
#' 
#' 
#' @param temp_file file containing temperature data. Should have columns
#' c("Year","Month","Day","Temp" - and "Hour" for hourly data).
#' @param remove_value numeric value indicating 'no data'.
#' @param running_mean_filter deviation from a running mean over all temperature
#' data that identifies a value as an erroneous outlier.
#' @param running_mean_length number of records to be included in a running mean.
#' @param min_extreme lowest plausible temperature on the record. All lower ones are removed.
#' @param max_extreme highest plausible temperature on the record. All higher ones are removed.
#' @param max_missing_in_window maximum share of values (0..1) in a running window of size
#' missing_window_size around each value that can be missing. If this is exceeded, the value
#' is removed.
#' @param missing_window_size size of the window used for checking for missing values.
#' @return filtered temperature dataset, from which records identified as erroneous were removed.
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#' 
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2009),])
#' 
#' hourtemps<-stack_hourly_temps(weather, latitude=50.4)
#' 
#' filtered<-filter_temperatures(hourtemps$hourtemps,remove_value=-99,
#'   running_mean_filter=3)
#' 
#' @export filter_temperatures
filter_temperatures <-
function(temp_file,remove_value=NA,running_mean_filter=NA,running_mean_length=3,
         min_extreme=NA,max_extreme=NA,max_missing_in_window=1,missing_window_size=9)
{
  
  runn_mean<-function(vec,runn_mean,na.rm=FALSE,exclude_central_value=FALSE,FUN=mean)
  {
    if(identical(FUN,mean)) if(na.rm==TRUE) FUN<-function(x) mean(x,na.rm=TRUE) else FUN<-function(x) mean(x,na.rm=FALSE)
    
    ww <- vec
  rr <- vec
  for (dd in 1:length(ww)) {
    if (dd < ceiling(runn_mean/2)) {
      if(!exclude_central_value) rr[dd] <- sapply(list(ww[1:(dd + floor(runn_mean/2))]),FUN)
      if(exclude_central_value) rr[dd] <- sapply(list(ww[(1:(dd + floor(runn_mean/2)))[which(!(1:(dd + floor(runn_mean/2)))==dd)]]),FUN)
    }
    if ((dd >= ceiling(runn_mean/2)) & (dd <= length(ww) - 
                                        ceiling(runn_mean/2))) {
      if(!exclude_central_value) rr[dd] <- sapply(list(ww[(dd - floor(runn_mean/2)):(dd + 
                                                                                floor(runn_mean/2))]),FUN)
      if(exclude_central_value) rr[dd] <- sapply(list(ww[((dd - floor(runn_mean/2)):(dd + 
                                                                                floor(runn_mean/2)))[
                                                                                  which(!((dd - floor(runn_mean/2)):(dd + 
                                                                                                                       floor(runn_mean/2)))==dd)]]),FUN)
    }
    if (dd > (length(ww) - ceiling(runn_mean/2))) {
      if(!exclude_central_value) rr[dd] <- sapply(list(ww[(dd - floor(runn_mean/2)):length(ww)]),FUN)
      if(exclude_central_value) rr[dd] <- sapply(list(ww[((dd - floor(runn_mean/2)):length(ww))[
        which(!((dd - floor(runn_mean/2)):length(ww))==dd)]]),FUN)
    }
  }
  return(rr)}
  
  
  
  
  
  if(!("Year" %in% colnames(temp_file)&
       "Month" %in% colnames(temp_file)&
       "Day" %in% colnames(temp_file)&
       "Temp" %in% colnames(temp_file)))
     stop("Error: required column missing: Year, Month, Day or Temp")
    
  if("Hour" %in% colnames(temp_file))
    alldays<-make_all_day_table(temp_file,timestep = "hour") else
      alldays<-make_all_day_table(temp_file)
    
    if("Minutes" %in% colnames(temp_file)) alldays[,"Minutes"]<-0
    
    if(!is.na(running_mean_filter))
    {runn_diff<-abs(alldays[,"Temp"]-runn_mean(alldays[,"Temp"],running_mean_length,na.rm=TRUE,exclude_central_value = TRUE))
     alldays[which(runn_diff>running_mean_filter),"Temp"]<-NA
    }
    
  if(!is.na(remove_value)) alldays[which(alldays[,"Temp"]==remove_value),"Temp"]<-NA
  if(!is.na(min_extreme)) alldays[which(alldays[,"Temp"]<min_extreme),"Temp"]<-NA
  if(!is.na(max_extreme)) alldays[which(alldays[,"Temp"]>max_extreme),"Temp"]<-NA
  
  #remove all values, for which more than 50% of the neighbors are missing in an 11-value window
  #(5 before, 5 after)
  if(max_missing_in_window<1)
    {missings<-runn_mean(alldays$Temp,missing_window_size,exclude_central_value=TRUE,FUN=function(x) length(which(is.na(x)))/length(x))
    alldays[which(missings>max_missing_in_window),"Temp"]<-NA}
  
  
  if(!"YEARMODA" %in% colnames(temp_file)) alldays<-alldays[,which(!colnames(alldays)=="YEARMODA")]
  if(!"YEARMODAHO" %in% colnames(temp_file)) alldays<-alldays[,which(!colnames(alldays)=="YEARMODAHO")]
  if(!"DATE" %in% colnames(temp_file)) alldays<-alldays[,which(!colnames(alldays)=="DATE")]
  
  
  
  return(alldays)
}
 