#' Patch gaps in daily weather records
#' 
#' This function uses auxiliary data sources to fill gaps in daily weather data.
#' It can accommodate multiple sources of auxiliary information, which are used
#' in the user-specified sequence. There have to be some overlapping records for
#' this to work, because without bias correction, this procedure could produce
#' erroneous records. Bias correction is done by computing the mean difference
#' between main and auxiliary data for each variable and adjusting for it in
#' filling the gaps. You can specify a maximum mean bias and a maximum standard
#' deviation of the bias to exclude unsuitable records that aren't similar
#' enough to the original data. 
#' 
#' @param weather chillR-compatible weather record to be patched
#' @param patch_weather list of chillR-compatible weather records to be used for
#' patching holes in weather. They are used sequentially, until all have been used
#' or until there are no holes left.
#' @param vars vector of column names to be considered in patching. Defaults to
#' c("Tmin","Tmax"), the most common variables in chillR applications.
#' @param max_mean_bias maximum mean bias of auxiliary data compared to the original
#' dataset (applied to all variables in vars). If this threshold is exceeded, the
#' respective variable from that particular dataset will not be used. Defaults to NA,
#' meaning no records are excluded.
#' @param max_stdev_bias maximum standard deviation of the bias in the auxiliary
#' data compared to the original dataset (applied to all variables in vars). If this
#' threshold is exceeded, the respective variable from that particular dataset will not
#' be used. Defaults to NA, meaning no records are excluded.
#' @return list of two elements: weather (the patched weather record, with additional
#' columns specifying the data source for each value) and statistics (containing
#' data.frames for each element of patch_weather that indicate the mean bias, the
#' number of values that were filled from this source and the number of missing records
#' that remained after exhausting this auxiliary data source.)
#' @author Eike Luedeling
#' @keywords temperature gap-filling
#' @examples
#' 
#' gap_weather<-KA_weather[1:100,]
#' gap_weather[c(3,4,7:15,20,22:25,27:28,35:45,55,67,70:75,80:88,95:97),"Tmin"]<-NA
#' gap_weather[c(10:25,30,36:44,50,57,65,70:80,86,91:94),"Tmax"]<-NA
#' p1<-KA_weather[65:95,]
#' p1$Tmin<-p1$Tmin-2
#' p2<-KA_weather[c(15:40,60:80),]
#' p2$Tmax<-p2$Tmax+3
#' p3<-KA_weather[12:35,]
#' p3$Tmax<-p3$Tmax-2
#' p4<-KA_weather
#' p4$Tmax<-p4$Tmax+0.5
#' patch_weather<-list(stat1=p1,st2=p2,home=p3,last=p4)
#' 
#' patched<-patch_daily_temperatures(gap_weather,patch_weather,max_mean_bias=1)
#' 
#' 
#' @export patch_daily_temperatures
patch_daily_temperatures<-function(weather,patch_weather,vars=c("Tmin","Tmax"),max_mean_bias=NA,
                                   max_stdev_bias=NA)
{  
  weather<-make_all_day_table(weather,no_variable_check = TRUE)
  if(!"YEARMODA" %in% colnames(weather))
    weather[,"YEARMODA"]<-weather[,"Year"]*10000+weather[,"Month"]*100+weather[,"Day"]
  
  if(is.data.frame(patch_weather)) daily<-list(patch_weather) else daily<-patch_weather
  
  statistics<-list()
  for(i in 1:length(daily))
    statistics[[length(statistics)+1]]<-data.frame(mean_bias=rep(NA,length(vars)),stdev_bias=NA,
                                                   filled=NA,gaps_remain=NA,row.names=vars)
  if(!is.null(names(daily))) names(statistics)<-names(daily)
  
  gaps<-length(which(is.na(weather[,vars])))
  
  for(dt in 1:length(daily))
  {
    if(gaps==0) 
    {if(dt==1) aux_weather<-weather} else
    {
      auxiliary<-make_all_day_table(daily[[dt]],no_variable_check = TRUE)
      auxiliary[,"YEARMODA"]<-auxiliary[,"Year"]*10000+auxiliary[,"Month"]*100+auxiliary[,"Day"]
      auxiliary<-auxiliary[,c("YEARMODA",vars)]
      vars2<-paste(vars,"temp",sep="")
      colnames(auxiliary)[which(colnames(auxiliary) %in% vars)]<-vars2
      
      if(dt==1) aux_weather<-merge(weather,auxiliary,by.x="YEARMODA",by.y="YEARMODA",all.x=TRUE) else
        aux_weather<-merge(aux_weather,auxiliary,by.x="YEARMODA",by.y="YEARMODA",all.x=TRUE)
      
      for(v in 1:length(vars))
      {bias<-mean(aux_weather[,vars[v]]-aux_weather[,vars2[v]],na.rm=TRUE)
      if(is.na(bias)) bias<-NA
      stdev<-sd(aux_weather[,vars[v]]-aux_weather[,vars2[v]],na.rm=TRUE)
      dont_use<-FALSE
      statistics[[dt]][vars[v],"mean_bias"]<-round(bias,3)
      statistics[[dt]][vars[v],"stdev_bias"]<-round(stdev,3)
      if(is.na(bias)) dont_use<-TRUE
      if(!is.na(bias))
      {if(!is.na(max_mean_bias)) if(abs(bias)>max_mean_bias) dont_use<-TRUE
      if(!is.na(max_stdev_bias)) if(stdev>max_stdev_bias) dont_use<-TRUE
      if(dt>1) if(statistics[[dt-1]][vars[v],"gaps_remain"]==0) dont_use<-TRUE
      }
      if(dont_use)
      {statistics[[dt]][vars[v],"filled"]<-0
      
      }
      if(!dont_use)
      {if(!is.null(names(daily)))
        aux_weather[which(is.na(aux_weather[,vars[v]])&
                            !is.na(aux_weather[,vars2[v]])),paste(vars[v],"_source",sep="")]<-
          paste("daily_",names(daily)[dt],sep="")
      if(is.null(names(daily)))
        aux_weather[which(is.na(aux_weather[,vars[v]])&
                            !is.na(aux_weather[,vars2[v]])),paste(vars[v],"_source",sep="")]<-
          paste("daily_records",dt,sep="")
      
      statistics[[dt]][vars[v],"filled"]<-
        length(which(is.na(aux_weather[,vars[v]])&
                       !is.na(aux_weather[,vars2[v]])))
      
      aux_weather[which(is.na(aux_weather[,vars[v]])&
                          !is.na(aux_weather[,vars2[v]])),vars[v]]<-
        aux_weather[which(is.na(aux_weather[,vars[v]])&
                            !is.na(aux_weather[,vars2[v]])),vars2[v]]+bias
      } else
        warnings("no overlap")
      
      gaps<-length(which(is.na(aux_weather[,vars[v]])))
      statistics[[dt]][vars[v],"gaps_remain"]<-gaps
      }
      aux_weather<-aux_weather[,which(!colnames(aux_weather) %in% c(vars2))]
    }
  }
  
  return(list(weather=aux_weather,statistics=statistics))
  
}

#KA_weather_gap<-rbind(KA_weather,c(Year=2011,Month=3,Day=3,Tmax=26,Tmin=14)) 

#weather<-KA_weather_gap
#p1<-KA_weather[95:105,]
#p1$Tmin<-p1$Tmin-2
#p2<-KA_weather[c(85:90,110:115),]
#p2$Tmax<-p2$Tmax+3
#p3<-KA_weather[120:125,]
#p3$Tmax<-p3$Tmax-2
#p4<-KA_weather
#p4$Tmax<-p4$Tmax+0.5
#patch_weather<-list(stat1=p1,st2=p4,home=p3,last=p4)

#patched<-patch_daily_temperatures(KA_weather_gap,patch_weather)
