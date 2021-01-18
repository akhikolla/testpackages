#' Compute what it takes to advance through development stages
#' 
#' Function to compute the thermal requirements of transitioning
#' through a series of developmental stages. 
#' 
#' 
#' @param observations data.frame containing observed developmental dates,
#' e.g. different stages of flower or leaf development. Should contain the
#' columns 'Stage' (containing the names of the development stages), 'Season'
#' (containing the 'development year' the observation belongs to, e.g. budbreak
#' for trees may be considered a stage of the 'dormancy year' that started
#' in the previous calendar year), 'Year' (the calendar year the observation
#' was made), 'JDay' (the Julian Date, a.k.a. day of the year, that the stage
#' was observed).
#' @param hourtemps a list of two elements, with element 'hourtemps' being a
#' dataframe of hourly temperatures (e.g. produced by stack_hourly_temps). This
#' data frame must have a column for Year, a column for JDay (Julian date, or
#' day of the year), a column for Hour and a column for Temp (hourly
#' temperature). The second (optional) element is QC, which is a data.frame
#' indicating completeness of the dataset. This is automatically produced by
#' stack_hourly_temps. This also works if only the 'hourtemps' dataframe is
#' passed to the function.
#' @param stages character vector containing the relevant development stages
#' in their order of occurrence.
#' @param models named list of models that should be applied to the hourly
#' temperature data. These should be functions that take as input a vector of
#' hourly temperatures. This defaults to list(Chill_Portions=Dynamic_Model,
#' GDH=GDH), models that are often used for describing chill and heat
#' accumulation in temperate fruit trees.
#' @param max_steps integer indicating the maximum number of stage steps
#' (i.e. transitions from one step to the next), for which thermal requirements
#' should be calculated. This defaults to length(stages), which is also the
#' maximum value. If only requirements between each stage and the following
#' stage are of interest, this should be set to 1.
#' @return data frame with rows for all transitions that occurred during the
#' observed records and the values of the metrics specified in 'models' that
#' accrued between the respective dates. Columns are c('Season','Stage',
#' 'to_Stage','stage_steps') and one column for each thermal metrics.
#' @author Eike Luedeling
#' @keywords stage transitions
#' @examples
#' hourtemps<-stack_hourly_temps(KA_weather)
#' observations<-data.frame(Stage=c("V1","V2","V3","V1","V2","V3","V1","V3"),
#'                         Season=c(2001,2001,2001,2002,2002,2002,2003,2003),
#'                         Year=c(2001,2001,2001,2002,2002,2002,2003,2003),
#'                         JDay=c(30,45,60,35,42,55,37,62))
#' stages<-c("V1","V2","V3") 
#' 
#' stage_transitions(observations,hourtemps,stages)
#' 
#' @export stage_transitions
stage_transitions <-
  function (observations,hourtemps,stages,
            models=list(Chill_Portions=Dynamic_Model,GDH=GDH),
            max_steps=length(stages))
  {
    
    if("hourtemps" %in% names(hourtemps))
      hourtemps<-hourtemps$hourtemps
    hourtemps[,"JDay"]<-make_JDay(hourtemps)$JDay
    seasons<-unique(observations$Season)
    trans<-expand.grid(to_Stage=stages,Stage=stages,Season=seasons)[,3:1]
    stage_steps<-sapply(trans$to_Stage,function(x) which(x==stages))-
      sapply(trans$Stage,function(x) which(x==stages))
    stage_steps[which(stage_steps<1)]<-stage_steps[which(stage_steps<1)]+length(stages)
    trans[,"stage_steps"]<-stage_steps
    trans<-trans[which(trans$stage_steps<=max_steps),]
    trans[,names(models)]<-NA
    
    for(tr in 1:nrow(trans))
    {start_obs=observations[which(observations$Season==trans$Season[tr]&
                                    as.character(observations$Stage)==
                                    as.character(trans$Stage[tr])),]
    if(which(trans$to_Stage[tr]==stages)<=which(trans$Stage[tr]==stages))
      to_season<-trans$Season[tr]+1 else to_season<-trans$Season[tr]
      end_obs=observations[which(observations$Season==to_season&
                                   as.character(observations$Stage)==
                                   as.character(trans$to_Stage[tr])),]
      if((nrow(start_obs)==1)&(nrow(end_obs)==1))
      {
        start_weather<-hourtemps[which(hourtemps$Year==start_obs$Year&
                                         hourtemps$JDay==start_obs$JDay),]
        end_weather<-hourtemps[which(hourtemps$Year==end_obs$Year&
                                       hourtemps$JDay==end_obs$JDay),]
        if(nrow(start_weather)>0&nrow(end_weather)>0)
        {tempR<-tempResponse(hourtemps=
                               hourtemps[which(hourtemps$Year==start_obs$Year&hourtemps$JDay==start_obs$JDay)[1]:
                                           which(hourtemps$Year==end_obs$Year&hourtemps$JDay==end_obs$JDay)[length(
                                             which(hourtemps$Year==end_obs$Year&hourtemps$JDay==end_obs$JDay)
                                           )],],
                             Start_JDay = start_obs$JDay,End_JDay = end_obs$JDay,models=models,
                             whole_record=TRUE)
        trans[tr,names(models)]<-tempR
        }
      }
    }   
    return(trans)
}
