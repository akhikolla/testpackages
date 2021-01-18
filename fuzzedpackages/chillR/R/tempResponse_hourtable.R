#' Add metric accumulation to table of hourly temperatures
#' 
#' This function calculates cumulative values for temperature response metrics
#' for every hour of an hourly temperature record. The count is
#' restarted on a specified date each year. The function is a generalized
#' version of chilling_hourtable, which only worked with three predefined
#' chilling one predefined heat metrics.
#' 
#'  
#' @param hourtemps a dataframe of stacked hourly temperatures (e.g. produced
#' by stack_hourly_temps). This data frame must have a column for Year, a
#' column for JDay (Julian date, or day of the year), a column for Hour and a
#' column for Temp (hourly temperature).
#' @param Start_JDay the start date (in Julian date, or day of the year) of the
#' calculation for the four metrics. The count is restarted on this date every
#' year.
#' @param models named list of models that should be applied to the hourly
#' temperature data. These should be functions that take as input a vector of
#' hourly temperatures. This defaults to c(Chill_Portions = Dynamic_Model, GDH
#' = GDH_model), which refer to the Dynamic chill model and the Growing Degree
#' Hours model functions contained in chillR.
#' @return data frame consisting of all the columns of the THourly input data
#' frame, plus one additional column for each model, which contains the
#' cumulative number of model metrics since the last Start_JDay).
#' @note After doing extensive model comparisons, and reviewing a lot of
#' relevant literature, I do not recommend using the Chilling Hours or Utah
#' Models, especially in warm climates! The Dynamic Model (Chill Portions),
#' though far from perfect, seems much more reliable.
#' @author Eike Luedeling
#' @references Model references:
#' 
#' 
#' Dynamic Model:
#' 
#' Erez A, Fishman S, Linsley-Noakes GC, Allan P (1990) The dynamic model for
#' rest completion in peach buds. Acta Hortic 276, 165-174
#' 
#' Fishman S, Erez A, Couvillon GA (1987a) The temperature dependence of
#' dormancy breaking in plants - computer simulation of processes studied under
#' controlled temperatures. J Theor Biol 126(3), 309-321
#' 
#' Fishman S, Erez A, Couvillon GA (1987b) The temperature dependence of
#' dormancy breaking in plants - mathematical analysis of a two-step model
#' involving a cooperative transition. J Theor Biol 124(4), 473-483
#' 
#' Growing Degree Hours:
#' 
#' Anderson JL, Richardson EA, Kesner CD (1986) Validation of chill unit and
#' flower bud phenology models for 'Montmorency' sour cherry. Acta Hortic 184,
#' 71-78
#' 
#' Review on chilling models in a climate change context:
#' 
#' Luedeling E, 2012. Climate change impacts on winter chill for temperate
#' fruit and nut production: a review. Scientia Horticulturae 144, 218-229
#' 
#' @keywords chill and heat calculation
#' @examples
#' 
#' weather<-fix_weather(KA_weather[which(KA_weather$Year>2008),])
#' 
#' hourtemps<-stack_hourly_temps(weather,latitude=50.4)
#' 
#' cht<-chilling_hourtable(hourtemps,20)
#' 
#' @export tempResponse_hourtable
tempResponse_hourtable <-
function (hourtemps,Start_JDay,
          models=c(Chill_Portions = Dynamic_Model, GDH = GDH_model))
{
  
  if((length(names(hourtemps))==2) & ("hourtemps" %in% names(hourtemps)) & ("QC" %in% names(hourtemps))) 
  {hourtemps<-hourtemps$hourtemps
  QC<-hourtemps$QC}
  
  cols<-colnames(hourtemps)
  
  hourtemps<-hourtemps[which(!is.na(hourtemps[,"Temp"])),]
  
  for(m in 1:length(models))
    hourtemps[,names(models)[m]]<-do.call(models[[m]], list(hourtemps[,"Temp"]),
                                          quote = FALSE, envir = parent.frame())
  
  if(!is.na(Start_JDay))
  { if(!"JDay" %in% colnames(hourtemps))
    hourtemps<-make_JDay(hourtemps)
    if(!"Season" %in% colnames(hourtemps))
    {hourtemps[which(hourtemps$JDay>=Start_JDay),"Season"]<-
      hourtemps[which(hourtemps$JDay>=Start_JDay),"Year"]
    hourtemps[which(hourtemps$JDay<Start_JDay),"Season"]<-
      hourtemps[which(hourtemps$JDay<Start_JDay),"Year"]-1}
    
    sea<-unique(hourtemps$Season)
    
    for(m in names(models))
      for (s in sea)
        hourtemps[which(hourtemps$Season==s),m]<-
          hourtemps[which(hourtemps$Season==s),m]-
           hourtemps[which(hourtemps$Season==s&hourtemps$JDay==round(Start_JDay)),m][1]
  }
  
  return(hourtemps[,c(cols,names(models))])
  
}
