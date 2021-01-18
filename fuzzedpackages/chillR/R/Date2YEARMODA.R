#' Date to YEARMODA conversion
#' 
#' Converts R dates to YEARMODA format
#' 
#' Converts R date to YEARMODA
#' 
#' @param Date Date in R date format
#' @param hours boolean variable indicating whether YEARMODAHO should be
#' calculated (YEARMODA + hours)
#' @return YEARMODA object (e.g. 20111224 for 24th December 2011)
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#' 
#' Date2YEARMODA(YEARMODA2Date(20001205))
#' Date2YEARMODA(YEARMODA2Date(19901003))
#' 
#'  
#' @export Date2YEARMODA
Date2YEARMODA<-function(Date,hours=FALSE)
{
  if(!hours) YEARMODA<-as.numeric(format(as.Date(Date),"%Y"))*10000+
                       as.numeric(format(as.Date(Date),"%m"))*100+
                       as.numeric(format(as.Date(Date),"%d"))
  if(hours) YEARMODA<-as.numeric(format(as.Date(Date),"%Y"))*1000000+
                      as.numeric(format(as.Date(Date),"%m"))*10000+
                      as.numeric(format(as.Date(Date),"%d"))*100+
                      as.numeric(format(Date,"%H"))
 return(YEARMODA) 
}
