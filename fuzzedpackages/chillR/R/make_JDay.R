#' Make Julian Day in dataframe
#' 
#' This function produced Julian Dates (days of the year) from columns "Day",
#' "Month" and "Year" in a dataframe.
#' 
#' @param dateframe data.frame, which should contain date information specified
#' as columns "Day", "Month" and "Year"
#' @return Returns the same data.frame, but with column "JDay" added. This then
#' contains the Julian Dates.
#' 
#' @author Eike Luedeling
#' @references The chillR package:
#' 
#' @keywords utilities
#' @examples
#' 
#' 
#' dates<-data.frame(Year=c(1977,1980,2004,2011,2016),Month=c(11,8,3,12,8),Day=c(1,21,2,24,2)) 
#' make_JDay(dates)
#' 
#' 
#' @export make_JDay
make_JDay<-function(dateframe)
{
  if(sum(is.element(c("Day","Month","Year"),colnames(dateframe)))==3)
   dateframe[,"JDay"]<-strptime(paste(dateframe$Month, 
                 "/", dateframe$Day, "/", dateframe$Year, sep = ""), "%m/%d/%Y")$yday + 1
  else print("Table is missing at least one required column ('Day','Month' or 'Year')")
  return(dateframe)
}