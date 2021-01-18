#' Check whether a Julian date is before or after another one
#'
#' For two Julian dates, this function checks whether the first date
#' is earlier than the second date within a user-defined phenological season.
#' This is particularly useful for seasons that start in one year and end
#' in the next, because simple > or < operations can produce wrong results
#' then.
#'
#' @param check_date integer ranging from 1 to 366, indicating a Julian date.
#' This is the date for which to check whether it is before the reference date.
#' If this is a vector, all elements are checked against the reference date.
#' @param ref_date integer ranging from 1 to 366, indicating a Julian date.
#' This is the reference date.
#' @param season integer vector of length 2, specifying the beginning and end
#' of the phenology season, respectivcely.
#' @return Boolean result (TRUE/FALSE) of the comparison.
#' @author Eike Luedeling
#' @keywords phenology season Julian date
#' @examples
#'
#' JDay_earlier(check_date=10,ref_date=365,season=c(305,59))
#'
#' @export JDay_earlier
JDay_earlier<-function(check_date,ref_date,season=c(1,366))
{
  if(!season[1] %in% 1:366) stop("Start of season is not a Julian date")
  if(!season[2] %in% 1:366) stop("End of season is not a Julian date")
  if(sum(check_date %in% 1:366)==0) stop("check_date is not a Julian date")
  if(!ref_date %in% 1:366) stop("ref_date is not a Julian date")
  
  if(season[2]>season[1]) alldays<-season[1]:season[2]
  if(season[2]<season[1]) alldays<-c(season[1]:366,1:season[2])
  
  #if(sum(check_date %in% alldays)==0)
  #  return(NA)
  if(!ref_date %in% alldays)    
    return(NA)
  
  check_dates<-
    unlist(lapply(check_date,
                  function(x) {xx<-which(x==alldays)
                  if(length(xx)==0)
                    return(NA) else return(xx)}))
  
  return(check_dates<which(ref_date==alldays))
}  
#' Check whether a Julian date is after another one
#'
#' For two Julian dates, this function checks whether the first date
#' is later than the second date within a user-defined phenological season.
#' This is particularly useful for seasons that start in one year and end
#' in the next, because simple > or < operations can produce wrong results
#' then.
#'
#' @param check_date integer ranging from 1 to 366, indicating a Julian date.
#' This is the date for which to check whether it is after the reference date.
#' If this is a vector, all elements are checked against the reference date.
#' @param ref_date integer ranging from 1 to 366, indicating a Julian date.
#' This is the reference date.
#' @param season integer vector of length 2, specifying the beginning and end
#' of the phenology season, respectivcely.
#' @return Boolean result (TRUE/FALSE) of the comparison.
#' @author Eike Luedeling
#' @keywords phenology season Julian date
#' @examples
#'
#' JDay_later(check_date=10,ref_date=365,season=c(305,59))
#'
#' @export JDay_later
JDay_later<-function(check_date,ref_date,season=c(1,366))
{
  if(!season[1] %in% 1:366) stop("Start of season is not a Julian date")
  if(!season[2] %in% 1:366) stop("End of season is not a Julian date")
  if(sum(check_date %in% 1:366)==0) stop("check_date is not a Julian date")
  if(!ref_date %in% 1:366) stop("ref_date is not a Julian date")
  
  if(season[2]>season[1]) alldays<-season[1]:season[2]
  if(season[2]<season[1]) alldays<-c(season[1]:366,1:season[2])
  
  #if(sum(check_date %in% alldays)==0)
  #  return(NA)
  if(!ref_date %in% alldays)    
    return(NA)
  
  check_dates<-
    unlist(lapply(check_date,
                  function(x) {xx<-which(x==alldays)
                  if(length(xx)==0)
                    return(NA) else return(xx)}))
  
  return(check_dates>which(ref_date==alldays))
}

#' Count days between two Julian dates
#'
#' This function counts the days between two Julian dates, taking into
#' account whether the season extends past the end of a calender year and
#' whether the count is to be done for a leap year.
#'
#' @param start_date integer ranging from 1 to 366, indicating a Julian date.
#' This is the start date of the interval of interest.
#' @param end_date integer ranging from 1 to 366, indicating a Julian date.
#' This is the end date of the interval of interest.
#' @param season integer vector of length 2, specifying the beginning and end
#' of the phenology season, respectivcely. If this is not specified, the
#' start_date and end_date are used to define the season.
#' @param leap_year either a Boolean parameter indicating whether the count
#' should be done for a leap year, or an integer specyfing the year, for
#' which the calculation is to be done. The function then determines
#' automatically, whether this is a leap year.
#' @return Boolean result (TRUE/FALSE) of the comparison.
#' @author Eike Luedeling
#' @keywords phenology season Julian date
#' @examples
#'
#' JDay_count(start_date=320,end_date=20,season=c(305,59),leap_year=2004)
#'
#' @export JDay_count
JDay_count<-function(start_date,end_date,season=NA,leap_year=FALSE)
{
  if(is.na(season[1])) season<-c(start_date,end_date)
  if(!season[1] %in% 1:366) stop("Start of season is not a Julian date")
  if(!season[2] %in% 1:366) stop("End of season is not a Julian date")
  if(!start_date %in% 1:366) stop("start_date is not a Julian date")
  if(!end_date %in% 1:366) stop("end_date is not a Julian date")
  
  if(is.numeric(leap_year))
    {if(leap_year==round(leap_year))
      leapyear<-leap_year(leap_year)} else leapyear<-leap_year
  
  
  if(season[2]>season[1]) alldays<-season[1]:season[2]
  if(season[2]<season[1]) alldays<-c(season[1]:366,1:season[2])
  
  if(!leapyear) alldays<-alldays[which(!alldays==366)]
  
  start_dates<-
    unlist(lapply(start_date,
                  function(x) {xx<-which(x==alldays)
                  if(length(xx)==0)
                    return(NA) else return(xx)}))
  
  end_dates<-
    unlist(lapply(end_date,
                  function(x) {xx<-which(x==alldays)
                  if(length(xx)==0)
                    return(NA) else return(xx)}))
  
  return(end_dates-start_dates)
}

