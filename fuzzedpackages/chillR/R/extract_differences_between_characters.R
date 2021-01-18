#' Identify shared leading or trailing character strings
#' 
#' For a vector of character strings, identify elements between shared leading and/or trailing substrings,
#' e.g. for a vector such as c("XXX01YYY",XXX02YYY") extract the numbers.
#' 
#' @param strings vector of character strings for elements to be extracted from.
#'  
#' @return vector of strings similar to the input vector but without shared leading and trailing characters.
#' 
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#'   extract_differences_between_characters(c("Temp_01","Temp_02","Temp_03"))
#'   extract_differences_between_characters(c("Temp_01_Tmin","Temp_02_Tmin","Temp_03_Tmin"))
#'   extract_differences_between_characters(c("a","b"))                                           
#'  
#' @export extract_differences_between_characters
extract_differences_between_characters<-function(strings)
{
  if(is.null(strings)) return(NA)
  if(length(strings)==1) return(strings)
  leader<-identify_common_string(strings,leading=TRUE)
  trailer<-identify_common_string(strings,leading=FALSE)
  #doesn't work at the moment if leader and trailer are the same (or if leading char occurs twice)
  #rather than extracting only the last element below, piece together all except the last, adding the leader between them - complicated.
  if(!is.na(leader))
    no_leader<-sapply(as.character(strings),function(x) substr(x,nchar(leader)+1,nchar(x))) else
       no_leader<-strings
  if(test_if_equal(no_leader)) return(NA)
  if(!is.na(trailer))
    no_leader_no_trailer<-sapply(as.character(no_leader),function(x) substr(x,1,nchar(x)-nchar(trailer))) else
      no_leader_no_trailer<-no_leader      

  return(as.character(no_leader_no_trailer))
}
