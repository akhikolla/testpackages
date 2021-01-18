#' Identify shared leading or trailing character strings
#' 
#' Compares all elements of a vector of numbers or character strings and returns TRUE
#' if they are all the same, FALSE otherwise.
#' 
#' @param strings vector of strings to be evaluated.
#' @param leading boolean variable indicating whether the function should look for common strings at the beginning
#' (leading==TRUE) or end (leading==FALSE) of the strings. Default is TRUE. 
#'  
#' @return if there is a leading (if leading==TRUE) or trailing (if leading==FALSE) string that all elements of
#' strings have in common, this string is returned; NA otherwise.
#' 
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#'   identify_common_string(c("Temp_01","Temp_02","Temp_03"))
#'   identify_common_string(c("Temp_01","Temp_02","Temp_03"),leading=FALSE)
#'   identify_common_string(c("file1.csv","file2.csv","file3.csv"),leading=FALSE)
#'     
#' @export identify_common_string
identify_common_string<-function(strings, leading=TRUE)
{
  if(is.null(strings)) return(NA)
  if(length(strings)==1)
    if(is.na(strings)) return(NA) else return(strings)
  allsame<-TRUE
  substs<-rep(NA,length(strings))
  if(leading) #this is for detecting the common leading characters
  {count<-1
  while(allsame==TRUE)
  {for(st in strings)
    substs[which(st==strings)]<-substr(st,1,count)
  if(!test_if_equal(substs)) allsame<-FALSE else
    if(count==min(sapply(strings,nchar))) {allsame<-FALSE; count<-count+1}
  if(allsame) count<-count+1
  if(!allsame) count<-count-1
  }
  if(count==0) return(NA)
  out_string<-substr(strings[1],1,count)}
  
  if(!leading) #this is for detecting the common trailing characters
  {backcount<-1
  while(allsame==TRUE)
  {for(st in strings)
  {end_string<-nchar(st)
  substs[which(st==strings)]<-substr(st,end_string-backcount+1,end_string)}
    if(!test_if_equal(substs)) {allsame<-FALSE; backcount<-backcount-1} else
      if(backcount==min(sapply(strings,nchar))) allsame<-FALSE
      if(allsame) backcount<-backcount+1}
  
  if(backcount==0) return(NA)
  end_string<-nchar(strings[1])
  out_string<-substr(strings[1],end_string-backcount+1,end_string)}
  
  return(out_string) 
}