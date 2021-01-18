#' Test if all character vectors in a string are equal
#' 
#' Compares all elements of a vector of numbers or character strings and returns TRUE
#' if they are all the same, FALSE otherwise.
#' 
#' @param test_vector vector of strings or numbers to be tested. 
#'  
#' @return TRUE if all elements of the vector are the same; FALSE otherwise.
#' 
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#'   test_if_equal(c(1,3,1))
#'   test_if_equal(c("a","a","a"))
#'   test_if_equal(c("a","b","a"))                                            
#'  
#' @export test_if_equal
test_if_equal<-function(test_vector)
{allsame<-TRUE
if(length(test_vector)==1) return(TRUE)
for(i in 2:length(test_vector))
  if(!test_vector[i]==test_vector[1]) allsame<-FALSE
  return(allsame)}
