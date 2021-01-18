#' Select string that end in a particular way (e.g. a certain file extension)
#' 
#' Sometimes it makes sense to apply a function to several files in a folder, but only to
#' those of a particular file type. This function can selects all elements in a vector of strings
#' that end in a particular way, e.g. on a common file extension.
#' 
#' @param strings vector of character strings for elements to be extracted from.
#' @param file_extension character string specifying the extension of the file type to be selected.
#' This can also be any other trailing string that marks all vector elements to be selected.
#'  
#' @return subset of the strings vector that only contains the elements that end on file_extension.
#' 
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#'   select_by_file_extension(c("Temp1.csv","Temp1.xls","Temp2.csv","Temp2.xls"),"csv")
#'   select_by_file_extension(c("red car","blue car","yellow duck"), "car")
#'  
#' @export select_by_file_extension
select_by_file_extension<-function(strings,file_extension)
 {ends<-sapply(strings,function(x) {ll<-nchar(x); substr(x,ll-nchar(file_extension)+1,ll)})
  outstrings<-strings[which(ends==file_extension)]
  if(length(outstrings)==0) return(NA)
  return(outstrings)
}
