#' Sort files in a folder, so that numbers are in ascending sequence
#' 
#' Sometimes lists of strings that contain numbers aren't listed automatically in the sequence
#' we would expect, e.g. because numbers below ten are lacking leading zeros (as in
#' c("a1","a10","a100","a11"...)). This function recognizes all shared leading and trailing
#' symbols around the numeric part of such strings and sorts the list according to the embedded
#' numbers.
#' 
#' @param strings vector of character strings to be sorted according to embedded numbers.
#' @param file_extension character string specifying the extension of the file type to be selected.
#' This can also be any other trailing string that marks all vector elements to be selected. This isn't
#' required for the function to run, but may be necessary if the string list of interest contains, for
#' instance, different file types, of which you only want to work with one.
#'  
#' @return subset of the strings vector that only contains the elements that end on file_extension and
#' are sorted in ascending order according to the numeric parts of the strings.
#' 
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#'   ordered_climate_list(c("Temp1_ws30.csv","Temp1_ws30.xls",
#'                          "Temp10_ws30.csv","Temp10_ws30.xls",
#'                          "Temp2_ws30.csv","Temp2_ws30.xls"),"csv")
#'   ordered_climate_list(c("Tx12", "Tx2","Tx4","Tx1"))
#'  
#' @export ordered_climate_list
ordered_climate_list<-function(strings, file_extension=NA)
 {if(!is.na(file_extension)) stringlist<-select_by_file_extension(strings, file_extension) else
  stringlist<-strings
  numbers<-extract_differences_between_characters(stringlist)
  stringlist<-stringlist[order(as.numeric(numbers))]
  return(stringlist)
}
