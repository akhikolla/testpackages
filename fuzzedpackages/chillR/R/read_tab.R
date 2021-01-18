#' Read csv table regardless of whether it is a true csv or the French type
#' 
#' csv tables are widely used for storing data as 'comma-separated values'.
#' This doesn't work, however, when the comma is also used as a decimal symbol,
#' as is practiced in French or German, for example. The separator symbol for csv
#' files then becomes a semi-colon. This is not problematic when you only work on one
#' machine, but it causes problems when you collaborate with people who use different
#' types of csv encoding.
#' 
#' This function overcomes this problem by checking first, which of the two
#' characters occurs most frequently in the table, assuming then that this is
#' the separator symbol. It then opens the table accordingly.
#' 
#' Currently limited to files that are either comma-separated with point as decimal
#' symbol or semicolon-separated with comma as decimal symbol. Files should also
#' have a header.
#' 
#' @param tab file name of a table to be read.
#' @return If the table is in one of the two formats described above, the stored table
#' is returned.
#' @author Eike Luedeling
#' @keywords utility
#' @examples
#' 
#' df<-data.frame(Var1=c(1,2,3.2,1.2),Var2=c(1.2,6,2.6,7))
#' write.csv(df,"filecsv.csv",row.names=FALSE)
#' read_tab("filecsv.csv")
#' write.table(df,"filesemicolon.csv",sep=";",dec=",")
#' read_tab("filesemicolon.csv")
#' file.remove("filecsv.csv")
#' file.remove("filesemicolon.csv")
#' 
#'  
#' @export read_tab
read_tab <-function(tab)
{
  content<-readChar(tab, file.info(tab)$size)
  freqcomma<-length(gregexpr(",", content)[[1]])
  freqsemicolon<-length(gregexpr(";", content)[[1]])
  
  if(freqcomma>freqsemicolon) return(read.table(tab,header=TRUE,sep=",",dec="."))
  if(freqcomma<freqsemicolon) return(read.table(tab,header=TRUE,sep=";",dec=","))}

