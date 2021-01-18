overmessage <- function(msg, prevmsglength) {
  msglength <- nchar(msg)
  if (prevmsglength>0) {message("\r", appendLF=F)}  	##FR: for backslash-b see ?Quotes ...
  message(msg, appendLF=F)
  return(msglength)
}
