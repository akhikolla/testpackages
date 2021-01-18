is.binaryGP <- function(object){
  if (inherits(object, "binaryGP")) return(TRUE)
  else return(FALSE)
}
