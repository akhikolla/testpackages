## internal function to get absolute path name under Unix and windows
## but only if absolute path has not already been specified. 
fullpath <- function(fname){
 ## check if full path has been given
  filen <- fname  ## initialize
  if (! (length(grep("/", fname)) > 0 || length(grep("[\\]", fname)) > 0 ) ){
     if(.Platform$OS.type == "unix") {
       filen <- paste(getwd(), "/", fname, sep="")
     } else {
      filen <- paste(getwd(), "\\", fname, sep="")
     }
  }
  return(filen)
}


