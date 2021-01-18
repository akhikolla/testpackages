
## to get the date NOT in the maintainer's French locale... 
.getDocDate <- function() {
  lct <- Sys.getlocale('LC_TIME') 
  Sys.setlocale('LC_TIME', 'C')
  dat <-format(Sys.time(), 'This documentation: %d %B %Y')
  Sys.setlocale('LC_TIME', lct)
  return(dat)
}