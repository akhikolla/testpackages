
"getskip" <- function(file){
  Line    <- readLines(file, n=20)
  nSpace  <- sapply(gregexpr("\\s+", Line), function(p) { sum(p>=0) } )
  isTable <- nSpace==nSpace[length(nSpace)]
  skip    <- length(Line) - which.max(rev(!c(FALSE,isTable)))+1
  cat(paste("Using  skip =",skip,"\n"))
  skip
}
