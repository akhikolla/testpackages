pairwiseDiffs <- function(x){
  x <- c(0,cumsum(x))
  x1 <- x[1:(length(x)-1)]
  x2 <- x[2:length(x)]
  (x1+x2)/2
}

trim.leading <- function (x)  sub("^\\s+", "", x)

chrOrder <- function(x){
  x.num <- suppressWarnings(as.numeric(x))
  nonnum.pos <- which(is.na(x.num))
  num.pos <- which(!is.na(x.num))
  
  x.num <- x[num.pos]
  x.nonnum <- x[nonnum.pos] 
  output <- c(sort(as.numeric(x.num)), sort(x.nonnum))
  
  output
}