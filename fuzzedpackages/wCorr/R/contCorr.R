contCorr <- function(x,y,w, method=c("Pearson", "Spearman")) {
  if(!is.numeric(x)) { 
    x <- as.numeric(x)
  }
  if(!is.numeric(y)) {
    y <- as.numeric(y)
  }
    if(!is.numeric(w)) { 
    w <- as.numeric(w)
  }
  pm <- pmatch(tolower(method[[1]]), tolower(c("Pearson", "Spearman")))
  if(pm == 2) {
    #x <- rank(x) # rank gives averages for ties
    #y <- rank(y)
    x <- wrank(x,w)
    y <- wrank(y,w)
  }
  xb <- sum(w*x)/sum(w)
  yb <- sum(w*y)/sum(w)
  numerator <- sum(w*(x-xb)*(y-yb))
  denom <- sqrt( sum(w*(x-xb)^2) * sum(w*(y-yb)^2))
  numerator/denom
}
