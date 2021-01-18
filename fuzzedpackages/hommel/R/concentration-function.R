concentration <- function(hommel, alpha = 0.05) {
  
  m <- length(hommel@p)
  
  h <- findHalpha(hommel@jumpalpha, alpha, m)
  simesfactor <- hommel@simesfactor[h+1]

  sortedp <- hommel@p[hommel@sorter]
  
  z <- findConcentration(sortedp, simesfactor, h, alpha, m)
  
  return(sortedp[z])

}
  