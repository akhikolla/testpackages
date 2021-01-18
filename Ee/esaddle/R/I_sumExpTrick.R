
.sumExpTrick <- function(x)
{
  z <- max(x)
  
  out <- sum( exp(x - z) ) * exp(z)
  
  return( out )
  
}

.meanExpTrick <- function(x)
{
  return( .sumExpTrick(x) / length(x) )
}