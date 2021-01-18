localtest <- function(hommel, ix, tdp=0)
{
  m <- length(hommel@p)
  if (missing(ix)) {
    p = hommel@p
    n = m
  }
  else {
    p <- hommel@p[ix]
    n <- length(p)
  }

  if (tdp<0 | tdp>=1)
    stop("'tdp' must be chosen from [0,1)")

  if (any(is.na(p)))
    stop('NAs produced by selecting with ix')

  if (n == 0)
    stop('empty selection')

  k <- floor(tdp*n)
  sortedp <- sort(p)
  pI <- sortedp[seq_len(n-k)+k]
  pI <- min(pI/rank(pI, ties.method='first'))

  adjustedp <- adjustedIntersection(pI, hommel@jumpalpha, m, hommel@simesfactor)

  return(adjustedp)
}
