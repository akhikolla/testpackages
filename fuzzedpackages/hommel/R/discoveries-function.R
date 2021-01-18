discoveries <- function(hommel, ix, incremental=FALSE, alpha=0.05)
{
  m <- length(hommel@p)
  if (missing(ix) & incremental==FALSE) {
    p = hommel@p
    k = m
  } 
  if (missing(ix) & incremental==TRUE) {
    stop('Found incremental=TRUE but missing ix.')
  }
  
  if (!missing(ix)) {
    p <- hommel@p[ix]
    k <- length(p)
  }

  if (any(is.na(p)))
    stop('NAs produced by selecting with ix.')

  if (k == 0) {
    warning('empty selection')
    return(0)
  }

  h <- findHalpha(hommel@jumpalpha, alpha, m)

  simesfactor <- hommel@simesfactor[h+1]

  allsortedp <- hommel@p[hommel@sorter]

  discoveries <- findDiscoveries(p, allsortedp, simesfactor, h, alpha, k, m)

  if(!incremental) 
    return(discoveries[k+1])
  else 
    return(discoveries[-1])

}


tdp <- function(hommel, ix, incremental=FALSE, alpha=0.05)
{
  m <- length(hommel@p)
  if (missing(ix)) {
    d <- discoveries(hommel, incremental=incremental, alpha=alpha)
    k <- m
  } else {
    p <- hommel@p[ix]
    k <- length(p)
    d <- discoveries(hommel, ix, incremental=incremental, alpha=alpha)
  }
  d/k
}

fdp <- function(hommel, ix, incremental = FALSE, alpha=0.05)
{
  1-tdp(hommel, ix, incremental=incremental, alpha=alpha)
}
