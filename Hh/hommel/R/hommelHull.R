hommelHull <- function(pvalues, simes = TRUE) {
  # TODO: check if arguments are sound

  m <- length(pvalues)
  if(m <= 0)
    return(pvalues)

  # sort pvalues
  perm <- order(pvalues)
  p <- pvalues[perm]
  
  # compute the lower convex hull of points (i,p[i]) by the monotone chain algorithm
  isconvex <- function(i,j,k,pi=p[i])  (j-i) * (p[k]-pi) > (k-i) * (p[j]-pi)
  hull <- numeric(m)
  n <- 0
  
  for(i in 1:m) {
    while(n >= 2 && !isconvex(hull[n-1], hull[n], i))
      n <- n-1
    hull[n <- n+1] <- i
  }
    
  # compute the description a[i] of h determined by h(a) >= i iff a < a[i]
  # we have a[i] = s[i] * slope[i], where s[i] is the Simes factor, and slope[i] = p{m-i+j} / j is minimal
  s <- ifelse(rep(simes,m), 1:m, 1:m*cumsum(1/1:m))
  a <- numeric(m)

  for(i in 1:m) {
    while(n >= 2 && isconvex(m-i,pi=0,hull[n-1],hull[n]))
      n <- n-1
    a[i] <- s[i] * p[hull[n]] / (hull[n] - m + i)
  }

  # calculate adjusted pvalues
  adj <- numeric(m)
  i <- 1
  j <- m+1
  a[j] <- s[j] <- 0
  while(i <= m && j >= 1) 
    if(j == 1 || s[j-1]*p[i] <= a[j]) {
      adj[i] <- min(s[j]*p[i], a[j])
      i <- i+1
    }
    else
      j <- j-1

  # back to original order
  adj[perm] <- adj
  return(adj)
}




# TODO: test corner cases like pvalues equal to zero, other distributions, large m
test <- function(pvals) {
  adj1 <- hommelHull(pvals) 
  adj2 <- p.adjust(pvals, method="hommel")
  
  if(!isTRUE(all.equal(adj1,adj2)))
    stop("Different adjusted pvalues for pvalues\n", paste(pvals,collapse=", "))
}

test(rep(0,100))
test(rep(0.5,100))
test(runif(1000)*1e-9)
test(0.5+runif(1000)*1e-9)

sizes <- c(0, sample.int(20, size=10000, replace=T), sample.int(1000, size=100, replace=T))
for(size in sizes) {
  pvals <- runif(size)
  test(pvals)
}

# only for a quick idea
for(size in 10^(2:6))
  print(system.time(hommelHull(runif(size))))

