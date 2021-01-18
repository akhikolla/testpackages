
# # example:
# n <- 10 # number of units to rank
# n2 <- 4 # the first n2 units are tied, demonstrating how ties are handled
# x <- rnorm(n) # make random data
# if(n2>1) { # implement replication of first n2 units
#   for(i in 2:n2) {
#     x[i] <- x[1]
#   }
# }
# rk <- rank(x, ties.method="average") # traditional rank with average
# rk2 <- wrank(x)
# all.equal(rk, rk2) # the ranks are equal when weights are all 1
# w <- rchisq(n, df=1) # make random ranks with strictly positive values
# rk3 <- wrank(x, w)
# cor(rk, rk3, method="spearman") # ordering is preserved
# cbind(x, rk, rk2, rk3, w)



# # example 2:

#x <- sample(runif(4), 8, replace=TRUE)
#w <- sample(runif(4), 8, replace=TRUE)
#cbind(x,w, wCorr:::wrank(x,w),wrank2(x,w), wrank2(x,w)- wCorr:::wrank(x,w))
#plot(wrank2(x,w), wCorr:::wrank(x,w))
#abline(0,1)
#all.equal(wrank2(x,w), wCorr:::wrank(x,w))
# 
#n <- 40000
#system.time(wrank2(1:n,rep(1,n)))
#system.time(wCorr:::wrank(1:n,rep(1,n)))


wrank <- function(x, w=rep(1,length(x))) {
  # sort by x so we can just traverse once
  ord <- order(x)
  rord <- (1:length(x))[order(ord)] # reverse order
  xp <- x[ord] # x, permuted
  wp <- w[ord] # weights, permuted
  rnk <- rep(NA, length(x)) # blank ranks vector
  # setup first itteration
  t1 <- 0 # total weight of lower ranked elements
  i <- 1 # index
  t2 <- 0 # total weight of tied elements (including self)
  n <- 0 # number of tied elements
  while(i < length(x)) {
    t2 <- t2 + wp[i] # tied weight increases by this unit
    n <- n + 1 # one additional tied unit
    if(xp[i+1] != xp[i]) { # the next one is not a tie
      # find the rank of all tied elements
      rnki <- t1 + (n+1)/(2*n)*t2
      # push that rank to all tied units
      for(ii in 1:n) {
        rnk[i-ii+1] <- rnki
      }
      # reset for next itteration
      t1 <- t1 + t2 # new total weight for lower values
      t2 <- 0 # new tied weight starts at 0
      n <- 0 # no tied units
    }
    i <- i + 1
  }
  # final row
  n <- n + 1 # add final tie
  t2 <- t2 + wp[i] # add final weight to tied weight
  rnki <- t1 + (n+1)/(2*n)*t2 # final rank
  # push that rank to all final tied units
  for(ii in 1:n) {
    rnk[i-ii+1] <- rnki
  }
  # order by incoming index, so put in the original order
  rnk[rord]
}

wrank_old <- function(x, w=rep(1,length(x))) {
  sapply(1:length(x), function(i) {
    t1 <- sum(w[x<x[i]]) # ranked below every unit below it
    t2 <- w[x==x[i]] #  ties
    # Note: when selecting the range to average over you have to figure out which unit is first.
    #       this method assumes that all of the units are exchangable and integrates over all
    #       units in the tie.
    # mean(t2) brings the unit up to the smallest rank
    # sum(t2) - mean(t2) /2 is then the middle of the ranks
    t1  + mean(t2) + (sum(t2) -mean(t2))/2
  })
}
