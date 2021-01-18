##  This file runs n=1000 trials on the assumption of equal player
##  strength and for each trial calculates the log-likelihood ratio
##  [for the hypothesis p=equal_strengths, vs the hypothesis
##  p=maxp(Hrand)].  It then draws a histogram of the loglikelihood
##  ratio thus obtained, and superimposes a line showig the observed
##  loglikelihood ratio.

source("counterstrike.R")
equalstrengths <- rep(0.1, 9)
datapoint <- loglik(H,indep(maxp(H))) - loglik(H,equalstrengths)


n <- 1000  # takes about 5 minutes for n=1000 trials
distrib <- rep(0,n)
for(i in seq_len(n)){
  Hrand <- hyper2(pnames=c(team1,team2))
  for(jj in 1:6){
    pp <- rdeath(team1,team2)
    Hrand <- Hrand + counterstrike(team1,team2, pp)
  }
  distrib[i] <- loglik(Hrand,indep(maxp(Hrand))) - loglik(Hrand,equalstrengths)
  print(i)
}

hist(distrib)
abline(v=datapoint,lwd=6)
