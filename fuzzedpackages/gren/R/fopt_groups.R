# objective function for penalty multiplier estimation
fopt_groups <- function(par, lambda2, nparts, partsind, partsmat, sizes, G, 
                        sum1) {
  
  s <- par[1]
  loglambdag <- par[-1]
  
  loglambdasum <- rowSums(sapply(1:nparts, function(part) {
    loglambdag[partsind==part][partsmat[, part]]}))
  
  partsum <- sum((0.5*lambda2*unlist(sapply(1:nparts, function(part) {
    tapply(exp(loglambdasum)*sum1, partsmat[, part], sum)})) + (s - 0.5)*sizes)^2)
  constr <- sum(loglambdag*sizes)^2
  magn <- partsum + constr
  return(magn)
  
}