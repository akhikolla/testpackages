Viterbi <- function(y, z, p){
  logA <- log(p@A[[z]])
  m <- matrix(0, length(y), p@M)
  whereNotZero <- which(y!=0)
  all.dgamma <- matrix(log(p@lambda$eps), length(y), p@M, byrow = T)
  if (any(y!=0)) all.dgamma[whereNotZero,] <- sapply(1:p@M, function(h) dgamma(y[whereNotZero], p@lambda$a[h], p@lambda$b[h], log = TRUE) + log(1-p@lambda$eps[h]))
  m[1, ] <- log(p@pi[z,]) + all.dgamma[1,]
  if (length(y)>1){
    for (t in 2:length(y)){
      for (l in 1:p@M){
        m[t,l] <- max(m[t-1,] + logA[,l]) + all.dgamma[t,l]
      }
    }
  }
  xMAP <- y * 0
  xMAP[length(y)] <- which.max(m[length(y), ])
  if (length(y)>1){
    for (t in (length(y)-1):1){
      xMAP[t] <- which.max(m[t,] + logA[,xMAP[t+1]])
    }
  }
  xMAP
}


