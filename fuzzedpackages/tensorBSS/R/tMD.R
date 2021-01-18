tMD <- function(W.hat, A){
  r <- length(W.hat)
  if(r == 1){
    return(MD(W.hat[[1]], A[[1]]))
  }
  
  Wk <- W.hat[[1]]
  Ak <- A[[1]]
  for(m in 2:r){
    Wk <- W.hat[[m]]%x%Wk
    Ak <- A[[m]]%x%Ak
  }
  
  MD(Wk, Ak)
}
