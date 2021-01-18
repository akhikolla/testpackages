"sim2dis"<-function(f, a=4.0, baseF=0.03, method=1){
  f  <- baseF + (1-baseF)*f

  if(method==1){
    D  <- (-log(f))^a
  # D1 <- diag(D)%*%t(rep(1,ncol(D)))
  # D  <- D-(t(D1)+D1)/2
  # D  <- pmax(D,0)
  }
  
  if(method==2){
    F1 <- diag(f)%*%t(rep(1,ncol(f)))
    D  <- sqrt(pmax((t(F1)+F1)/2 - f,0))^a
  }

  D
}