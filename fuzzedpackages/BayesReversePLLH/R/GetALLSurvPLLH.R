#' Computes posterior distribution of survival probabilities for a vector x for the Piecewise Linear Log Hazard model (PLLH)
#' @param x Vector of times to compute the posterior mean survival probability.
#' @param G1 List of posterior samples from the BayesPiecewiseLinearLogHazard function.
#' @return Matrix containing the posterior distribution survival probabilities S(x)
#'@export
GetALLSurvPLLH = function(x,G1){
  
  MeanHold=rep(0,length(x))
  SurvHold = matrix(nrow=nrow(G1[[1]]),ncol=length(x))
  for(b in 1:nrow(G1[[1]])){
    s=G1[[1]][b,]
    lam=(G1[[2]])[b,]
    J = G1[[3]][b]
    
    SurvHold[b,]=SurvPLLH(x ,  s,  lam,  J)
    
  }
  
  
  
  return(SurvHold)  
}
