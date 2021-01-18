#' Computes the posterior mean survival probabilities for a vector x for the Piecewise Linear Log Hazard model (PLLH)
#' @param x Vector of times to compute the posterior mean survival probability.
#' @param G1 List of posterior samples from the BayesPiecewiseLinearLogHazard function.
#' @return Vector containing the posterior mean survival probabilities S(x)
#'@export
PostMeanSurvPLLH = function(x,G1){
SurvHold=rep(0,length(x))

for(b in 1:nrow(G1[[1]])){
  s=G1[[1]][b,]
  lam=(G1[[2]])[b,]
  J = G1[[3]][b]
  
  
  SurvHold=SurvHold+  SurvPLLH(x ,  s,  lam,  J)
  
}

return(SurvHold/nrow(G1[[1]]))  
}
