#' Computes the posterior mean hazard values for a vector x for the Piecewise Exponential Hazard model (PEH)
#' @param x Vector of times to compute the posterior mean hazard.
#' @param G1 List of posterior samples from the BayesPiecewiseHazard function.
#' @return Vector containing the posterior mean hazard values h(x)
#'@export
PostMeanHazPiece = function(x,G1){
  GetHazPEH = function(x,s,lam,J){
    
    y=x
    
    
    
    for(m in 1:length(x)){
      for(k in 1:(J+1)){
        if((x[m]>s[k]) && (x[m]<=s[k+1])){
          y[m]=lam[k]
          
          
        }
      }
      
    }
    
    return(y)
  }
  
  Store=rep(NA,nrow(G1[[1]]))
  
  y1=rep(0,length(x))
  y=x
  
  for(b in 1:nrow(G1[[1]])){
    s=G1[[1]][b,]
    lam=(G1[[2]])[b,]
    J = G1[[3]][b]
    
    
    y=GetHazPEH(x,s,lam,J)
    
    
    
    
    
    y1=y1+    y
  }
  
  return(y1/nrow(G1[[1]]))
  
}
