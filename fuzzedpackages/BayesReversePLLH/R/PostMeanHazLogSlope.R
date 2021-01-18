#' Computes the posterior mean hazard value for a vector x for the Piecewise Linear Log Hazard model (PLLH)
#' @param x Vector of times to compute the posterior mean hazard function
#' @param G1 List of posterior samples from the BayesPiecewiseLinearLogHazard function.
#' @return Vector containing the posterior mean hazard values h(x)
#'@export
PostMeanHazLogSlope = function(x,G1){
  
  
  GetHazPLLH = function(x,s,lam,J){
    y=x
    
    slopes = diff(lam)/diff(s)
    slopes=slopes[1:(J+1)]
    
    for(m in 1:length(x)){
      for(k in 1:(J+1)){
        if((x[m]>s[k]) && (x[m]<=s[k+1])){
          y[m]=(x[m]-s[k])*slopes[k]+lam[k]
          
          
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
    lam=G1[[2]][b,]
    J = G1[[3]][b]
    
    y=GetHazPLLH(x,s,lam,J)
    
    ##Add up hazard
    y1=y1+y
    
  }
  
  return(y1/nrow(G1[[1]]))
  
}
