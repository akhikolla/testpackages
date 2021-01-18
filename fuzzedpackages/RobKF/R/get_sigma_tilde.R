get_sigma_tilde = function(Final_Sigma,C,A,Sigma_Add,Sigma_Inn){
  
  sigmas = rep(NA,nrow(Sigma_Add))
  
  Sigma_Pred = C %*% ( A %*% Final_Sigma %*% t(A) + Sigma_Inn) %*% t(C) + Sigma_Add
  
  Sigma_Pred_Inverse = solve(Sigma_Pred)
  
  for (ii in 1:length(sigmas)){
    
    sigmas[ii] = Sigma_Add[ii,ii]*Sigma_Pred_Inverse[ii,ii]
    
  }
  
  return(sigmas)
}