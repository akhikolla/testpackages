get_sigma_hat = function(Final_Sigma,C,A,Sigma_Add,Sigma_Inn,Sigma_Inn_list,horizon_matrix,p,q){
  
  candidates = matrix(0,nrow = nrow(horizon_matrix),ncol=ncol(horizon_matrix))
  
  for (jj in 1:nrow(horizon_matrix)){
    
    candidates[jj,] = diag((  t(C[[jj]]) %*%  solve(C[[jj]] %*% ( A %*% Final_Sigma %*% t(A) ) %*% t(C[[jj]]) + Sigma_Add[[jj]] + Sigma_Inn_list[[jj]]) %*% C[[jj]] ))
    
    candidates[jj,which(horizon_matrix[jj,] == 0)] = 0 
  
  }
    
  sigmas = diag(Sigma_Inn)*apply(candidates,max,MARGIN=2)
  
  return(sigmas)
  
}