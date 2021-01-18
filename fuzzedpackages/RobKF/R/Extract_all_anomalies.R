Extract_anomaly = function(x,A,C,Sigma_Add,Sigma_Inn){
  
  Innovation = x$Y - C %*% A %*% x$mu 
  
  Predictive_Var = Sigma_Add + C %*% ( A %*% x$Sigma %*% t(A) + Sigma_Inn ) %*% t(C)
  
  score = t(Innovation) %*% solve(Predictive_Var) %*% Innovation
  
  return(score)
  
}


Extract_all_anomalies = function(x){
  
  new_list = list()
  
  for (ii in 1:length(x[["Y"]])){
    
    add = list()
    
    add[["Y"]]     = x[["Y"]][[ii]]
    
    add[["mu"]]    = x[["States"]][[ii]][[1]]
 
    add[["Sigma"]] = x[["States"]][[ii]][[2]]
       
    new_list[[ii]] = add
    
  }
  
  scores = unlist(lapply(new_list,Extract_anomaly,A=x[["A"]],C=x[["C"]],Sigma_Add=x[["Sigma_Add"]],Sigma_Inn=x[["Sigma_Inn"]]))
  
  return(scores)
  
}