change.param <- function(coefficients,design_list){
  
  ## get basic parameters
  m <- design_list$m
  n <- design_list$n
  I <- design_list$I
  n_sigma <- design_list$n_sigma
  q <- design_list$q
  
  ## extract different parameter matrices
  ## basic model parameters, called delta 
  p.delta <- I*q
  if(design_list$RSM){
    p.delta <- I+q
  }
  delta <- coefficients[,1:p.delta,drop=FALSE]
  
  ## DIF parameters, called gamma
  if(design_list$DSF){
    gamma <- coefficients[,(p.delta+1):(p.delta+I*m*q),drop=FALSE]
    start.sigma <- p.delta+I*m*q+1
  }else{
    gamma <- coefficients[,(p.delta+1):(p.delta+I*m),drop=FALSE]
    start.sigma <- p.delta+I*m+1
  }
  
  ## Discrimination parameters, called sigma
  sigma <- coefficients[,start.sigma:(start.sigma+n_sigma-1),drop=FALSE]
  
}