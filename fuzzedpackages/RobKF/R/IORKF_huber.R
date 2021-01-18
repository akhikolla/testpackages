#' A huberisation based innovative outlier robust Kalman filter
#'
#' @name IORKF_huber
#'
#' @description An innovative outlier robust Kalman filter, based on the work by Ruckdeschel et al. (2014).
#' This function assumes that the innovations are potentially polluted by a heavy tailed process.
#' The update equations are made robust to these via huberisation.
#' 
#' @param Y A list of matrices containing the observations to be filtered.
#' @param mu_0 A matrix indicating the mean of the prior for the hidden states. 
#' @param Sigma_0 A matrix indicating the variance of the prior for the hidden states. It defaults to the limit of the variance of the Kalman filter.
#' @param A A matrix giving the updates for the hidden states. 
#' @param C A matrix mapping the hidden states to the observed states.
#' @param Sigma_Add A positive definite matrix giving the additive noise covariance.
#' @param Sigma_Inn  A positive definite matrix giving the innovative noise covariance.
#' @param h A numeric giving the huber threshold. It defaults to 2. 
#' @param epsilon A positive numeric giving the precision to which the limit of the covariance is to be computed. It defaults to 0.000001.
#' @return An rkf S3 class. 
#' @references \insertRef{ruckdeschel2014robust}{RobKF}
#'
#' @examples
#' library(RobKF)
#' 
#' set.seed(2019)
#'
#' A = matrix(c(1), nrow = 1, ncol = 1)
#' C = matrix(c(1), nrow = 1, ncol = 1)
#'
#' Sigma_Inn = diag(1,1)*0.01
#' Sigma_Add = diag(1,1)
#'
#' mu_0 = matrix(0,nrow=1,ncol=1)
#'
#' Y_list = Generate_Data(1000,A,C,Sigma_Add,Sigma_Inn,mu_0,anomaly_loc = c(100,400,700),
#'                        anomaly_type = c("Inn","Inn","Inn"),anomaly_comp = c(1,1,1),
#'                        anomaly_strength = c(50,80,-100))
#'
#' Output = IORKF_huber(Y_list,mu_0,Sigma_0=NULL,A,C,Sigma_Add,Sigma_Inn,h=2)
#'
#' plot(Output,conf_level = 0.9999)
#' 
#' @export
IORKF_huber = function(Y,mu_0,Sigma_0=NULL,A,C,Sigma_Add,Sigma_Inn,h=2,epsilon=0.000001)
{
  
  p = nrow(Sigma_Add)
  q = nrow(Sigma_Inn)
  
  if (is.null(p)){
    stop("Sigma_Add must be a matrix.")
  }
  
  if (is.null(q)){
    stop("Sigma_Inn must be a matrix.")
  }
  
  for (ii in 1:length(Y)){
    
    Y[[ii]] = as.matrix(Y[[ii]])
    
    if (nrow( Y[[ii]]) != p){
      stop("all observations must be of the same number of rows as Sigma_Add.")
    }
    if (ncol( Y[[ii]]) != 1){
      stop("all observations must have exactly one column,")
    }
    
  }
  
  mu_0 = as.matrix(mu_0)
  
  if (nrow(mu_0) != q){
    stop("mu_0 must be of the same number of rows as Sigma_Add.")
  }
  
  if (ncol(mu_0) != 1){
    stop("mu_0 must have 1 column.")
  }
  
  if (ncol(Sigma_Add) != p){
    stop("Sigma_Add needs to be a square matrix.")
  }
  
  if (ncol(Sigma_Inn) != q){
    stop("Sigma_Inn needs to be a square matrix.")
  }
  
  if (nrow(C) != p){
    stop("C must have the same number of rows as the observations")
  }
  
  if (ncol(C) != q){
    stop("The number of columns of C must equal the dimensions of the hidden state")
  }
  
  tryCatch(
    expr = {
      solve(C)
    },
    error = function(e){ 
      stop("C must be invertible")
    }
  )
  
  if (nrow(A) != q){
    stop("A must have the same number of rows as the hidden states")
  }
  
  if (ncol(A) != q){
    stop("The number of columns of A must equal the dimensions of the hidden state")
  }
  
  if (epsilon <= 0){
    stop ("epsilon must be greater than 0.")
  }
  
  if (!isSymmetric(Sigma_Inn)){
    stop("Sigma_Inn must be positive definite.")
  }
  
  if (sum(eigen(Sigma_Inn)$values <= 0) > 0 ){
    stop("Sigma_Inn must be positive definite.")
  }
  
  if (!isSymmetric(Sigma_Add)){
    stop("Sigma_Add must be positive definite.")
  }
  
  if (sum(eigen(Sigma_Add)$values <= 0) > 0 ){
    stop("Sigma_Add must be positive definite.")
  }
  
  if(is.null(Sigma_0)){
    
    New_Matrix = C
    
    Full_Matrix = C
    
    Rank = 0
    
    repeat {
      
      Rank = Rank+1
      
      if(Rank > q+5) {
        break
      }
      if( rankMatrix(Full_Matrix)  == q ) {
        break
      }
      
      New_Matrix =  New_Matrix %*% A
      
      Full_Matrix = rbind(Full_Matrix,New_Matrix)
      
    }
    
    if (Rank >q){
      stop("The system has to be observable to infer Sigma_0.")
    }
    
    
    Final_Sigma = Sigma_Limit(diag(1,nrow = q),C,A,Sigma_Inn,Sigma_Add,epsilon)
    
    Sigma_0 = Final_Sigma
    
  } 
  
  if (nrow(Sigma_0) != q){
    stop("Sigma_0 must have the same number of rows as the hidden states.")
  }
  
  if (ncol(Sigma_0) != q){
    stop("The number of columns of Sigma_0 must equal the dimensions of the hidden state.")
  }
  
  if (!isSymmetric(Sigma_0)){
    stop("Sigma_0 must be positive definite.")
  }
  
  if (sum(eigen(Sigma_0)$values <= 0) > 0 ){
    stop("Sigma_0 must be positive definite.")
  }
  
  if (h <= 0){
    stop ("h must be positive.")
  }
  
  b = matrix(rep(0,q),ncol=1)
  d = matrix(rep(0,p),ncol=1)
  
  algo_output  = iorkf_huber_list(mu_0,Sigma_0,Y,A,b,C,d,Sigma_Add,Sigma_Inn,h)
  
  if (is.null(algo_output)){
    stop("User interrupt.")
  }
  
  output = list()
  
  output[["Y"]] = Y
  
  output[["A"]] = A
  
  output[["C"]] = C
  
  output[["Sigma_Inn"]] = Sigma_Inn
  
  output[["Sigma_Add"]] = Sigma_Add
  
  output[["States"]] = algo_output
  
  output["Type"] = "IO"
  
  return(structure(output,class="rkf"))
  
}