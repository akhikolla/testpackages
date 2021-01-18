#' Simulate data from a Kalman model
#'
#' @name Generate_Data 
#'
#' @description This function simulates data obeying a Kalman model whilst allowing the user to add innovative and additive anomalies. 
#' 
#' @param n A positive integer giving the number of observations desired
#' @param A A matrix giving the updates for the hidden states. 
#' @param C A matrix mapping the hidden states to the observed states.
#' @param Sigma_Add A positive definite diagonal matrix giving the additive noise covariance.
#' @param Sigma_Inn  A positive definite diagonal matrix giving the innovative noise covariance.
#' @param mu_0 A matrix indicating the mean of the prior for the hidden states. It defaults to a zero-vector.
#' @param anomaly_loc A vector of integers giving the locations of anomalies.
#' @param anomaly_type A vector of strings, either "Add" or "Inn" indicating whether the anomaly is additive or innovative.
#' @param anomaly_comp A vector of integers giving the component affected by the anomalies.
#' @param anomaly_strength  A vector of numerics giving the strength of the anomalies (in sigmas).
#' @return A list of matrices, each corresponding to an observation. 
#'
#' @examples
#' 
#' library(RobKF)
#' library(ggplot2)
#' 
#' set.seed(2018)
#' 
#' A = diag(2)*0.99
#' A[1,2] = -0.05
#' C = matrix(c(10,0.1),nrow=1)
#' mu = matrix(c(0,0),nrow=2)
#' Sigma_Inn = diag(c(1,0.01)*0.00001,nrow=2)
#' Sigma_Add = diag(c(1)*0.1,nrow=1)
#' 
#' Y_list = Generate_Data(100,A,C,Sigma_Add,Sigma_Inn, mu_0 = mu,  anomaly_loc = c(10,30,50), 
#'                       anomaly_type = c("Inn","Add","Inn"), 
#'                       anomaly_comp = c(1,1,2),  anomaly_strength = c(400,-10,3000))
#'                       
#' qplot(1:100,unlist(Y_list),xlab="time",ylab="observation")+theme_minimal()
#' 
#' 
#' @export
#' 
#' 
#' @export
Generate_Data = function(n,A,C,Sigma_Add,Sigma_Inn,mu_0 = NULL,  anomaly_loc = integer(0), anomaly_type = character(0), 
                         anomaly_comp = integer(0),  anomaly_strength = NULL){
  
  p = nrow(Sigma_Add)
  q = nrow(Sigma_Inn)
  
  n = as.integer(n)
  
  if (n < 1){
    stop("n must be positive.")
  }
  
  if (is.null(mu_0)){
    mu_0 = rep(0, nrow(Sigma_Inn))
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
  
  if (sum(Sigma_Add) != sum(abs(Sigma_Add))){
    stop("All entried of Sigma_Add need to be positive")
  }
  
  if (sum(Sigma_Add) != sum(diag(Sigma_Add))){
    stop("Sigma_Add must be diagonal")
  }
  
  if (sum(Sigma_Inn) != sum(abs(Sigma_Inn))){
    stop("All entried of Sigma_Inn need to be positive")
  }
  
  if (sum(Sigma_Inn) != sum(diag(Sigma_Inn))){
    stop("Sigma_Inn must be diagonal")
  }
  
  if (nrow(C) != p){
    stop("C must have the same number of rows as the observations")
  }
  
  if (ncol(C) != q){
    stop("The number of columns of C must equal the dimensions of the hidden state")
  }
  
  if (nrow(A) != q){
    stop("A must have the same number of rows as the hidden states")
  }
  
  if (ncol(A) != q){
    stop("The number of columns of A must equal the dimensions of the hidden state")
  }
  
  if (is.null(anomaly_strength)){
    anomaly_strength = rep(1000,length(anomaly_comp))
  }
  
  if (length(anomaly_strength) != length(anomaly_comp)){
    stop("anomaly_strength, anomaly_comp, anomaly_type, and anomaly_loc must be of the same length")
  }
  
  if (length(anomaly_strength) != length(anomaly_type)){
    stop("anomaly_strength, anomaly_comp, anomaly_type, and anomaly_loc must be of the same length")
  }
  
  if (length(anomaly_strength) != length(anomaly_loc)){
    stop("anomaly_strength, anomaly_comp, anomaly_type, and anomaly_loc must be of the same length")
  }
  
  anomaly_loc = as.integer(anomaly_loc)
  
  if (length(anomaly_loc) > 0){
    
    if (max(anomaly_loc)>n){
      stop("The entries of anomaly_loc must be between 1 and n")
    }
  
    if (min(anomaly_loc)<1){
      stop("The entries of anomaly_loc must be between 1 and n")
    }
    
  }
  
  if (sum(anomaly_type %in% c("Add","Inn")) < length(anomaly_type)){
    stop("The entries of anomaly_type must be either Add or Inn")
  }
  
  anomaly_comp = as.integer(anomaly_comp)
  
  if (sum(anomaly_comp[which(anomaly_type == "Add")] %in% 1:p) != sum(anomaly_type == "Add")){
    stop("Out of bounds for anomaly_comp ")
  }
  
  if (sum(anomaly_comp[which(anomaly_type == "Inn")] %in% 1:q) != sum(anomaly_type == "Inn")){
    stop("Out of bounds for anomaly_comp ")
  }
  
  Innovations =  Sigma_Inn^(1/2) %*% matrix(rnorm(n*nrow(Sigma_Inn)),nrow = nrow(Sigma_Inn))
  Additions =  Sigma_Add^(1/2) %*% matrix(rnorm(n*nrow(Sigma_Add)),nrow = nrow(Sigma_Add))
  
  if (length(anomaly_loc)>0){
  
    for (ii in 1:length(anomaly_loc)){
    
    
      if (anomaly_type[ii] == "Inn"){
        Innovations[anomaly_comp[ii],anomaly_loc[ii]] =  (Sigma_Inn[anomaly_comp[ii],anomaly_comp[ii]]^(1/2))*anomaly_strength[ii]
      }
    
      if (anomaly_type[ii] == "Add"){
        Additions[anomaly_comp[ii],anomaly_loc[ii]] = (Sigma_Add[anomaly_comp[ii],anomaly_comp[ii]]^(1/2))*anomaly_strength[ii]
      }
    
    }
    
  }
  
  Y_list = list()
  
  hidden_stats = matrix(0,ncol = ncol(Innovations)+1,nrow = nrow(Innovations))
  hidden_stats[,1] = mu_0
  
  for (ii in 1:ncol(Innovations)) {
    hidden_stats[,ii+1] = A %*% matrix(hidden_stats[,ii],ncol=1) + matrix(Innovations[,ii],ncol = 1)
  }
  
  for (ii in 1:ncol(Additions)) {
    Y_list = c(Y_list,list(C %*% matrix(hidden_stats[,ii+1],ncol=1) + matrix(Additions[,ii],ncol = 1) ))
  }
  
  return(Y_list)
  
}