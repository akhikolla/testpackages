#' An innovative and additive outlier robust Kalman filter
#'
#' @name IOAORKF 
#'
#' @description An implementation of Computationally Efficient Bayesian Anomaly detection by Sequential Sampling (CE-BASS) by Fisch et al. (2020).
#' This function assumes that both the innovations and additions are potentially polluted by a heavy tailed process, which is approximated by a t-distribution.
#' To approximate the posterior, particles for the precision (inverse variance) are sampled using a robust approximation to the posterior. Conditionally on those samples, the classical Kalman updates are used.
#' 
#' 
#' @param Y A list of matrices containing the observations to be filtered.
#' @param mu_0 A matrix indicating the mean of the prior for the hidden states. 
#' @param Sigma_0 A matrix indicating the variance of the prior for the hidden states. It defaults to the limit of the variance of the Kalman filter.
#' @param A A matrix giving the updates for the hidden states. 
#' @param C A matrix mapping the hidden states to the observed states.
#' @param Sigma_Add A positive definite diagonal matrix giving the additive noise covariance.
#' @param Sigma_Inn  A positive definite diagonal matrix giving the innovative noise covariance.
#' @param Particles An integer giving the number of particles to be maintained at each step. More particles lead to more accuracy, but also require more memory and CPU time. The parameter should be at least p + q + 1, where p s the dimension of the observations and q the dimension of the hidden states.
#' @param Descendants An integer giving the number of descendants to be sampled for each of the possible anomalies. Increasing Descendants leads to higher accuracy but also higher memory and CPU requirements. The default value is 1.
#' @param anom_add_prob A vector of probabilities with length equal to the dimension of the observations giving the probabilities of additive outliers in each of the components. It defaults to 1/10000.
#' @param anom_inn_prob A vector of probabilities with length equal to the dimension of the hidden state giving the probabilities of innovative outliers in each of the components. It defaults to 1/10000.
#' @param s A numeric giving the shape of the t-distribution to be considered. It defaults to 2. 
#' @param epsilon A positive numeric giving the precision to which the limit of the covariance is to be computed. It defaults to 0.000001.
#' @param horizon_matrix A matrix of 0s and 1s giving the horizon's at which innovative particles are to be resampled. It defaults to a k by q matrix, where k is the number of observations required for observability of the system and q is the dimension of the hidden states.
#' @return An ioaorkf S3 class. 
#' @references \insertRef{fisch2020innovative}{RobKF}
#'
#' @examples
#' 
#' library(RobKF)
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
#' horizon_matrix = matrix(1,nrow = 3 ,ncol = 2)
#' 
#' Particle_List = IOAORKF(Y_list,mu,Sigma_0=NULL,A,C,Sigma_Add,Sigma_Inn,Particles=20,
#'                         horizon_matrix=horizon_matrix)
#' 
#' plot(Particle_List)
#' summary(Particle_List)
#' 
#' 
#' @export
IOAORKF = function(Y,mu_0,Sigma_0=NULL,A,C,Sigma_Add,Sigma_Inn,Particles,Descendants=1,s=2,anom_add_prob=NULL,anom_inn_prob=NULL,epsilon=0.000001,horizon_matrix=NULL)
{
  
  p = nrow(Sigma_Add)
  q = nrow(Sigma_Inn)
  
  if (is.null(p)){
    stop("Sigma_Add must be a matrix.")
  }
  
  if (is.null(q)){
    stop("Sigma_Inn must be a matrix.")
  }
  
  Particles = as.integer(Particles)
  
  Descendants = as.integer(Descendants)
  
  if(Particles < 0.5){
    stop("Particles must be positive!")
  }
  
  if(Descendants < 0.5){
    stop("Descendants must be positive!")
  }
  
  if(Particles < p+q+1){
    warning("Particles should be stricly greater than sum of the dimensions of the observed and hidden states.")
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
  
  if (is.null(anom_add_prob)){
    anom_add_prob = rep(0.0001,p)
  } 
  
  if (length(anom_add_prob) == 1){
    anom_add_prob = rep(anom_add_prob,p)
  }
  
  if (length(anom_add_prob) != p){
    stop("anom_add_prob must be of length equal to the dimensions of the observations")
  }
  
  if (sum(anom_add_prob <= 0) + sum(anom_add_prob >= 1) >0){
    stop("anom_add_prob must be contained in (0,1)" )
  }
  
  if (is.null(anom_inn_prob)){
    anom_inn_prob = rep(0.0001,q)
  }
  
  if (length(anom_inn_prob) == 1){
    anom_inn_prob = rep(anom_inn_prob,nrow(mu_0))
  }
  
  if (length(anom_inn_prob) != q){
    stop("anom_inn_prob must be of length equal to the dimensions of the hidden states")
  }
  
  if (sum(anom_inn_prob <= 0) + sum(anom_inn_prob >= 1) >0){
    stop("anom_inn_prob must be contained in (0,1)" )
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
    
    New_Matrix = New_Matrix %*% A
    
    Full_Matrix = rbind(Full_Matrix,New_Matrix)
    
  }
  
  if (Rank >q){
    stop("The system has to be observable.")
  }
  
  if (epsilon <= 0){
    stop ("epsilon must be greater than 0.")
  }
  
  Final_Sigma = Sigma_Limit(diag(1,nrow = q),C,A,Sigma_Inn,Sigma_Add,epsilon)
  
  if(is.null(Sigma_0)){
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
  
  if (s <= 1){
    stop ("s must be greater than 1.")
  }
  
  
  if (is.null(horizon_matrix)){
    
    horizon_matrix = matrix(1,nrow=Rank,ncol=q)
    
  }
  
  if (nrow(horizon_matrix) < Rank){
    warning("The number of rows of horizon_matrix is less than the Rank of the system. This is not advised.")
  }
  
  if (ncol(horizon_matrix) != q){
    stop("The number of columns of horizon_matrix must be equal to the dimension of the hidden states.")
  }
  
  if (sum(horizon_matrix %in% c(0,1)) < nrow(horizon_matrix)*ncol(horizon_matrix)){
    stop("horizon_matrix must only contain zeros and ones.")
  }
  
  Obs_Matrix = C
  
  for (ii in 1:nrow(horizon_matrix)){
    
    observable = 1 * (colSums(Obs_Matrix != 0) > 0)
    
    horizon_matrix[ii,] = horizon_matrix[ii,] * observable
    
    Obs_Matrix = Obs_Matrix %*% A
  
  } 
  
  if (nrow(horizon_matrix) > length(Y)){
    
    warning("the number of rows horizon_matrix shouldn't exceed the number of observations")
    
    horizon_matrix = matrix(horizon_matrix[1:length(Y),],nrow=length(Y))
    
  }
  
  
  Y_Full_list     = Y
  Particle_Number = Particles
  
  Num_Descendants = Descendants
  Num_Particles   = Particles
  
  prob_add        = anom_add_prob
  
  prob_inn        = anom_inn_prob
  
  horizon = nrow(horizon_matrix)
  
  Number_of_resamples = colSums(horizon_matrix)
  
  if (sum(Number_of_resamples<0.5) > 0){
    stop("horizon_matrix must contain at least one 1 in each column.")
  }
  
  to_sample = list()
  
  for (ii in 1:horizon){
    to_sample[[ii]] = horizon_matrix[ii,]
  }
  
  C_matrix_list = list()
  Sigma_Add_matrix_list = list()
  Sigma_Inn_matrix_list = list()
  
  p = nrow(Sigma_Add)
  q = nrow(Sigma_Inn)
  
  New_matrix  = C
  Full_matrix = C
  New_Inn_Matrix = C %*%  Sigma_Inn %*% t(C)
  
  for (jj in 1:nrow(horizon_matrix)){
    
    C_matrix_list[[jj]] = Full_matrix
    New_matrix          = New_matrix%*%A
    Full_matrix         = rbind(Full_matrix,New_matrix)
    
    Sigma_Add_matrix_list[[jj]] = diag(rep(diag(Sigma_Add),jj),ncol=p*jj)
    Sigma_Inn_matrix_list[[jj]] = New_Inn_Matrix
    
    Even_Newer_New_Matrix = matrix(0,ncol=ncol(New_Inn_Matrix)+nrow(C),nrow=ncol(New_Inn_Matrix)+nrow(C)) 
    Even_Newer_New_Matrix[(1+nrow(C)):nrow(Even_Newer_New_Matrix),(1+nrow(C)):nrow(Even_Newer_New_Matrix)] = New_Inn_Matrix
    
    New_Inn_Matrix = Full_matrix  %*%  Sigma_Inn   %*% t(Full_matrix) + Even_Newer_New_Matrix
    
  }
  
  Y_expanded = list()
  
  for (jj in 1:nrow(horizon_matrix)){
    
    Considered_Y = Y_Full_list[jj:1]
    
    Full_Y_list = list()
    
    New_Y = matrix(0,ncol=1,nrow=0)
    
    for (ii in 1:length(Considered_Y)){
      
      New_Y = rbind(Considered_Y[[ii]],New_Y)
      
      Full_Y_list[[ii]] = New_Y
      
    }
    
    Y_expanded[[jj]] = Full_Y_list
    
  } 
  
  if (length(Y_Full_list) > nrow(horizon_matrix)){
  
    for (jj in (nrow(horizon_matrix)+1):length(Y_Full_list) ){
    
      Considered_Y = Y_Full_list[jj:(jj-nrow(horizon_matrix)+1)]
    
      Full_Y_list = list()
    
      New_Y = matrix(0,ncol=1,nrow=0)
    
      for (ii in 1:length(Considered_Y)){
      
        New_Y = rbind(Considered_Y[[ii]],New_Y)
      
        Full_Y_list[[ii]] = New_Y
      
      }
    
      Y_expanded[[jj]] = Full_Y_list
      
    }
    
  }
  
  sigma_tilde = get_sigma_tilde(Final_Sigma,C,A,Sigma_Add,Sigma_Inn)
  
  sigma_hat   = get_sigma_hat(Final_Sigma,C_matrix_list,A,Sigma_Add_matrix_list,Sigma_Inn,Sigma_Inn_matrix_list,horizon_matrix)
  
  #### forget abut the stuff below:
  
  precision = solve( C %*% ( A %*%  Sigma_0  %*% t(A)  +  Sigma_Inn ) %*%  t(C) + Sigma_Add )
  
  
  tryCatch(
  {
    Pre_Out = Robust_filter(Y_expanded, C_matrix_list, Sigma_Add_matrix_list, Sigma_Inn_matrix_list, A, Sigma_Inn, Sigma_Add, s, Num_Descendants, Num_Particles, to_sample, Number_of_resamples, sigma_tilde, sigma_hat, mu_0, Sigma_0, horizon, prob_inn, prob_add, Particle_Number, Y_Full_list)
  },error = function(e) {print(e$message);stop();})
  
  Transform_particle = function(x){
    
    out = list()
    
    out[["mu"]]    = x[[1]]
    out[["Sigma"]] = x[[2]]
    
    out[["strength"]]   = as.numeric(x[[8]])
    if (x[[3]] < 0.5){
      out[["which_type"]] = "None"
    } else {
      if (x[[3]] < 1.5){
        out[["which_type"]] = "W"
      } else {
        out[["which_type"]] = "V"
      }
    }
    
    out[["component"]]    = as.numeric(x[[4]])+1
    out[["horizon"]]      = as.numeric(x[[7]])
    out[["ancestor_id"]]  = as.numeric(x[[6]])+1 
    out[["id"]]           = as.numeric(x[[5]])+1
      
    return(out)
  }
  
  Transform_list = function(x){
    lapply(x,Transform_particle)
  }
  
  Particles = lapply(Pre_Out, Transform_list)
  
  output = list()
  
  output[["particles"]] = Particles
  output[["p"]]         = p
  output[["q"]]         = q
  
  output[["Y"]] = Y
  
  output[["horizon"]] = horizon
  
  return(structure(output,class="ioaorkf"))  
  
}