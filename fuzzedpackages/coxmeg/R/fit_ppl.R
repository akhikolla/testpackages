
#' Estimate HRs using PPL given a known variance component (tau)
#'
#' \code{fit_ppl} returns estimates of HRs and their p-values given a known variance component (tau).
#' 
#' @section About \code{type}:
#' 'bd' is used for a block-diagonal relatedness matrix, or a sparse matrix the inverse of which is also sparse. 'sparse' is used for a general sparse relatedness matrix the inverse of which is not sparse. 
#' @section About \code{solver}:
#' When \code{solver=1,3}/\code{solver=2}, Cholesky decompositon/PCG is used to solve the linear system. When \code{solver=3}, the solve function in the Matrix package is used, and when \code{solver=1}, it uses RcppEigen:LDLT to solve linear systems. 
#' 
#' @param tau A positive scalar. A variance component given by the user. Default is 0.5.
#' @param X A matrix of the preidctors. Can be quantitative or binary values. Categorical variables need to be converted to dummy variables. Each row is a sample, and the predictors are columns.
#' @param outcome A matrix contains time (first column) and status (second column). The status is a binary variable (1 for failure / 0 for censored).
#' @param corr A relatedness matrix. Can be a matrix or a 'dgCMatrix' class in the Matrix package. Must be symmetric positive definite or symmetric positive semidefinite.
#' @param type A string indicating the sparsity structure of the relatedness matrix. Should be 'bd' (block diagonal), 'sparse', or 'dense'. See details.
#' @param FID An optional string vector of family ID. If provided, the data will be reordered according to the family ID.
#' @param eps An optional positive value indicating the tolerance in the optimization algorithm. Default is 1e-6.
#' @param spd An optional logical value indicating whether the relatedness matrix is symmetric positive definite. Default is TRUE. 
#' @param solver An optional bianry value that can be either 1 (Cholesky Decomposition using RcppEigen), 2 (PCG) or 3 (Cholesky Decomposition using Matrix). Default is NULL, which lets the function select a solver. See details.
#' @param verbose An optional logical value indicating whether to print additional messages. Default is TRUE.
#' @param order An optional integer value starting from 0. Only valid when dense=FALSE. It specifies the order of approximation used in the inexact newton method. Default is 1.
#' @return beta: The estimated coefficient for each predictor in X.
#' @return HR: The estimated HR for each predictor in X.
#' @return sd_beta: The estimated standard error of beta.
#' @return p: The p-value.
#' @return iter: The number of iterations until convergence.
#' @return ppl: The PPL when the convergence is reached.
#' @keywords Cox mixed-effects model
#' @export fit_ppl
#' @examples
#' library(Matrix)
#' library(MASS)
#' library(coxmeg)
#' 
#' ## simulate a block-diagonal relatedness matrix
#' tau_var <- 0.2
#' n_f <- 100
#' mat_list <- list()
#' size <- rep(10,n_f)
#' offd <- 0.5
#' for(i in 1:n_f)
#' {
#'   mat_list[[i]] <- matrix(offd,size[i],size[i])
#'   diag(mat_list[[i]]) <- 1
#' }
#' sigma <- as.matrix(bdiag(mat_list))
#' n <- nrow(sigma)
#' 
#' ## simulate random effexts and outcomes
#' x <- mvrnorm(1, rep(0,n), tau_var*sigma)
#' myrates <- exp(x-1)
#' y <- rexp(n, rate = myrates)
#' cen <- rexp(n, rate = 0.02 )
#' ycen <- pmin(y, cen)
#' outcome <- cbind(ycen,as.numeric(y <= cen))
#' 
#' ## fit the ppl
#' re = fit_ppl(x,outcome,sigma,type='bd',tau=0.5,order=1)
#' re

fit_ppl <- function(X,outcome,corr,type,tau=0.5,FID=NULL,eps=1e-06,order=1,solver=NULL,spd=TRUE,verbose=TRUE){

  if(eps<0)
  {eps <- 1e-06}
  
  if(!(type %in% c('bd','sparse','dense')))
  {stop("The type argument should be 'bd', 'sparse' or 'dense'.")}
  
  ## family structure
  if(is.null(FID)==FALSE)
  {
    ord <- order(FID)
    FID <- as.character(FID[ord])
    X <- as.matrix(X[ord,,drop = FALSE])
    outcome <- as.matrix(outcome[ord,,drop = FALSE])
    corr <- corr[ord,ord,drop = FALSE]
  }else{
    X <- as.matrix(X)
    outcome <- as.matrix(outcome)
  }
  
  min_d <- min(outcome[which(outcome[,2]==1),1])
  rem <- which((outcome[,2]==0)&(outcome[,1]<min_d))
  if(length(rem)>0)
  {
    outcome <- outcome[-rem, ,drop = FALSE]
    X <- as.matrix(X[-rem,,drop = FALSE])
    corr <- corr[-rem,-rem,drop = FALSE]
  }
  
  if(verbose==TRUE)
  {message(paste0('Remove ', length(rem), ' subjects censored before the first failure.'))}
  
  x_sd = which(as.vector(apply(X,2,sd))>0)
  x_ind = length(x_sd)
  if(x_ind==0)
  {stop("The predictors are all constants after the removal of subjects.")}else{
    k <- ncol(X)
    if(x_ind<k)
    {
      warning(paste0(k-x_ind," predictor(s) is/are removed because they are all constants after the removal of subjects."))
      X = X[,x_sd,drop=FALSE]
      k <- ncol(X)
    }
  }
  
  n <- nrow(outcome)
  if(min(outcome[,2] %in% c(0,1))<1)
  {stop("The status should be either 0 (censored) or 1 (failure).")}
  u <- rep(0,n)
  beta <- rep(0,k)
  
  d_v <- outcome[,2]
  
  ## risk set matrix
  ind <- order(outcome[,1])
  ind <- as.matrix(cbind(ind,order(ind)))
  rk <- rank(outcome[ind[,1],1],ties.method='min')
  n1 <- sum(d_v>0)
  
  rs <- rs_sum(rk-1,d_v[ind[,1]])
  if(spd==FALSE)
  {
    rk_cor = matrix.rank(as.matrix(corr),method='chol')
    spsd = FALSE
    if(rk_cor<n)
    {spsd = TRUE}
    if(verbose==TRUE)
    {message(paste0('The sample size included is ',n,'. The rank of the relatedness matrix is ', rk_cor))}
    
  }else{
    spsd = FALSE
    rk_cor = n
    if(verbose==TRUE)
    {message(paste0('The sample size included is ',n,'.'))}
  }
  
  nz <- nnzero(corr)
  if( nz > ((as.double(n)^2)/2) )
  {type <- 'dense'}
  inv = NULL
  
  eigen = TRUE
  if(type=='dense')
  {
    if(verbose==TRUE)
    {message('The relatedness matrix is treated as dense.')}
    corr = as.matrix(corr)
    if(spsd==FALSE)
    {
      corr = chol(corr)
      corr = as.matrix(chol2inv(corr))
    }else{
      ei = eigen(corr)
      ei$values[ei$values<1e-10] = 1e-6
      corr = ei$vectors%*%diag(1/ei$values)%*%t(ei$vectors)
      # corr <- ginv(corr)
      rk_cor = n
      spsd = FALSE
    }
    inv <- TRUE
    sigma_i_s = corr
    corr <- s_d <- NULL
    si_d <- as.vector(diag(sigma_i_s))
    
    if(is.null(solver))
    {solver = 2}else{
      if(solver==3)
      {solver = 1}
    }
    
  }else{
    
    if(verbose==TRUE)
    {message('The relatedness matrix is treated as sparse.')}
    
    corr <- as(corr, 'dgCMatrix')
    si_d = s_d = NULL
    
    if(spsd==FALSE)
    {
      if(type=='bd')
      {
        sigma_i_s <- Matrix::chol2inv(Matrix::chol(corr))
        inv = TRUE
        si_d <- as.vector(Matrix::diag(sigma_i_s))
      }else{
        sigma_i_s = NULL
        inv = FALSE
        s_d <- as.vector(Matrix::diag(corr))
        
      }
    }else{
      sigma_i_s = eigen(corr)
      if(min(sigma_i_s$values) < -1e-10)
      {
        stop("The relatedness matrix has negative eigenvalues.")
      }
      # sigma_i_s = sigma_i_s$vectors%*%(c(1/sigma_i_s$values[1:rk_cor],rep(0,n-rk_cor))*t(sigma_i_s$vectors))
      sigma_i_s$values[sigma_i_s$values<1e-10] = 1e-6
      sigma_i_s = sigma_i_s$vectors%*%diag(1/sigma_i_s$values)%*%t(sigma_i_s$vectors)
      rk_cor = n
      spsd = FALSE
      
      inv = TRUE
      si_d <- as.vector(Matrix::diag(sigma_i_s))
    }
    
    if(is.null(solver))
    {
      if(type=='bd')
      {
        solver = 1
        if(n>5e4)
        {
          eigen = FALSE
        }
      }else{solver = 2}
    }else{
      if(solver==3)
      {
        eigen = FALSE
        solver = 1
      }
    }
    
    if(inv==TRUE)
    {
      sigma_i_s <- as(sigma_i_s,'dgCMatrix')
      if(eigen==FALSE)
      {
        sigma_i_s = Matrix::forceSymmetric(sigma_i_s)
      }
      corr <- s_d <- NULL
    }
  }
  
  if(verbose==TRUE)
  {
    if(inv==TRUE)
    {message('The relatedness matrix is inverted.')}
    
    if(type=='dense')
    {
      switch(
        solver,
        '1' = message('Solver: solve (base).'),
        '2' = message('Solver: PCG (RcppEigen:dense).')
      )
    }else{
      switch(
        solver,
        '1' = message('Solver: Cholesky decomposition (RcppEigen=',eigen,').'),
        '2' = message('Solver: PCG (RcppEigen:sparse).')
      )
    }
  }
  
  if(type=='dense')
  {
    res <- irls_ex(beta, u, tau, si_d, sigma_i_s, X, eps, d_v, ind, rs$rs_rs, rs$rs_cs,rs$rs_cs_p,det=FALSE,detap='exact',solver=solver)
  }else{
    res <- irls_fast_ap(beta, u, tau, si_d, sigma_i_s, X, eps, d_v, ind, rs$rs_rs, rs$rs_cs,rs$rs_cs_p,order,det=FALSE,detap='exact',sigma_s=corr,s_d=s_d,eigen=eigen,solver=solver)
  }
  
  res_beta = as.vector(res$beta)
  res_var = diag(as.matrix(res$v11))
  p = pchisq(res_beta^2/res_var,1,lower.tail=FALSE)
  
  re = list(beta=res_beta,HR=exp(res_beta),sd_beta=sqrt(res_var),p=as.vector(p),iter=res$iter,ppl=res$ll)
  return(re)
}

