
#' Fit a Cox mixed-effects model
#'
#' \code{coxmeg} returns estimates of the variance component, the HRs and p-values for the predictors.
#' 
#' @section About \code{type}:
#' 'bd' is used for a block-diagonal relatedness matrix, or a sparse matrix the inverse of which is also sparse. 'sparse' is used for a general sparse relatedness matrix the inverse of which is not sparse. 
#' @section About \code{spd}:
#' When \code{spd=TRUE}, the relatedness matrix is treated as SPD. If the matrix is SPSD or not sure, use \code{spd=FALSE}.
#' @section About \code{solver}:
#' When \code{solver=1,3}/\code{solver=2}, Cholesky decompositon/PCG is used to solve the linear system. When \code{solver=3}, the solve function in the Matrix package is used, and when \code{solver=1}, it uses RcppEigen:LDLT to solve linear systems. When \code{type='dense'}, it is recommended to set \code{solver=2} to have better computational performance. 
#' @section About \code{detap}:
#' When \code{detap='exact'}, the exact log-determinant is computed for estimating the variance component. Specifying \code{detap='diagonal'} uses diagonal approximation, and is only effective for a sparse relatedness matrix. Specifying \code{detap='slq'} uses stochastic lanczos quadrature approximation.   
#' 
#' @param outcome A matrix contains time (first column) and status (second column). The status is a binary variable (1 for events / 0 for censored).
#' @param corr A relatedness matrix. Can be a matrix or a 'dgCMatrix' class in the Matrix package. Must be symmetric positive definite or symmetric positive semidefinite.
#' @param type A string indicating the sparsity structure of the relatedness matrix. Should be 'bd' (block diagonal), 'sparse', or 'dense'. See details.
#' @param X An optional matrix of the preidctors with fixed effects. Can be quantitative or binary values. Categorical variables need to be converted to dummy variables. Each row is a sample, and the predictors are columns. 
#' @param FID An optional string vector of family ID. If provided, the data will be reordered according to the family ID.
#' @param eps An optional positive scalar indicating the tolerance in the optimization algorithm. Default is 1e-6.
#' @param min_tau An optional positive scalar indicating the lower bound in the optimization algorithm for the variance component \code{tau}. Default is 1e-4.
#' @param max_tau An optional positive scalar indicating the upper bound in the optimization algorithm for the variance component \code{tau} Default is 5.
#' @param opt An optional logical scalar for the Optimization algorithm for tau. Can have the following values: 'bobyqa', 'Brent' or 'NM'. Default is 'bobyqa'.
#' @param spd An optional logical value indicating whether the relatedness matrix is symmetric positive definite. Default is TRUE. See details.
#' @param detap An optional string indicating whether to use approximation for log-determinant. Can be 'exact', 'diagonal' or 'slq'. Default is NULL, which lets the function select a method based on 'type' and other information. See details.
#' @param solver An optional bianry value that can be either 1 (Cholesky Decomposition using RcppEigen), 2 (PCG) or 3 (Cholesky Decomposition using Matrix). Default is NULL, which lets the function select a solver. See details.
#' @param order An optional integer starting from 0. Only valid when \code{dense=FALSE}. It specifies the order of approximation used in the inexact newton method. Default is 1.
#' @param verbose An optional logical scalar indicating whether to print additional messages. Default is TRUE.
#' @param mc An optional integer scalar specifying the number of Monte Carlo samples used for approximating the log-determinant. Only valid when \code{dense=TRUE} and \code{detap='slq'}. Default is 100.
#' @return beta: The estimated coefficient for each predictor in X.
#' @return HR: The estimated HR for each predictor in X.
#' @return sd_beta: The estimated standard error of beta.
#' @return p: The p-value.
#' @return iter: The number of iterations until convergence.
#' @return tau: The estimated variance component.
#' @return int_ll: The marginal likelihood (-2*log(lik)) of tau evaluated at the estimate of tau.
#' @return rank: The rank of the relatedness matrix.
#' @return nsam: Actual sample size.
#' @keywords Cox mixed-effects model
#' @export coxmeg
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
#' ## simulate random effects and outcomes
#' x <- mvrnorm(1, rep(0,n), tau_var*sigma)
#' myrates <- exp(x-1)
#' y <- rexp(n, rate = myrates)
#' cen <- rexp(n, rate = 0.02 )
#' ycen <- pmin(y, cen)
#' outcome <- cbind(ycen,as.numeric(y <= cen))
#' 
#' ## fit a Cox mixed-effects model
#' re = coxmeg(outcome,sigma,type='bd',detap='diagonal',order=1)
#' re

coxmeg <- function(outcome,corr,type,X=NULL,FID=NULL,eps=1e-6, min_tau=1e-04,max_tau=5,order=1,detap=NULL,opt='bobyqa',solver=NULL,spd=TRUE,verbose=TRUE, mc=100){
  
  if(eps<0)
  {eps <- 1e-6}
  
  if(!(type %in% c('bd','sparse','dense')))
  {stop("The type argument should be 'bd', 'sparse' or 'dense'.")}
  
  if(!is.null(detap))
  {
    if(!(detap %in% c('exact','slq','diagonal')))
    {stop("The detap argument should be 'exact', 'diagonal' or 'slq'.")}
  }
    
  ## family structure
  if(is.null(FID)==FALSE)
  {
    ord <- order(FID)
    FID <- as.character(FID[ord])
    X <- as.matrix(X[ord,])
    outcome <- as.matrix(outcome[ord,])
    corr <- corr[ord,ord]
  }else{
    if(is.null(X)==FALSE)
    {
      X <- as.matrix(X)
    }
    outcome <- as.matrix(outcome)
  }
  
  min_d <- min(outcome[which(outcome[,2]==1),1])
  rem <- which((outcome[,2]==0)&(outcome[,1]<min_d))
  if(length(rem)>0)
  {
    outcome <- outcome[-rem,,drop = FALSE]
    if(is.null(X)==FALSE)
    {
      X <- as.matrix(X[-rem,,drop = FALSE])
    }
    corr <- corr[-rem,-rem,drop = FALSE]
  }
  
  if(verbose==TRUE)
  {message(paste0('Remove ', length(rem), ' subjects censored before the first failure.'))}
  
  n <- nrow(outcome)
  if(min(outcome[,2] %in% c(0,1))<1)
  {stop("The status should be either 0 (censored) or 1 (failure).")}
  
  u <- rep(0,n)
  if(is.null(X)==FALSE)
  {
    x_sd = which(as.vector(apply(X,2,sd))>0)
    x_ind = length(x_sd)
    if(x_ind==0)
    {
      warning("The predictors are all constants after the removal of subjects.")
      k <- 0
      beta <- numeric(0)
      X <- matrix(0,0,0)
    }else{
      k <- ncol(X)
      if(x_ind<k)
      {
        warning(paste0(k-x_ind," predictor(s) is/are removed because they are all constants after the removal of subjects."))
        X = X[,x_sd,drop=FALSE]
        k <- ncol(X)
      }
      beta <- rep(0,k)
    }
  }else{
    k <- 0
    beta <- numeric(0)
    X <- matrix(0,0,0)
  }
  
  tau <- 0.5
  
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
    {message(paste0('There is/are ',k,' predictors. The sample size included is ',n,'. The rank of the relatedness matrix is ', rk_cor))}
    
  }else{
    spsd = FALSE
    rk_cor = n
    if(verbose==TRUE)
    {message(paste0('There is/are ',k,' predictors. The sample size included is ',n,'.'))}
  }
  
  nz <- nnzero(corr)
  if( nz > ((as.double(n)^2)/2) )
  {type <- 'dense'}
  
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
    corr = s_d = NULL
    
    si_d <- as.vector(diag(sigma_i_s))
    
    if(is.null(solver))
    {solver = 2}else{
      if(solver==3)
      {solver = 1}
    }
    if(is.null(detap))
    {
      if(n<3000)
      {detap = 'exact'}else{
        detap = 'slq'
      }
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
        if(is.null(detap))
        {detap = 'diagonal'}
      }else{
        sigma_i_s = NULL
        inv = FALSE
        s_d <- as.vector(Matrix::diag(corr))
        if(is.null(detap))
        {
          detap = 'slq'
        }else{
          if(detap=='exact')
          {
            detap = 'slq'
            if(verbose==TRUE)
            {message("detap=exact is not supported under this setting. The detap argument is changed to 'slq'.")}
          }
        }
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
      if(is.null(detap))
      {
        if(type=='bd')
        {detap = 'diagonal'}else{
          detap = 'slq'
        }
      }
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
    
    message("The method for computing the determinant is '", detap, "'.")
    
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
  
  rad = NULL
  if(detap=='slq')
  {
    rad = rbinom(n*mc,1,0.5)
    rad[rad==0] = -1
    rad = matrix(rad,n,mc)/sqrt(n)
  }
  
  marg_ll = 0
  new_t = switch(
    opt,
    'bobyqa' = bobyqa(tau, mll, type=type, beta=beta,u=u,si_d=si_d,sigma_i_s=sigma_i_s,X=X, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order,det=TRUE,detap=detap,inv=inv,sigma_s=corr,s_d=s_d,eps=eps,lower=min_tau,upper=max_tau,eigen=eigen,solver=solver,rad=rad),
    'Brent' = optim(tau, mll, type=type, beta=beta,u=u,si_d=si_d,sigma_i_s=sigma_i_s,X=X, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order,det=TRUE,detap=detap,inv=inv,sigma_s=corr,s_d=s_d,eps=eps,lower=min_tau,upper=max_tau,method='Brent',eigen=eigen,solver=solver,rad=rad),
    'NM' = optim(tau, mll, type=type, beta=beta,u=u,si_d=si_d,sigma_i_s=sigma_i_s,X=X, d_v=d_v, ind=ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,rk_cor=rk_cor,order=order,det=TRUE,detap=detap,inv=inv,sigma_s=corr,s_d=s_d,eps=eps,method='Nelder-Mead',eigen=eigen,solver=solver,rad=rad),
    stop("The argument opt should be bobyqa, Brent or NM.")
  )
  marg_ll = new_t$value
  if(opt=='bobyqa')
  {iter <- new_t$iter}else{
    iter <- new_t$counts
  }
  
  tau_e <- new_t$par
  if(tau_e==min_tau)
  {warning(paste0("The estimated variance component equals the lower bound (", min_tau, "), probably suggesting no random effects."))}
  
  if(type=='dense')
  {
    re <- irls_ex(beta, u, tau_e, si_d, sigma_i_s, X, eps, d_v, ind, rs$rs_rs, rs$rs_cs,rs$rs_cs_p,det=FALSE,detap=detap,solver=solver)
  }else{
    re <- irls_fast_ap(beta, u, tau_e, si_d, sigma_i_s, X, eps, d_v, ind, rs_rs=rs$rs_rs, rs_cs=rs$rs_cs,rs_cs_p=rs$rs_cs_p,order,det=FALSE,detap=detap,sigma_s=corr,s_d=s_d,eigen=eigen,solver=solver)
  }
  
  if(k>0)
  {
    res_beta = as.vector(re$beta)
    res_var = diag(as.matrix(re$v11))
    HR = exp(res_beta)
    sdb = sqrt(res_var)
    p = as.vector(pchisq(res_beta^2/res_var,1,lower.tail=FALSE))
  }else{
    res_beta=HR=sdb=p=NULL
  }
  
  res <- list(beta=res_beta,HR=HR,sd_beta=sdb,p=p,tau=tau_e,iter=iter,rank=rk_cor,nsam=n,int_ll=marg_ll)
  return(res)
}

