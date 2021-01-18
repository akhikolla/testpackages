

irls_ex <- function(beta, u, tau,si_d, sigma_i_s, X, eps=1e-6, d_v, ind, rs_rs, rs_cs, rs_cs_p,det=FALSE,detap='slq',solver=1,rad=NULL)
{
  n <- length(u)
  n_c <- length(beta)
  dim_v <- n_c + n
  brc <- (n_c+1):dim_v
  n1_ind <- which(d_v>0)
  tol = eps*1e-3
  maxiter = 200
  
  ## rs_rs, rs_cs only for c++
  rs_rs = rs_rs - 1
  rs_cs = rs_cs - 1
  
  u_new <- u
  beta_new <- vi11 <- NULL
  if(n_c>0)
  {
    beta_new <- beta
    eta_v <- X%*%beta_new+u_new
  }else{
    eta_v <- u_new
  }
  
  # sigma_i_s = sigma_i_s/tau
  
  w_v <- as.vector(exp(eta_v))
  ## (-1) for cpp
  s <- as.vector(cswei(w_v,rs_rs,ind-1,1))
  
  loglik <- 0
  newloglik <- 0
  siu <- as.vector(sigma_i_s%*%u_new)/tau
  newloglik <- sum(eta_v[n1_ind]) - sum(log(s)[n1_ind]) - 0.5*t(u_new)%*%siu
  lik_dif <- 0
  iter <- 0
  eps_s <- 0
  
  v <- matrix(NA,dim_v,dim_v)
  
  while(((lik_dif>eps_s)||(iter<1))&&(iter<maxiter))
  {
    damp = 1
    loglik <- newloglik
    a_v <- d_v/s
    
    bw_v <- w_v*as.vector(cswei(a_v,rs_cs,ind-1,0))
    deriv <- d_v - bw_v
    deriv_full <- c(as.vector(t(X)%*%deriv),deriv)-c(rep(0,n_c),siu)
    
    a_v_p <- a_v[ind[,1]]
    a_v_2 <- as.vector(a_v_p*a_v_p)
    a_v_p <- a_v_p[a_v_p>0]
    
    v[brc,brc] = sigma_i_s/tau - wma_cp(w_v,rs_cs_p-1,ind-1,a_v_p)
    diag(v[brc,brc]) = diag(v[brc,brc]) + bw_v
    
    if(n_c>0)
    {
      v[1:n_c,(n_c+1):dim_v] <- t(bw_v*X - csqei(w_v,X,rs_rs,rs_cs,ind-1,a_v_2))
      # v[1:n_c,(n_c+1):dim_v] <- t(X)%*%v[brc,brc]
      v[1:n_c,1:n_c] <- v[1:n_c,(n_c+1):dim_v]%*%X
      v[(n_c+1):dim_v,1:n_c] <- t(v[1:n_c,(n_c+1):dim_v])
    }
    if(solver<2)
    {
      new <- solve(v,deriv_full)
    }else{
      new = pcg_dense(v, as.matrix(deriv_full),tol)
    }
    
    u_new <- u_new + new[(n_c+1):dim_v]
    if(n_c>0)
    {
      beta_new <- beta_new + new[1:n_c]
      eta_v <- X%*%beta_new+u_new
    }else{
      eta_v <- u_new
    }
    
    w_v <- as.vector(exp(eta_v))
    s <- as.vector(cswei(w_v,rs_rs,ind-1,1))
    siu <- as.vector(sigma_i_s%*%u_new)/tau
    newloglik <- sum(eta_v[n1_ind]) - sum(log(s)[n1_ind]) - 0.5*t(u_new)%*%siu
    eps_s = eps*(-1)*loglik
    lik_dif <- as.numeric(newloglik - loglik)
    
    while(lik_dif<(-eps_s))
    {
      damp = damp/2
      if(damp<1e-2)
      {
        warning(paste0("The optimization of PPL may not converge."))
        lik_dif = 0
        break
      }
      new_d = damp*new
      u_new <- u_new - new_d[(n_c+1):dim_v]
      if(n_c>0)
      {
        beta_new <- beta_new - new_d[1:n_c]
        eta_v <- X%*%beta_new+u_new
      }else{
        eta_v <- u_new
      }
      
      w_v <- as.vector(exp(eta_v))
      s <- as.vector(cswei(w_v,rs_rs,ind-1,1))
      siu <- as.vector(sigma_i_s%*%u_new)/tau
      newloglik <- sum(eta_v[n1_ind]) - sum(log(s)[n1_ind]) - 0.5*t(u_new)%*%siu
      lik_dif <- as.numeric(newloglik - loglik)
    }
    
    iter <- iter + 1
  }
  
  if(iter==maxiter)
  {warning(paste0("The number of iterations reaches maxiter (", maxiter,"). The optimization of PPL has likely not converged."))}
  
  logdet <- 0
  
  if(det==FALSE)
  {
    if(n_c>0)
    {
      a_v <- d_v/s
      bw_v <- w_v*as.vector(cswei(a_v,rs_cs,ind-1,0))
      a_v_p <- a_v[ind[,1]]
      a_v_2 <- as.vector(a_v_p*a_v_p)
      a_v_p <- a_v_p[a_v_p>0]
      # hx = as.matrix(v[brc,brc]%*%X)
      hx = bw_v*X - csqei(w_v,X,rs_rs,rs_cs,ind-1,a_v_2)
      # v[brc,brc] = diag(bw_v) - wma_cp(w_v,rs_cs_p-1,ind-1,a_v_p)+sigma_i_s
      v[brc,brc] = sigma_i_s/tau - wma_cp(w_v,rs_cs_p-1,ind-1,a_v_p)
      diag(v[brc,brc]) = diag(v[brc,brc]) + bw_v
      vi11 = pcg_dense(v[brc,brc],hx,tol)
      vi11 <- solve(t(hx)%*%(X - vi11))
    }
  }else{
    a_v <- d_v/s
    bw_v <- w_v*as.vector(cswei(a_v,rs_cs,ind-1,0))
    
    a_v_p <- a_v[ind[,1]]
    a_v_p <- a_v_p[a_v_p>0]
    
    if(detap=='slq')
    {
      v[brc,brc] = sigma_i_s/tau - wma_cp(w_v,rs_cs_p-1,ind-1,a_v_p)
      diag(v[brc,brc]) = diag(v[brc,brc]) + bw_v
      
      v[brc,brc] = v[brc,brc]*tau
      logdet = logdet_lanczos(v[brc,brc], rad, 8) - n*log(tau)
    }else{
      logdet <- logdeth(as(sigma_i_s,'dgCMatrix'),si_d,bw_v, w_v,rs_cs_p-1,ind-1,a_v_p,tau,1,0)
    }
    
  }
  
  return(list(beta=beta_new,u=u_new,v11=vi11,iter=iter,ll=newloglik,logdet=logdet))
}
