

irls_fast_ap <- function(beta, u, tau,si_d, sigma_i_s, X, eps=1e-6, d_v, ind, rs_rs, rs_cs, rs_cs_p, order=1,det=FALSE,detap='diagonal',sigma_s=NULL,s_d=NULL,eigen=TRUE,solver=1,rad=NULL)
{
  n <- length(u)
  n_c <- length(beta)
  dim_v <- n_c + n
  n1_ind <- which(d_v>0)
  inv <- TRUE
  if(is.null(s_d)==FALSE)
  {inv <- FALSE}
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
  w_v <- as.vector(exp(eta_v))
  ## (-1) for cpp
  s <- as.vector(cswei(w_v,rs_rs,ind-1,1))
  loglik <- newloglik <- 0
  
  if(inv==FALSE)
  {
    siu = as.vector(pcg_sparse(sigma_s,as.matrix(u_new),eps*1e-3)) 
  }else{
    siu = as.vector(sigma_i_s%*%u_new)
  }
  newloglik <- sum(eta_v[n1_ind]) - sum(log(s)[n1_ind]) - 0.5*t(u_new)%*%siu/tau
  lik_dif <- iter <- eps_s <- logdet <- 0
  maxiter = 200
  
  if(eigen==FALSE)
  {
    if(inv==TRUE)
    {
      sigma_i_s_bk <- as(sigma_i_s,'symmetricMatrix')
    }else{
      sigma_i_s_bk <- as(sigma_s,'symmetricMatrix')
    }
  }
  
  while(((lik_dif>eps_s)||(iter<1))&&(iter<maxiter))
  {
    damp = 1
    loglik <- newloglik
    a_v <- d_v/s
    bw_v <- w_v*as.vector(cswei(a_v,rs_cs,ind-1,0))
    deriv <- d_v - bw_v
    if(n_c>0)
    {
      deriv_full <- c(as.vector(t(X)%*%deriv),deriv)-c(rep(0,n_c),siu/tau)
    }else{
      deriv_full <- deriv - siu/tau
    }
    der_t <- deriv_full[(n_c+1):dim_v]	
    
    if(eigen==TRUE)
    {
      if(inv==TRUE)
      {
        re <- invsph(sigma_i_s, der_t, si_d,w_v,X,rs_rs,rs_cs,ind-1,a_v[ind[,1]], bw_v, order, 1,tau,solver)
      }else{
        re <- invsph(sigma_s, der_t, s_d,w_v,X,rs_rs,rs_cs,ind-1,a_v[ind[,1]], bw_v, order, 0,tau,solver)
      }
      wb_sig_i_hx_der <- re$wb_sig_i_hx_der
      xh <- re$xh
    }else{
      a_v_2 <- (as.vector(a_v)*as.vector(a_v))[ind[,1]]
      wb_sig_i_hx_der <- der_t
      if(n_c>0){
        xh <- t(bw_v*X - csqei(w_v,X,rs_rs,rs_cs,ind-1,a_v_2))
        wb_sig_i_hx_der <- cbind(wb_sig_i_hx_der, t(xh))
      }
      if(inv==TRUE)
      {
        diag(sigma_i_s_bk) <- si_d + bw_v*tau
        wb_sig_i_hx_der <- tau*Matrix::solve(sigma_i_s_bk,wb_sig_i_hx_der)
      }else{
        bw_i <- 1/bw_v
        diag(sigma_i_s_bk) <- s_d + bw_i/tau
        wb_sig_i_hx_der <- bw_i*Matrix::solve(sigma_i_s_bk,sigma_s%*%wb_sig_i_hx_der)
      }
      if(order>0)
      {
        wb_sig_i_hx_der_2 <- wb_sig_i_hx_der
        if(inv==TRUE)
        {
          for(i in 1:order)
          {
            wb_sig_i_hx_der_2 <- csqei(w_v,as.matrix(wb_sig_i_hx_der_2),rs_rs,rs_cs,ind-1,a_v_2)
            wb_sig_i_hx_der_2 <- tau*Matrix::solve(sigma_i_s_bk,wb_sig_i_hx_der_2)
            wb_sig_i_hx_der <- wb_sig_i_hx_der+wb_sig_i_hx_der_2
          }
        }else{
          for(i in 1:order)
          {
            wb_sig_i_hx_der_2 <- csqei(w_v,as.matrix(wb_sig_i_hx_der_2),rs_rs,rs_cs,ind-1,a_v_2)
            wb_sig_i_hx_der_2 <- bw_i*Matrix::solve(sigma_i_s_bk,sigma_s%*%wb_sig_i_hx_der_2)
            wb_sig_i_hx_der <- wb_sig_i_hx_der+wb_sig_i_hx_der_2
          }
        }
      }
    }
    
    if(n_c>0)
    {
      wb_sig_i_hx_der_p1 <- wb_sig_i_hx_der[,(2:(n_c+1))]
      vi11 <- solve(xh%*%(X - wb_sig_i_hx_der_p1))
    
      vi12 <- (-1)*vi11%*%Matrix::t(wb_sig_i_hx_der_p1)
      vi12der = vi12%*%der_t
      new_x <- vi11%*%deriv_full[1:n_c] + vi12der
      new_z <- Matrix::t(vi12)%*%deriv_full[1:n_c] + wb_sig_i_hx_der[,1] - wb_sig_i_hx_der_p1%*%vi12der
      beta_new <- beta_new + new_x
      u_new <- u_new + new_z
      eta_v <- X%*%beta_new+u_new
    }else{
      new_z <- wb_sig_i_hx_der[,1]
      u_new <- u_new + new_z
      eta_v <- u_new
    }
    
    w_v <- as.vector(exp(eta_v))
    s <- as.vector(cswei(w_v,rs_rs,ind-1,1))
    
    # siu <- as.vector(sigma_i_s%*%u_new)
    if(inv==FALSE)
    {
      siu = as.vector(pcg_sparse(sigma_s,as.matrix(u_new),eps*1e-3)) 
    }else{
      siu = as.vector(sigma_i_s%*%u_new)
    }
    usiu <- Matrix::t(u_new)%*%siu
    newloglik <- sum(eta_v[n1_ind]) - sum(log(s)[n1_ind]) - 0.5*usiu/tau
    lik_dif <- as.numeric(newloglik) - as.numeric(loglik)
    eps_s = eps*(-1)*as.numeric(loglik)
    
    while(lik_dif<(-eps_s))
    {
      damp = damp/2
      if(damp<1e-2)
      {
        warning(paste0("The optimization of PPL may not converge."))
        lik_dif = 0
        break
      }
      
      u_new <- u_new - damp*new_z
      if(n_c>0)
      {
        beta_new <- beta_new - damp*new_x
        eta_v <- X%*%beta_new+u_new
      }else{
        eta_v <- u_new
      }
      
      w_v <- as.vector(exp(eta_v))
      s <- as.vector(cswei(w_v,rs_rs,ind-1,1))
      if(inv==FALSE)
      {
        siu = as.vector(pcg_sparse(sigma_s,as.matrix(u_new),eps*1e-3)) 
      }else{
        siu = as.vector(sigma_i_s%*%u_new)
      }
      usiu <- Matrix::t(u_new)%*%siu
      newloglik <- sum(eta_v[n1_ind]) - sum(log(s)[n1_ind]) - 0.5*usiu/tau
      lik_dif <- as.numeric(newloglik - loglik)
    }
    
    iter <- iter + 1
  }
  
  if(iter==maxiter)
  {warning(paste0("The number of iterations reaches maxiter (", maxiter,"). The optimization of PPL has likely not converged."))}
  
  if(det==FALSE)
  {
    if(n_c>0)
    {
      a_v <- d_v/s
      bw_v <- w_v*as.vector(cswei(a_v,rs_cs,ind-1,0))
      
      if(eigen==TRUE)
      {
        if(inv==TRUE)
        {
          re <- invsph(sigma_i_s, der_t, si_d,w_v,X,rs_rs,rs_cs,ind-1,a_v[ind[,1]], bw_v, order, 1,tau,solver)
        }else{
          re <- invsph(sigma_s, der_t, s_d,w_v,X,rs_rs,rs_cs,ind-1,a_v[ind[,1]], bw_v, order, 0,tau,solver)
        }
        vi11 <- solve(re$xh%*%(X - re$wb_sig_i_hx_der[,(2:(n_c+1))]))
      }else{
        a_v_2 <- (as.vector(a_v)*as.vector(a_v))[ind[,1]]
        
        # w_v_x <- csqei(w_v,X,rs_rs-1,rs_cs-1,ind-1,a_v_2)
        xh <- t(bw_v*X - csqei(w_v,X,rs_rs,rs_cs,ind-1,a_v_2))
        wb_sig_i_hx_der <- t(xh)
        
        if(inv==TRUE)
        {
          diag(sigma_i_s_bk) <- si_d + bw_v*tau
          wb_sig_i_hx_der <- tau*Matrix::solve(sigma_i_s_bk,wb_sig_i_hx_der)
        }else{
          bw_i <- 1/bw_v
          diag(sigma_i_s_bk) <- s_d + bw_i/tau
          wb_sig_i_hx_der <- bw_i*Matrix::solve(sigma_i_s_bk,sigma_s%*%wb_sig_i_hx_der)
        }
        if(order>0)
        {
          wb_sig_i_hx_der_2 <- wb_sig_i_hx_der
          if(inv==TRUE)
          {
            for(i in 1:order)
            {
              wb_sig_i_hx_der_2 <- csqei(w_v,as.matrix(wb_sig_i_hx_der_2),rs_rs,rs_cs,ind-1,a_v_2)
              wb_sig_i_hx_der_2 <- tau*Matrix::solve(sigma_i_s_bk,wb_sig_i_hx_der_2)
              wb_sig_i_hx_der <- wb_sig_i_hx_der+wb_sig_i_hx_der_2
            }
          }else{
            for(i in 1:order)
            {
              wb_sig_i_hx_der_2 <- csqei(w_v,as.matrix(wb_sig_i_hx_der_2),rs_rs,rs_cs,ind-1,a_v_2)
              wb_sig_i_hx_der_2 <- bw_i*Matrix::solve(sigma_i_s_bk,sigma_s%*%wb_sig_i_hx_der_2)
              wb_sig_i_hx_der <- wb_sig_i_hx_der+wb_sig_i_hx_der_2
            }
          }
        }
        
        vi11 <- solve(xh%*%(X - wb_sig_i_hx_der))
      }
      
    }
  }else{
    a_v <- d_v/s
    # b_v <- as.vector(cswei(a_v,rs_cs-1,ind-1,0))
    bw_v <- w_v*as.vector(cswei(a_v,rs_cs,ind-1,0))
    
    a_v_p <- a_v[ind[,1]]
    a_v_p <- a_v_p[a_v_p>0]
    
    if(detap=='diagonal')
    {
      if(eigen==TRUE)
      {
        if(inv==TRUE)
        {
          logdet <- logdeth(sigma_i_s,si_d,bw_v, w_v,rs_cs_p-1,ind-1,a_v_p,tau,1,1)
          logdet <- logdet - n*log(tau)
        }else{
          logdet <- logdeth(sigma_s,s_d,bw_v, w_v,rs_cs_p-1,ind-1,a_v_p,tau,0,1)
          # logdet <- logdet + si_der
        }
      }else{
        wv2 <- w_v*w_v
        avp2 <- a_v_p*a_v_p
        qd <- sapply(1:n, function(x) wv2[x]*sum(avp2[1:rs_cs_p[ind[x,2]]]))
        if(inv==TRUE)
        {
          diag(sigma_i_s_bk) <- si_d + tau*(bw_v - qd)
          logdet <- abs(Matrix::determinant(Matrix::Cholesky(sigma_i_s_bk))$modulus*2)
          logdet <- logdet - n*log(tau)
        }else{
          dd <- bw_v - qd
          diag(sigma_i_s_bk) <- s_d + 1/tau/dd
          logdet <- abs(Matrix::determinant(Matrix::Cholesky(sigma_i_s_bk))$modulus*2) + sum(log(dd))
        }
      }
    }else{
      if(detap=='exact')
      {
        logdet <- logdeth(as(sigma_i_s,'dgCMatrix'),si_d,bw_v, w_v,rs_cs_p-1,ind-1,a_v_p,tau,1,0)
      }else{
        wv2 <- w_v*w_v
        avp2 <- a_v_p*a_v_p
        qd <- sapply(1:n, function(x) wv2[x]*sum(avp2[1:rs_cs_p[ind[x,2]]]))
        
        if(inv==TRUE)
        {
          v = sigma_i_s/tau
          diag(v) = diag(v) + bw_v - qd
          v = v*tau
          logdet = logdet_lanczos_sp(as(v,'dgCMatrix'), rad, 8) - n*log(tau)
        }else{
          dd <- bw_v - qd
          v = sigma_s%*%Diagonal(n,dd)
          diag(v) = diag(v) + 1/tau
          v = v%*%t(v)
          logdet = logdet_lanczos_sp(as(v,'dgCMatrix'), rad, 8)/2
        }
      }
        
    }
    
  }
  
  return(list(beta=beta_new,u=u_new,v11=vi11,iter=iter,ll=newloglik,logdet=logdet))
}
