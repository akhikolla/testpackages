test_acetp_2 <-
  function(acetp, comp)
  {
    
      # re <- acetp_mcmc_2(acetp,iter_num=iter_num, sd = 0.05, burnin=1000)
      order <- 3

    penal_a <- 2
    penal_c <- 2
    penal_e <- 2
    if(acetp$mod[1]=='c')
    {penal_a <- 1}
    if(acetp$mod[2]=='c')
    {penal_c <- 1}
    if(acetp$mod[3]=='c')
    {penal_e <- 1}
    
    num_m <- length(acetp$T_m)
    num_d <- length(acetp$T_d)
    
    delta_a <- matrix(0, length(acetp$knot_a)-4+order-2-penal_a, length(acetp$knot_a)-4+order-2)
    for(i in 1:nrow(delta_a))
    {
	    if(penal_a==2)
	    {delta_a[i, i:(i+2)] <- c(1,-2,1)}else{
		    delta_a[i, i:(i+1)] <- c(1,-1)
	    }
    }
    D_a <- t(delta_a)%*%delta_a

    B_des_a_m <- splineDesign(acetp$knot_a, x=acetp$T_m, ord=order)
    B_des_a_d <- splineDesign(acetp$knot_a, x=acetp$T_d, ord=order)
    ei_a <- eigen(D_a)
    B_des_a_m <- B_des_a_m%*%ei_a$vectors
    B_des_a_d <- B_des_a_d%*%ei_a$vectors
    D_a <- diag(c(ei_a$values[1:(length(ei_a$values)-2)],0,0))
    if(acetp$mod[1]=='l')
    {
      D_a <- matrix(0,2,2)
      B_des_a_m <- B_des_a_m[,(ncol(B_des_a_m)-1):ncol(B_des_a_m)]
      B_des_a_d <- B_des_a_d[,(ncol(B_des_a_d)-1):ncol(B_des_a_d)]

    }
    if(acetp$mod[1]=='c')
    {
        D_a <- matrix(0,1,1)
        B_des_a_m <- matrix(1, num_m, 1)
        B_des_a_d <- matrix(1, num_d, 1)
    }

    
    delta_c <- matrix(0, length(acetp$knot_c)-4+order-2-penal_c, length(acetp$knot_c)-4+order-2)
    for(i in 1:nrow(delta_c))
    {
	    if(penal_c==2)
	    {delta_c[i, i:(i+2)] <- c(1,-2,1)}else{
		    delta_c[i, i:(i+1)] <- c(1,-1)
	    }
    }
    D_c <- t(delta_c)%*%delta_c

    B_des_c_m <- splineDesign(acetp$knot_c, x=acetp$T_m, ord=order)
    B_des_c_d <- splineDesign(acetp$knot_c, x=acetp$T_d, ord=order)
    ei_c <- eigen(D_c)
    B_des_c_m <- B_des_c_m%*%ei_c$vectors
    B_des_c_d <- B_des_c_d%*%ei_c$vectors
    D_c <- diag(c(ei_c$values[1:(length(ei_c$values)-2)],0,0))
    if(acetp$mod[2]=='l')
    {
      D_c <- matrix(0,2,2)
      B_des_c_m <- B_des_c_m[,(ncol(B_des_c_m)-1):ncol(B_des_c_m)]
      B_des_c_d <- B_des_c_d[,(ncol(B_des_c_d)-1):ncol(B_des_c_d)]

    }
    if(acetp$mod[2]=='c')
    {
        D_c <- matrix(0,1,1)
        B_des_c_m <- matrix(1, num_m, 1)
        B_des_c_d <- matrix(1, num_d, 1)
    }
    
    
    delta_e <- matrix(0, length(acetp$knot_e)-4+order-2-penal_e, length(acetp$knot_e)-4+order-2)
    for(i in 1:nrow(delta_e))
    {
	    if(penal_e==2)
	    {delta_e[i, i:(i+2)] <- c(1,-2,1)}else{
		    delta_e[i, i:(i+1)] <- c(1,-1)
	    }
    }
    D_e <- t(delta_e)%*%delta_e


    B_des_e_m <- splineDesign(acetp$knot_e, x=acetp$T_m, ord=order)
    B_des_e_d <- splineDesign(acetp$knot_e, x=acetp$T_d, ord=order)
    ei_e <- eigen(D_e)
    B_des_e_m <- B_des_e_m%*%ei_e$vectors
    B_des_e_d <- B_des_e_d%*%ei_e$vectors
    D_e <- diag(c(ei_e$values[1:(length(ei_e$values)-2)],0,0))
    if(acetp$mod[3]=='l')
    {
      D_e <- matrix(0,2,2)
      B_des_e_m <- B_des_e_m[,(ncol(B_des_e_m)-1):ncol(B_des_e_m)]
      B_des_e_d <- B_des_e_d[,(ncol(B_des_e_d)-1):ncol(B_des_e_d)]

    }
    if(acetp$mod[3]=='c')
    {
        D_e <- matrix(0,1,1)
        B_des_e_m <- matrix(1, num_m, 1)
        B_des_e_d <- matrix(1, num_d, 1)
    }

    
      n_a <- length(acetp$beta_a)
      n_c <- length(acetp$beta_c)
      n_e <- length(acetp$beta_e)
      
      low_a <- -12
      upp_a <- 12
      low_c <- -12
      upp_c <- 12
      low_e <- rep(-8,n_e)
      #low_e[c(n_e-1,n_e)] <- rep(-8,2)
      #upp_e <- 10
      upp_e <- rep(8,n_e)
      #upp_e[c(n_e-1,n_e)] <- rep(8,2)
      
      if(acetp$var_b_a>0)
      {acetp$var_b_a <- 100000000}
      if(acetp$var_b_c>0)
      {acetp$var_b_c <- 100000000}
      if(acetp$var_b_e>0)
      {acetp$var_b_e <- 100000000}
    
      
      result <- optim(runif(n_a+n_c+n_e,min=-0.5,max=0.5), loglik_AtCtEt_epsp_g, gr_AtCtEt_epsp_g, pheno_m = acetp$pheno_m, pheno_d = acetp$pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_c_m=B_des_c_m, B_des_c_d=B_des_c_d, B_des_e_m=B_des_e_m, B_des_e_d=B_des_e_d, var_b_a=acetp$var_b_a, var_b_c=acetp$var_b_c, var_b_e=acetp$var_b_e, D_a=D_a, D_c=D_c, D_e=D_e, lower = c(rep(low_a,n_a),rep(low_c,n_c),low_e), upper = c(rep(upp_a,n_a),rep(upp_c,n_c),upp_e), method = "L-BFGS-B", control=list(maxit = 3000), hessian = TRUE)
    
      if(comp=='a')
      {
        index <- 1:(n_a-2)
      }
      if(comp=='c')
      {
        index <- (1:(n_c-2))+n_a
      }
      if(comp=='e')
      {
        index <- (1:(n_e-2))+n_c+n_a
        
      }
      beta_t <- result$par[index]
        
      sigma <-  solve(result$hessian)[index,index]
    
      chisq_t <- t(beta_t)%*%solve(sigma)%*%beta_t
      p_1 <- pchisq(chisq_t,length(index),lower.tail=FALSE)
      test <- list(p = p_1, chisq = chisq_t)
      return(test)
    
  }