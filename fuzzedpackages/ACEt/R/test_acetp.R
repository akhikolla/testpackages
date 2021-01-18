test_acetp <-
function(acetp, comp, sim = 100, robust = 0, pe = TRUE, verbose = TRUE)
{
	if(!(class(acetp) %in% c('AtCtEtp_model')))
	{
		stop('The first parameter must be an acetp object.')
	}

	if(!(comp %in% c('a','c','e')))
	{
		stop('The variable \'comp\' must be \'a\',\'c\' or \'e\' to specify which component to test linearity.')
	}

	if(acetp$mod[match(comp,c('a','c','e'))]=='c')
	{
		stop('The component to test is a constant.')
	}
  
  if(acetp$mod[match(comp,c('a','c','e'))]=='d')
  {
    if(verbose == TRUE)
	{
		cat("Model comparison: \n")
		print("Log-linear (null) vs. Spline")
	}
    
	mod_a <- acetp$mod
	mod_n <- mod_a
	mod_n[match(comp,c('a','c','e'))] <- 'l' 
	
	k_a <- ifelse((length(acetp$beta_a)-1)<2,8,length(acetp$beta_a)-1)
	k_c <- ifelse((length(acetp$beta_c)-1)<2,8,length(acetp$beta_c)-1)
	k_e <- ifelse((length(acetp$beta_e)-1)<2,8,length(acetp$beta_e)-1)
	data_m <- cbind(acetp$pheno_m[seq(from=1,to=nrow(acetp$pheno_m),by=2)],acetp$pheno_m[seq(from=2,to=nrow(acetp$pheno_m),by=2)])
	data_m <- cbind(data_m,acetp$T_m[seq(from=1,to=length(acetp$T_m),by=2)])
	data_d <- cbind(acetp$pheno_d[seq(from=1,to=nrow(acetp$pheno_d),by=2)],acetp$pheno_d[seq(from=2,to=nrow(acetp$pheno_d),by=2)])
	data_d <- cbind(data_d,acetp$T_d[seq(from=1,to=length(acetp$T_d),by=2)])	
	

	m_a <- AtCtEtp_2(data_m, data_d, knot_a=k_a, knot_c=k_c, knot_e=k_e, mod=mod_a, robust=robust)
	
	if(pe==FALSE)
	{
      p <- test_acetp_2(m_a, comp)
      return(p)
	}

	m_n <- AtCtEtp_2(data_m, data_d, knot_a=k_a, knot_c=k_c, knot_e=k_e, mod=mod_n, robust=robust)

	llr <- m_n$lik - m_a$lik

	
	order <- 3

	num_m <- nrow(data_m)
	num_d <- nrow(data_d)
	D_a_n <- m_n$D_a 
	D_c_n <- m_n$D_c
	D_e_n <- m_n$D_e
	
	if(m_n$mod[1]!='c')
	{
	  delta_a <- matrix(0, k_a+3-2-2, k_a+3-2)
	  for(i in 1:nrow(delta_a))
	  {
		  delta_a[i, i:(i+2)] <- c(1,-2,1)
	  }
	  D_a_n <- t(delta_a)%*%delta_a
	  ei_a <- eigen(D_a_n)
	  bb_a_m <- splineDesign(m_n$knot_a, x = data_m[,3], ord=order, outer.ok = TRUE)
	  bb_a_d <- splineDesign(m_n$knot_a, x = data_d[,3], ord=order, outer.ok = TRUE)
	  bb_a_m <- bb_a_m%*%ei_a$vectors
	  bb_a_d <- bb_a_d%*%ei_a$vectors
	}
	if(m_n$mod[2]!='c')
	{
	  delta_c <- matrix(0, k_c+3-2-2, k_c+3-2)
	  for(i in 1:nrow(delta_c))
	  {
		  delta_c[i, i:(i+2)] <- c(1,-2,1)
	  }
	  D_c_n <- t(delta_c)%*%delta_c
	  ei_c <- eigen(D_c_n)
	  bb_c_m <- splineDesign(m_n$knot_c, x = data_m[,3], ord=order, outer.ok = TRUE)
	  bb_c_d <- splineDesign(m_n$knot_c, x = data_d[,3], ord=order, outer.ok = TRUE)
	  bb_c_m <- bb_c_m%*%ei_c$vectors
	  bb_c_d <- bb_c_d%*%ei_c$vectors
	}
	if(m_n$mod[3]!='c')
	{
	  delta_e <- matrix(0, k_e+3-2-2, k_e+3-2)
	  for(i in 1:nrow(delta_e))
	  {
		  delta_e[i, i:(i+2)] <- c(1,-2,1)
	  }
	  D_e_n <- t(delta_e)%*%delta_e
	  ei_e <- eigen(D_e_n)
	  bb_e_m <- splineDesign(m_n$knot_e, x = data_m[,3], ord=order, outer.ok = TRUE)
	  bb_e_d <- splineDesign(m_n$knot_e, x = data_d[,3], ord=order, outer.ok = TRUE)
	  bb_e_m <- bb_e_m%*%ei_e$vectors
	  bb_e_d <- bb_e_d%*%ei_e$vectors
	}
	
	sim_m <- data_m
	sim_d <- data_d
 
	llr_sim <- rep(NA, sim)
 
	for(i in 1:sim)
	{
	  beta_a_n <- m_n$beta_a
	  beta_c_n <- m_n$beta_c
	  beta_e_n <- m_n$beta_e
	  var_b_a_n <- m_n$var_b_a
	  var_b_c_n <- m_n$var_b_c
	  var_b_e_n <- m_n$var_b_e
	  beta_a_1 <- m_a$beta_a
	  beta_c_1 <- m_a$beta_c
	  beta_e_1 <- m_a$beta_e
	  var_b_a_1 <- m_a$var_b_a
	  var_b_c_1 <- m_a$var_b_c
	  var_b_e_1 <- m_a$var_b_e
	  
	  
	if(m_n$mod[1]=='d')
	{
	  
	  beta_a_n[1:(ncol(bb_a_m)-2)] <- rnorm(ncol(bb_a_m)-2, mean=0, sd=sqrt(var_b_a_n/ei_a$values[1:(ncol(bb_a_m)-2)]))
	  a_m <- exp(bb_a_m%*%beta_a_n)
		a_d <- exp(bb_a_d%*%beta_a_n)
	}else{
	  	if(m_n$mod[1]=='l')
		  {			
			  a_m <- exp(bb_a_m[,(ncol(bb_a_m)-1):ncol(bb_a_m)]%*%beta_a_n)
			  a_d <- exp(bb_a_d[,(ncol(bb_a_d)-1):ncol(bb_a_d)]%*%beta_a_n)			
		  }else{			
			  a_m <- exp(beta_a_n*matrix(1, num_m, 1))
			  a_d <- exp(beta_a_n*matrix(1, num_d, 1))			
		  }
	}
	
	if(m_n$mod[2]=='d')
	{
	    beta_c_n[1:(ncol(bb_c_m)-2)] <- rnorm(ncol(bb_c_m)-2, mean=0, sd=sqrt(var_b_c_n/ei_c$values[1:(ncol(bb_c_m)-2)]))
	    c_m <- exp(bb_c_m%*%beta_c_n)
		  c_d <- exp(bb_c_d%*%beta_c_n)
		
	  }else{
		
	  	if(m_n$mod[2]=='l')
		  {	
			  c_m <- exp(bb_c_m[,(ncol(bb_c_m)-1):ncol(bb_c_m)]%*%beta_c_n)
			  c_d <- exp(bb_c_d[,(ncol(bb_c_d)-1):ncol(bb_c_d)]%*%beta_c_n)			
		  }else{			
			  c_m <- exp(beta_c_n*matrix(1, num_m, 1))
			  c_d <- exp(beta_c_n*matrix(1, num_d, 1))			
		  }
	}
	
	if(m_n$mod[3]=='d')
	{
	  beta_e_n[1:(ncol(bb_e_m)-2)] <- rnorm(ncol(bb_e_m)-2, mean=0, sd=sqrt(var_b_e_n/ei_e$values[1:(ncol(bb_e_m)-2)]))
	  
	  e_m <- exp(bb_e_m%*%beta_e_n)
		e_d <- exp(bb_e_d%*%beta_e_n)
		
	}else{
		
	  	if(m_n$mod[3]=='l')
		{
			
			e_m <- exp(bb_e_m[,(ncol(bb_e_m)-1):ncol(bb_e_m)]%*%beta_e_n)
			e_d <- exp(bb_e_d[,(ncol(bb_e_d)-1):ncol(bb_e_d)]%*%beta_e_n)
			
		}else{
			
			e_m <- exp(beta_e_n*matrix(1, num_m, 1))
			e_d <- exp(beta_e_n*matrix(1, num_d, 1))
		}
	}
		for(j in 1:nrow(sim_m))
		{
			sigma <- matrix(c(a_m[j]+c_m[j]+e_m[j],a_m[j]+c_m[j],a_m[j]+c_m[j],a_m[j]+c_m[j]+e_m[j]),2,2)
			sim_m[j,1:2] <- mvrnorm(1, rep(0,2), sigma)
		}
		for(j in 1:nrow(sim_d))
		{
			sigma <- matrix(c(a_d[j]+c_d[j]+e_d[j],0.5*a_d[j]+c_d[j],0.5*a_d[j]+c_d[j],a_d[j]+c_d[j]+e_d[j]),2,2)
			sim_d[j,1:2] <- mvrnorm(1, rep(0,2), sigma)
		}

		s_a <- AtCtEtp_2(sim_m, sim_d, knot_a=k_a, knot_c=k_c, knot_e=k_e, mod=mod_a, robust=robust)

	  s_n <- AtCtEtp_2(sim_m, sim_d, knot_a=k_a, knot_c=k_c, knot_e=k_e, mod=mod_n, robust=robust)

		llr_sim[i] <- s_n$lik-s_a$lik
	}

	p <- sum(llr_sim>llr)/sim
	test <- list(p = p, llr = llr, llr_sim=llr_sim)
  }else{
  
	if(verbose == TRUE)
	{
		cat("Model comparison: \n")
		print("Constancy (null) vs. Log-linear")
	}
    
    re <- acetp_mcmc(acetp,iter_num=10000)
    num_v <- 0
    if(sum(acetp$mod=='d')>0)
    {
      ind <- c()
      if(acetp$mod[1]=='d')
      {ind <- c(ind, 1)}
      if(acetp$mod[2]=='d')
      {ind <- c(ind, 2)}
      if(acetp$mod[3]=='d')
      {ind <- c(ind, 3)}
      hessian <- solve(acetp$hessian[ind,ind])
      num_v <- nrow(hessian)
      for(k in 1:num_v)
      {
        if(hessian[k,k]<0)
        {
          hessian[k,] <- 0
          hessian[,k] <- 0
        }
      }
      var <- c()
      if(acetp$mod[1]=='d')
       {var <- c(var, acetp$var_b_a)}
      if(acetp$mod[2]=='d')
        {var <- c(var, acetp$var_b_c)}
      if(acetp$mod[3]=='d')
        {var <- c(var, acetp$var_b_e)}
    }
    index <- c(1,2)
    if(comp=='a')
    {beta_t <- re$beta_a_mc}
    if(comp=='c')
    {
      beta_t <- re$beta_c_mc
      index <- index+length(re$beta_a_mc)
    }
    if(comp=='e')
    {
      beta_t <- re$beta_e_mc
      index <- index+length(re$beta_a_mc)+length(re$beta_c_mc)
    }
    cov_cm <- matrix(0,2,2)
    if(sum(acetp$mod=='d')>0)
    {
    cov_cor <- matrix(0,20,2) 
    for(j in 1:20)
    {
      var_t <- mvrnorm(1, var, hessian)
      if(acetp$mod[1]=='d')
      {
        acetp$var_b_a <- var_t[1]
        var_t <- var_t[-1]
      }
      if(acetp$mod[2]=='d')
      {acetp$var_b_c <- var_t[1]
       var_t <- var_t[-1]}
      if(acetp$mod[3]=='d')
      {acetp$var_b_e <- var_t[1]
       var_t <- var_t[-1]}
      re_b <- acetp_mcmc(acetp, iter_num = 4000, burnin = 500)
      if(comp=='a')
      {cov_cor[j,] <- re_b$beta_a_mc}
      if(comp=='c')
      {cov_cor[j,] <- re_b$beta_c_mc}
      if(comp=='e')
      {cov_cor[j,] <- re_b$beta_e_mc}
    }
    cov_cm <- cov(cov_cor)
    }
    sigma<-c(1,-1)%*%(re$cov_mc[index,index]+cov_cm)%*%c(1,-1)
    
    p <- pchisq((beta_t[1]-beta_t[2])^2/sigma,1,lower.tail=FALSE)
    test <- list(p = p, chisq = (beta_t[1]-beta_t[2])^2/sigma)
  }
  
	if(verbose == TRUE)
	{
		return(test)
	}else{
		return(invisible(test))
	}
	
}