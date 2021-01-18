acetp_mcmc <- function(acetp, iter_num = 10000, sd = 0.1, burnin =1000)
{
	if(!(class(acetp) %in% c('AtCtEp_model','AtEtp_model','AtCtEtp_model')))
	{
		stop('The first parameter must be an acetp object.')
	}

	
	if(burnin >= iter_num)
	{
		stop('The number of burnins must be smaller than the number of MCMC iterations.')
	}


	if(class(acetp)=='AtCtEtp_model')
	{
		res <- AtCtEtp_mcmc(acetp, iter_num, sd, burnin)
		return(res)
	}

	#if(class(acetp)=='AtEtp_model')
	#{
	#	res <- AtEtp_mcmc(acetp, iter_num, sd, burnin)
	#	return(res)
	#}
	
}




AtCtEtp_mcmc <-
function(AtCtEtp, iter_num = 10000, sd = 0.1, burnin =1000)
{

if(class(AtCtEtp)!='AtCtEtp_model')
{
	stop('The first parameter must be an object obtained from the AtCtEtp function.')
}

T_m <- AtCtEtp$T_m
num_m <- length(T_m)
T_d <- AtCtEtp$T_d
num_d <- length(T_d)

max_t <- max(c(T_m,T_d))
min_t <- min(c(T_m,T_d))

t_int <- max_t-min_t
l_m_1 <- (max_t-T_m)/t_int
l_m_2 <- (T_m-min_t)/t_int
l_d_1 <- (max_t-T_d)/t_int
l_d_2 <- (T_d-min_t)/t_int

order <- 3
if(length(AtCtEtp$beta_a)>2)
{
	B_des_a_m <- splineDesign(AtCtEtp$knot_a, x=AtCtEtp$T_m, ord=order)
	B_des_a_d <- splineDesign(AtCtEtp$knot_a, x=AtCtEtp$T_d, ord=order)
}else{
	if(length(AtCtEtp$beta_a)==2)
	{
		B_des_a_m <- matrix(NA, num_m, 2)
		B_des_a_m[,1] <- l_m_1
		B_des_a_m[,2] <- l_m_2
		B_des_a_d <- matrix(NA, num_d, 2)
		B_des_a_d[,1] <- l_d_1
		B_des_a_d[,2] <- l_d_2
	}else{
		B_des_a_m <- matrix(1, num_m, 1)
		B_des_a_d <- matrix(1, num_d, 1)
	}
}
if(length(AtCtEtp$beta_c)>2)
{
	B_des_c_m <- splineDesign(AtCtEtp$knot_c, x=AtCtEtp$T_m, ord=order)
	B_des_c_d <- splineDesign(AtCtEtp$knot_c, x=AtCtEtp$T_d, ord=order)
}else{
	if(length(AtCtEtp$beta_c)==2)
	{
		B_des_c_m <- matrix(NA, num_m, 2)
		B_des_c_m[,1] <- l_m_1
		B_des_c_m[,2] <- l_m_2
		B_des_c_d <- matrix(NA, num_d, 2)
		B_des_c_d[,1] <- l_d_1
		B_des_c_d[,2] <- l_d_2
	}else{
		B_des_c_m <- matrix(1, num_m, 1)
		B_des_c_d <- matrix(1, num_d, 1)
	}
}
if(length(AtCtEtp$beta_e)>2)
{
	B_des_e_m <- splineDesign(AtCtEtp$knot_e, x=AtCtEtp$T_m, ord=order)
	B_des_e_d <- splineDesign(AtCtEtp$knot_e, x=AtCtEtp$T_d, ord=order)
}else{
	if(length(AtCtEtp$beta_e)==2)
	{
		B_des_e_m <- matrix(NA, num_m, 2)
		B_des_e_m[,1] <- l_m_1
		B_des_e_m[,2] <- l_m_2
		B_des_e_d <- matrix(NA, num_d, 2)
		B_des_e_d[,1] <- l_d_1
		B_des_e_d[,2] <- l_d_2
	}else{
		B_des_e_m <- matrix(1, num_m, 1)
		B_des_e_d <- matrix(1, num_d, 1)
	}
}

result <- mcmc_epsp_AtCtEt(AtCtEtp$pheno_m, AtCtEtp$pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d, AtCtEtp$var_b_a, AtCtEtp$var_b_c, AtCtEtp$var_b_e, AtCtEtp$beta_a, AtCtEtp$beta_c, AtCtEtp$beta_e, AtCtEtp$D_a, AtCtEtp$D_c, AtCtEtp$D_e, iter_num, burnin, sd)

AtCtEtp_mc_mod <- list(beta_a_mc=result$beta_a_mc, beta_c_mc=result$beta_c_mc, beta_e_mc=result$beta_e_mc, cov_mc = result$cov, knots_a=AtCtEtp$knot_a, knots_c=AtCtEtp$knot_c, knots_e=AtCtEtp$knot_e, min_t = min(AtCtEtp$T_m, AtCtEtp$T_d), max_t = max(AtCtEtp$T_m, AtCtEtp$T_d))

class(AtCtEtp_mc_mod) <- 'AtCtEtp_mc_model'

return(AtCtEtp_mc_mod)
}


acetp_mcmc_2 <- function(AtCtEtp, iter_num = 5000, sd = 0.1, burnin =500)
  {
    
    if(class(AtCtEtp)!='AtCtEtp_model')
    {
      stop('The first parameter must be an object obtained from the AtCtEtp function.')
    }
    
    T_m <- AtCtEtp$T_m
    num_m <- length(T_m)
    T_d <- AtCtEtp$T_d
    num_d <- length(T_d)
    
    t_int <- max(c(T_m,T_d))-min(c(T_m,T_d))
    
    order <- 3
    if(length(AtCtEtp$beta_a)>2)
    {
      ei_a <- eigen(AtCtEtp$D_a)
      B_des_a_m <- splineDesign(AtCtEtp$knot_a, x=AtCtEtp$T_m, ord=order)
      B_des_a_d <- splineDesign(AtCtEtp$knot_a, x=AtCtEtp$T_d, ord=order)
      B_des_a_m <- B_des_a_m%*%ei_a$vectors
      B_des_a_d <- B_des_a_d%*%ei_a$vectors
      D_a <- diag(c(ei_a$values[1:(length(ei_a$values)-2)],0,0))
    }else{
      if(length(AtCtEtp$beta_a)==2)
      {
        k_a <- length(AtCtEtp$knot_a)-4
        delta_a <- matrix(0, k_a+3-2-2, k_a+3-2)
        for(i in 1:nrow(delta_a))
        {
          delta_a[i, i:(i+2)] <- c(1,-2,1)
        }
        D_a_n <- t(delta_a)%*%delta_a
        ei_a <- eigen(D_a_n)
        B_des_a_m <- splineDesign(AtCtEtp$knot_a, x=AtCtEtp$T_m, ord=order)
        B_des_a_d <- splineDesign(AtCtEtp$knot_a, x=AtCtEtp$T_d, ord=order)
        B_des_a_m <- B_des_a_m%*%ei_a$vectors
        B_des_a_d <- B_des_a_d%*%ei_a$vectors
        B_des_a_m <- B_des_a_m[,(ncol(B_des_a_m)-1):ncol(B_des_a_m)]
        B_des_a_d <- B_des_a_d[,(ncol(B_des_a_d)-1):ncol(B_des_a_d)]
        D_a <- AtCtEtp$D_a
      }else{
        B_des_a_m <- matrix(1, num_m, 1)
        B_des_a_d <- matrix(1, num_d, 1)
        D_a <- AtCtEtp$D_a
      }
    }
    
    if(length(AtCtEtp$beta_c)>2)
    {
      ei_c <- eigen(AtCtEtp$D_c)
      B_des_c_m <- splineDesign(AtCtEtp$knot_c, x=AtCtEtp$T_m, ord=order)
      B_des_c_d <- splineDesign(AtCtEtp$knot_c, x=AtCtEtp$T_d, ord=order)
      B_des_c_m <- B_des_c_m%*%ei_c$vectors
      B_des_c_d <- B_des_c_d%*%ei_c$vectors
      D_c <- diag(c(ei_c$values[1:(length(ei_c$values)-2)],0,0))
    }else{
      if(length(AtCtEtp$beta_c)==2)
      {
        k_c <- length(AtCtEtp$knot_c)-4
        delta_c <- matrix(0, k_c+3-2-2, k_c+3-2)
        for(i in 1:nrow(delta_c))
        {
          delta_c[i, i:(i+2)] <- c(1,-2,1)
        }
        D_c_n <- t(delta_c)%*%delta_c
        ei_c <- eigen(D_c_n)
        B_des_c_m <- splineDesign(AtCtEtp$knot_c, x=AtCtEtp$T_m, ord=order)
        B_des_c_d <- splineDesign(AtCtEtp$knot_c, x=AtCtEtp$T_d, ord=order)
        B_des_c_m <- B_des_c_m%*%ei_c$vectors
        B_des_c_d <- B_des_c_d%*%ei_c$vectors
        B_des_c_m <- B_des_c_m[,(ncol(B_des_c_m)-1):ncol(B_des_c_m)]
        B_des_c_d <- B_des_c_d[,(ncol(B_des_c_d)-1):ncol(B_des_c_d)]
        D_c <- AtCtEtp$D_c
      }else{
        B_des_c_m <- matrix(1, num_m, 1)
        B_des_c_d <- matrix(1, num_d, 1)
        D_c <- AtCtEtp$D_c
      }
    }
    if(length(AtCtEtp$beta_e)>2)
    {
      ei_e <- eigen(AtCtEtp$D_e)
      B_des_e_m <- splineDesign(AtCtEtp$knot_e, x=AtCtEtp$T_m, ord=order)
      B_des_e_d <- splineDesign(AtCtEtp$knot_e, x=AtCtEtp$T_d, ord=order)
      B_des_e_m <- B_des_e_m%*%ei_e$vectors
      B_des_e_d <- B_des_e_d%*%ei_e$vectors
      D_e <- diag(c(ei_e$values[1:(length(ei_e$values)-2)],0,0))
    }else{
      if(length(AtCtEtp$beta_e)==2)
      {
        k_e <- length(AtCtEtp$knot_e)-4
        delta_e <- matrix(0, k_e+3-2-2, k_e+3-2)
        for(i in 1:nrow(delta_e))
        {
          delta_e[i, i:(i+2)] <- c(1,-2,1)
        }
        D_e_n <- t(delta_e)%*%delta_e
        ei_e <- eigen(D_e_n)
        B_des_e_m <- splineDesign(AtCtEtp$knot_e, x=AtCtEtp$T_m, ord=order)
        B_des_e_d <- splineDesign(AtCtEtp$knot_e, x=AtCtEtp$T_d, ord=order)
        B_des_e_m <- B_des_e_m%*%ei_e$vectors
        B_des_e_d <- B_des_e_d%*%ei_e$vectors
        B_des_e_m <- B_des_e_m[,(ncol(B_des_e_m)-1):ncol(B_des_e_m)]
        B_des_e_d <- B_des_e_d[,(ncol(B_des_e_d)-1):ncol(B_des_e_d)]
        D_e <- AtCtEtp$D_e
      }else{
        B_des_e_m <- matrix(1, num_m, 1)
        B_des_e_d <- matrix(1, num_d, 1)
        D_e <- AtCtEtp$D_e
      }
    }
    
    result <- mcmc_epsp_AtCtEt_2(AtCtEtp$pheno_m, AtCtEtp$pheno_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d, AtCtEtp$var_b_a, AtCtEtp$var_b_c, AtCtEtp$var_b_e, D_a, D_c, D_e, iter_num, burnin, sd)
    
    AtCtEtp_mc_mod <- list(beta_a_mc=result$beta_a_mc, beta_c_mc=result$beta_c_mc, beta_e_mc=result$beta_e_mc, cov_mc = result$cov, knots_a=AtCtEtp$knot_a, knots_c=AtCtEtp$knot_c, knots_e=AtCtEtp$knot_e, min_t = min(AtCtEtp$T_m, AtCtEtp$T_d), max_t = max(AtCtEtp$T_m, AtCtEtp$T_d))
    
    class(AtCtEtp_mc_mod) <- 'AtCtEtp_mc_model'
    
    return(AtCtEtp_mc_mod)
  }

