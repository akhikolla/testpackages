AtCtEtp_2 <-
function(data_m, data_d, knot_a=8, knot_c=8, knot_e=8, eps = 0.1, mod=c('d','d','d'), robust=2)
{

if((is.vector(mod)==FALSE) | (length(mod)!=3) )
{stop('The model parameter must be a vector of length 3.')}

if(!(mod[1] %in% c('d','c','l')))
{stop('The \'mod\' parameter for the A component must be \'d\'(dynamic), \'c\'(constant) or \'l\'(linear).')}

if(!(mod[2] %in% c('d','c','l')))
{stop('The \'mod\' parameter for the C component must be \'d\'(dynamic), \'c\'(constant) or \'l\'(linear).')}

if(!(mod[3] %in% c('d','c','l')))
{stop('The \'mod\' parameter for the E component must be \'d\'(dynamic), \'c\'(constant) or \'l\'(linear).')}

if((knot_a<3)|(knot_c<3)|(knot_e<3))
{stop('The number of knots must be larger than 2.')}

num_m <- nrow(data_m)*2
num_d <- nrow(data_d)*2
pheno_m <- matrix(NA, num_m, 1)
pheno_d <- matrix(NA, num_d, 1)
pheno_m[seq(from=1, to=num_m, by=2),1] <- data_m[,1]
pheno_m[seq(from=2, to=num_m, by=2),1] <- data_m[,2]
pheno_d[seq(from=1, to=num_d, by=2),1] <- data_d[,1]
pheno_d[seq(from=2, to=num_d, by=2),1] <- data_d[,2]
T_m <- rep(data_m[,3],each=2)
T_d <- rep(data_d[,3],each=2)

mag <- var(pheno_m)
init_max <- log(mag)
init_min <- log(mag) - abs(log(mag))*1.2
limit <- 12
limit_e <- 10

low_var <- 1e-06
upp_var <- 100

eps <- eps*2

order <- 3

penal_a <- 2
penal_c <- 2
penal_e <- 2
if(mod[1]=='c')
{penal_a <- 1}
if(mod[2]=='c')
{penal_c <- 1}
if(mod[3]=='c')
{penal_e <- 1}

t_int <- max(c(T_m,T_d))-min(c(T_m,T_d))

delta_a <- matrix(0, knot_a+order-2-penal_a, knot_a+order-2)
for(i in 1:nrow(delta_a))
{
	if(penal_a==2)
	{delta_a[i, i:(i+2)] <- c(1,-2,1)}else{
		delta_a[i, i:(i+1)] <- c(1,-1)
	}
}
D_a <- t(delta_a)%*%delta_a
knots_a <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=knot_a)
interval_a <- knots_a[2] - knots_a[1]
knots_a <- c(c(min(T_m, T_d)-interval_a*2,min(T_m, T_d)-interval_a), knots_a)
knots_a <- c(knots_a, c(max(T_m, T_d)+interval_a,max(T_m, T_d)+interval_a*2))
B_des_a_m <- splineDesign(knots_a, x=T_m, ord=order)
B_des_a_d <- splineDesign(knots_a, x=T_d, ord=order)
ei_a <- eigen(D_a)
B_des_a_m <- B_des_a_m%*%ei_a$vectors
B_des_a_d <- B_des_a_d%*%ei_a$vectors
D_a <- diag(c(ei_a$values[1:(length(ei_a$values)-2)],0,0))
if(mod[1]=='l')
{
D_a <- matrix(0,2,2)
B_des_a_m <- B_des_a_m[,(ncol(B_des_a_m)-1):ncol(B_des_a_m)]
B_des_a_d <- B_des_a_d[,(ncol(B_des_a_d)-1):ncol(B_des_a_d)]
# knots_a <- c(min(T_m, T_d),max(T_m, T_d))
}
if(mod[1]=='c')
{
D_a <- matrix(0,1,1)
B_des_a_m <- matrix(1, num_m, 1)
B_des_a_d <- matrix(1, num_d, 1)
# knots_a <- c(min(T_m, T_d))
}

delta_c <- matrix(0, knot_c+order-2-penal_c, knot_c+order-2)
for(i in 1:nrow(delta_c))
{
	if(penal_c==2)
	{delta_c[i, i:(i+2)] <- c(1,-2,1)}else{
		delta_c[i, i:(i+1)] <- c(1,-1)
	}
}
D_c <- t(delta_c)%*%delta_c

knots_c <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=knot_c)
interval_c <- knots_c[2] - knots_c[1]
knots_c <- c(c(min(T_m, T_d)-interval_c*2,min(T_m, T_d)-interval_c), knots_c)
knots_c <- c(knots_c, c(max(T_m, T_d)+interval_c,max(T_m, T_d)+interval_c*2))
B_des_c_m <- splineDesign(knots_c, x=T_m, ord=order)
B_des_c_d <- splineDesign(knots_c, x=T_d, ord=order)
ei_c <- eigen(D_c)
B_des_c_m <- B_des_c_m%*%ei_c$vectors
B_des_c_d <- B_des_c_d%*%ei_c$vectors
D_c <- diag(c(ei_c$values[1:(length(ei_c$values)-2)],0,0))

if(mod[2]=='l')
{
	D_c <- matrix(0,2,2)
	B_des_c_m <- B_des_c_m[,(ncol(B_des_c_m)-1):ncol(B_des_c_m)]
  B_des_c_d <- B_des_c_d[,(ncol(B_des_c_d)-1):ncol(B_des_c_d)]
	# knots_c <- c(min(T_m, T_d),max(T_m, T_d))
}
if(mod[2]=='c')
{
	B_des_c_m <- matrix(1, num_m, 1)
	B_des_c_d <- matrix(1, num_d, 1)
	D_c <- matrix(0,1,1)
	# knots_c <- c(min(T_m, T_d))
}


delta_e <- matrix(0, knot_e+order-2-penal_e, knot_e+order-2)
for(i in 1:nrow(delta_e))
{
	if(penal_e==2)
	{delta_e[i, i:(i+2)] <- c(1,-2,1)}else{
		delta_e[i, i:(i+1)] <- c(1,-1)
	}
}
D_e <- t(delta_e)%*%delta_e

knots_e <- seq(from=min(T_m, T_d), to=max(T_m, T_d), length.out=knot_e)
interval_e <- knots_e[2] - knots_e[1]
knots_e <- c(c(min(T_m, T_d)-interval_e*2,min(T_m, T_d)-interval_e), knots_e)
knots_e <- c(knots_e, c(max(T_m, T_d)+interval_e,max(T_m, T_d)+interval_e*2))
B_des_e_m <- splineDesign(knots_e, x=T_m, ord=order)
B_des_e_d <- splineDesign(knots_e, x=T_d, ord=order)
ei_e <- eigen(D_e)
B_des_e_m <- B_des_e_m%*%ei_e$vectors
B_des_e_d <- B_des_e_d%*%ei_e$vectors
D_e <- diag(c(ei_e$values[1:(length(ei_e$values)-2)],0,0))

if(mod[3]=='l')
{
	D_e <- matrix(0,2,2)
	B_des_e_m <- B_des_e_m[,(ncol(B_des_e_m)-1):ncol(B_des_e_m)]
  B_des_e_d <- B_des_e_d[,(ncol(B_des_e_d)-1):ncol(B_des_e_d)]
	# knots_e <- c(min(T_m, T_d),max(T_m, T_d))
}
if(mod[3]=='c')
{
	B_des_e_m <- matrix(1, num_m, 1)
	B_des_e_d <- matrix(1, num_d, 1)
	D_e <- matrix(0,1,1)
	# knots_e <- c(min(T_m, T_d))
}

n_a <- ncol(B_des_a_m)
n_c <- ncol(B_des_c_m)
n_e <- ncol(B_des_e_m)

lower <- 0

var_b_a <- runif(1,min=0.5*abs(log(mag)),max=1.5*abs(log(mag)))
if(mod[1] %in% c('l','c'))
{
	var_b_a <- lower
	n_a <- ifelse(mod[1]=='l',2,1) 
}
var_b_c <- runif(1,min=0.5*abs(log(mag)),max=1.5*abs(log(mag)))
if(mod[2] %in% c('l','c'))
{
	var_b_c <- lower
	n_c <- ifelse(mod[2]=='l',2,1) 
}
var_b_e <- runif(1,min=0.5*abs(log(mag)),max=1.5*abs(log(mag)))
if(mod[3] %in% c('l','c'))
{
	var_b_e <- lower
	n_e <- ifelse(mod[3]=='l',2,1) 
}

beta_a <- runif(n_a,min=-1,max=1)
beta_c <- runif(n_c,min=-1,max=1)
beta_e <- runif(n_e,min=-1,max=1)

lik <- 100000
lik_pre <- 200000

liks <- c()
betas <- matrix(0,0,n_a+n_c+n_e)
vars <- matrix(0,0,3)
if((mod[1]!='d')&(mod[2]!='d')&(mod[3]!='d'))
{
  low_a <- -15
  upp_a <- 15
  low_c <- -15
  upp_c <- 15
  low_e <- -9
  upp_e <- 9
  result <- optim(c(beta_a,beta_c,beta_e), loglik_AtCtEt_epsp_g, gr_AtCtEt_epsp_g, pheno_m = pheno_m, pheno_d = pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_c_m=B_des_c_m, B_des_c_d=B_des_c_d, B_des_e_m=B_des_e_m, B_des_e_d=B_des_e_d, var_b_a=0, var_b_c=0, var_b_e=0, D_a=D_a, D_c=D_c, D_e=D_e, lower = c(rep(low_a,n_a),rep(low_c,n_c),rep(low_e,n_e)), upper = c(rep(upp_a,n_a),rep(upp_c,n_c),rep(upp_e,n_e)), method = "L-BFGS-B", control=list(maxit = 3000))
  if(robust>0)
  {
    for(j in 1:ceiling(robust))
    {
      init <- runif(n_a+n_c+n_e, min=init_min, max=init_max)
      result_r <- optim(init, loglik_AtCtEt_epsp_g, gr_AtCtEt_epsp_g, pheno_m = pheno_m, pheno_d = pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_c_m=B_des_c_m, B_des_c_d=B_des_c_d, B_des_e_m=B_des_e_m, B_des_e_d=B_des_e_d, var_b_a=0, var_b_c=0, var_b_e=0, D_a=D_a, D_c=D_c, D_e=D_e, lower = c(rep(low_a,n_a),rep(low_c,n_c),rep(low_e,n_e)), upper = c(rep(upp_a,n_a),rep(upp_c,n_c),rep(upp_e,n_e)), method = "L-BFGS-B", control=list(maxit = 3000))
      if(result_r$value < result$value)
      {
        result <- result_r
      }
    }
  }
  beta_a <- result$par[1:n_a]
  beta_c <- result$par[(1+n_a):(n_a+n_c)]
  beta_e <- result$par[(1+n_a+n_c):(n_a+n_c+n_e)]
  lik <- loglik_AtCtEt_epsp(c(0,0,0), pheno_m=pheno_m, pheno_d=pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, beta_a=beta_a, D_a=D_a, B_des_c_m=B_des_c_m, B_des_c_d=B_des_c_d, beta_c=beta_c, D_c=D_c, B_des_e_m=B_des_e_m, B_des_e_d=B_des_e_d, beta_e=beta_e, D_e=D_e)		
  AtCtEtp_model <- list(D_a = D_a, D_c = D_c, D_e = D_e, pheno_m = pheno_m, pheno_d = pheno_d, T_m = T_m, T_d = T_d, knot_a=knots_a, knot_c=knots_c, knot_e=knots_e, beta_a=beta_a, beta_c=beta_c, beta_e=beta_e, con=result$convergence, lik=lik/2, iter=lik/2, var_b_a=lower, var_b_c=lower, var_b_e=lower, mod=mod, bf = lik/2)
  
}else{
  while(abs(lik-lik_pre)>eps)
  {
    lik_pre <- lik
    if((mod[1]=='d')&(mod[2]=='d')&(mod[3]=='d'))
    {   
      result <- optim(c(beta_a,beta_c,beta_e), loglik_AtCtEt_epsp_g, gr_AtCtEt_epsp_g, pheno_m = pheno_m, pheno_d = pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_c_m=B_des_c_m, B_des_c_d=B_des_c_d, B_des_e_m=B_des_e_m, B_des_e_d=B_des_e_d, var_b_a=var_b_a, var_b_c=var_b_c, var_b_e=var_b_e, D_a=D_a, D_c=D_c, D_e=D_e, lower = c(rep((-1)*limit,n_a+n_c),rep((-1)*limit_e,n_e)), upper = c(rep(limit,n_a+n_c),rep(limit_e,n_e)), method = "L-BFGS-B", control=list(maxit = 3000))
      betas <- rbind(betas, result$par)
      beta_a <- result$par[1:n_a]
      beta_c <- result$par[(1+n_a):(n_a+n_c)]
      beta_e <- result$par[(1+n_a+n_c):(n_a+n_c+n_e)]
      
      result <- optim(c(var_b_a,var_b_c,var_b_e), loglik_AtCtEt_epsp, gr_AtCtEt_epsp, pheno_m = pheno_m, pheno_d = pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_c_m=B_des_c_m, B_des_c_d=B_des_c_d, B_des_e_m=B_des_e_m, B_des_e_d=B_des_e_d, beta_a=beta_a, beta_c=beta_c, beta_e = beta_e, D_a=D_a, D_c=D_c, D_e=D_e, lower = rep(low_var,3), upper = rep(upp_var,3), method = "L-BFGS-B", control=list(maxit = 3000))
    }else
    {	
      v_a_t <- var_b_a
      v_c_t <- var_b_c
      v_e_t <- var_b_e
      
      if(mod[1]!='d')
      {v_a_t <- 0}
      if(mod[2]!='d')
      {v_c_t <- 0}
      if(mod[3]!='d')
      {v_e_t <- 0}
      
      low_a <- (-1)*limit
      upp_a <- limit
      low_c <- (-1)*limit
      upp_c <- limit
      low_e <- rep((-1)*limit_e,n_e)
      upp_e <- rep(limit_e,n_e)
      
      result <- optim(c(beta_a,beta_c,beta_e), loglik_AtCtEt_epsp_g, gr_AtCtEt_epsp_g, pheno_m = pheno_m, pheno_d = pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_c_m=B_des_c_m, B_des_c_d=B_des_c_d, B_des_e_m=B_des_e_m, B_des_e_d=B_des_e_d, var_b_a=v_a_t, var_b_c=v_c_t, var_b_e=v_e_t, D_a=D_a, D_c=D_c, D_e=D_e, lower = c(rep(low_a,n_a),rep(low_c,n_c),low_e), upper = c(rep(upp_a,n_a),rep(upp_c,n_c),upp_e), method = "L-BFGS-B", control=list(maxit = 3000))
      beta_a <- result$par[1:n_a]
      beta_c <- result$par[(1+n_a):(n_a+n_c)]
      beta_e <- result$par[(1+n_a+n_c):(n_a+n_c+n_e)]
      
      betas <- rbind(betas,c(beta_a,beta_c,beta_e))
      
      low_a <- low_c <-  low_var
      low_e <- low_var
      upp_a <- upp_c <- upp_e <- upp_var
      if((mod[1]=='l')|(mod[1]=='c'))
      {low_a <- upp_a <- lower}
      if((mod[2]=='l')|(mod[2]=='c'))
      {low_c <- upp_c <- lower}
      if((mod[3]=='l')|(mod[3]=='c'))
      {low_e <- upp_e <- lower}
      
      result <- optim(c(var_b_a,var_b_c,var_b_e), loglik_AtCtEt_epsp, gr_AtCtEt_epsp, pheno_m = pheno_m, pheno_d = pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_c_m=B_des_c_m, B_des_c_d=B_des_c_d, B_des_e_m=B_des_e_m, B_des_e_d=B_des_e_d, beta_a=beta_a, beta_c=beta_c, beta_e=beta_e, D_a=D_a, D_c=D_c, D_e=D_e, lower = c(low_a, low_c, low_e), upper = c(upp_a, upp_c, upp_e), method = "L-BFGS-B", control=list(maxit = 3000))	
    }
    vars <- rbind(vars, result$par)
    var_b_a <- result$par[1]
    var_b_c <- result$par[2]
    var_b_e <- result$par[3]
    lik <- result$value
    liks <- c(liks, result$value)
  }
  
  min_i <- match(min(liks), liks)
  
  if(robust>0)
  {
    for(rob in 1:ceiling(robust))
    {
      lik <- 100000
      lik_pre <- 200000
      
      liks_r <- c()
      betas_r <- matrix(0,0,n_a+n_c+n_e)
      vars_r <- matrix(0,0,3)
      
      beta_a <- runif(n_a,min=init_min,max=init_max)
      beta_c <- runif(n_c,min=init_min,max=init_max)
      beta_e <- runif(n_e,min=init_min,max=init_max)
      if(var_b_a!=0)
      {var_b_a <- runif(1,min=0.5*abs(log(mag)),max=1.5*abs(log(mag)))}
      if(var_b_c!=0)
      {var_b_c <- runif(1,min=0.5*abs(log(mag)),max=1.5*abs(log(mag)))}
      if(var_b_e!=0)
      {var_b_e <- runif(1,min=0.5*abs(log(mag)),max=1.5*abs(log(mag)))}
      while(abs(lik-lik_pre)>eps)
      {
        lik_pre <- lik
        if((mod[1]=='d')&(mod[2]=='d')&(mod[3]=='d'))
        {
          result <- optim(c(beta_a,beta_c,beta_e), loglik_AtCtEt_epsp_g, gr_AtCtEt_epsp_g, pheno_m = pheno_m, pheno_d = pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_c_m=B_des_c_m, B_des_c_d=B_des_c_d, B_des_e_m=B_des_e_m, B_des_e_d=B_des_e_d, var_b_a=var_b_a, var_b_c=var_b_c, var_b_e=var_b_e, D_a=D_a, D_c=D_c, D_e=D_e, lower = c(rep(-12,n_a+n_c),rep(-9,n_e)), upper = c(rep(12,n_a+n_c),rep(9,n_e)), method = "L-BFGS-B", control=list(maxit = 3000))
          betas_r <- rbind(betas_r, result$par)
          beta_a <- result$par[1:n_a]
          beta_c <- result$par[(1+n_a):(n_a+n_c)]
          beta_e <- result$par[(1+n_a+n_c):(n_a+n_c+n_e)]
          result <- optim(c(var_b_a,var_b_c,var_b_e), loglik_AtCtEt_epsp, gr_AtCtEt_epsp, pheno_m = pheno_m, pheno_d = pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_c_m=B_des_c_m, B_des_c_d=B_des_c_d, B_des_e_m=B_des_e_m, B_des_e_d=B_des_e_d, beta_a=beta_a, beta_c=beta_c, beta_e = beta_e, D_a=D_a, D_c=D_c, D_e=D_e, lower = rep(low_var,3), upper = rep(upp_var,3), method = "L-BFGS-B", control=list(maxit = 3000))
        }else
        {	
          v_a_t <- var_b_a
          v_c_t <- var_b_c
          v_e_t <- var_b_e
          
          if(mod[1]!='d')
          {v_a_t <- 0}
          if(mod[2]!='d')
          {v_c_t <- 0}
          if(mod[3]!='d')
          {v_e_t <- 0}
          
          low_a <- (-1)*limit
          upp_a <- limit
          low_c <- (-1)*limit
          upp_c <- limit
          low_e <- rep((-1)*limit_e,n_e)
          upp_e <- rep(limit_e,n_e)
          
          result <- optim(c(beta_a,beta_c,beta_e), loglik_AtCtEt_epsp_g, gr_AtCtEt_epsp_g, pheno_m = pheno_m, pheno_d = pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_c_m=B_des_c_m, B_des_c_d=B_des_c_d, B_des_e_m=B_des_e_m, B_des_e_d=B_des_e_d, var_b_a=v_a_t, var_b_c=v_c_t, var_b_e=v_e_t, D_a=D_a, D_c=D_c, D_e=D_e, lower = c(rep(low_a,n_a),rep(low_c,n_c),low_e), upper = c(rep(upp_a,n_a),rep(upp_c,n_c),upp_e), method = "L-BFGS-B", control=list(maxit = 3000))
          beta_a <- result$par[1:n_a]
          beta_c <- result$par[(1+n_a):(n_a+n_c)]
          beta_e <- result$par[(1+n_a+n_c):(n_a+n_c+n_e)]
          
          betas_r <- rbind(betas_r,c(beta_a,beta_c,beta_e))
          
          low_a <- low_c <-  low_var
          low_e <- low_var
          upp_a <- upp_c <- upp_e <- upp_var
          if((mod[1]=='l')|(mod[1]=='c'))
          {low_a <- upp_a <- lower}
          if((mod[2]=='l')|(mod[2]=='c'))
          {low_c <- upp_c <- lower}
          if((mod[3]=='l')|(mod[3]=='c'))
          {low_e <- upp_e <- lower}
          
          result <- optim(c(var_b_a,var_b_c,var_b_e), loglik_AtCtEt_epsp, gr_AtCtEt_epsp, pheno_m = pheno_m, pheno_d = pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_c_m=B_des_c_m, B_des_c_d=B_des_c_d, B_des_e_m=B_des_e_m, B_des_e_d=B_des_e_d, beta_a=beta_a, beta_c=beta_c, beta_e=beta_e, D_a=D_a, D_c=D_c, D_e=D_e, lower = c(low_a, low_c, low_e), upper = c(upp_a, upp_c, upp_e), method = "L-BFGS-B", control=list(maxit = 3000))	
        }
        vars_r <- rbind(vars_r, result$par)
        var_b_a <- result$par[1]
        var_b_c <- result$par[2]
        var_b_e <- result$par[3]
        lik <- result$value
        liks_r <- c(liks_r, result$value)
      }
      if(min(liks_r)<min(liks))
      {
        liks <- liks_r
        vars <- vars_r
        betas <- betas_r
      }
    }
    
    min_i <- match(min(liks), liks)
  }	
  
  
  if((mod[1]=='d')&(mod[2]=='d')&(mod[3]=='d'))
  {
    result <- optim(betas[min_i,], loglik_AtCtEt_epsp_g, gr_AtCtEt_epsp_g, pheno_m = matrix(pheno_m), pheno_d = matrix(pheno_d), B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_c_m=B_des_c_m, B_des_c_d=B_des_c_d, B_des_e_m=B_des_e_m, B_des_e_d=B_des_e_d, var_b_a=vars[min_i,1], var_b_c=vars[min_i,2], var_b_e=vars[min_i,3], D_a=D_a, D_c=D_c, D_e=D_e, lower = c(rep(-12,n_a+n_c),rep(-10,n_e)), upper = c(rep(12,n_a+n_c),rep(10,n_e)), method = "L-BFGS-B", control=list(maxit = 3000), hessian = TRUE)
  }else{
    v_a_t <- vars[min_i,1]
    v_c_t <- vars[min_i,2]
    v_e_t <- vars[min_i,3]
    
    if(mod[1]!='d')
    {v_a_t <- lower}
    if(mod[2]!='d')
    {v_c_t <- lower}
    if(mod[3]!='d')
    {v_e_t <- lower}
    
    low_a <- (-1)*limit
      upp_a <- limit
      low_c <- (-1)*limit
      upp_c <- limit
      low_e <- rep((-1)*limit_e,n_e)
      upp_e <- rep(limit_e,n_e)
    
    result <- optim(betas[min_i,], loglik_AtCtEt_epsp_g, gr_AtCtEt_epsp_g, pheno_m = pheno_m, pheno_d = pheno_d, B_des_a_m=B_des_a_m, B_des_a_d=B_des_a_d, B_des_c_m=B_des_c_m, B_des_c_d=B_des_c_d, B_des_e_m=B_des_e_m, B_des_e_d=B_des_e_d, var_b_a=v_a_t, var_b_c=v_c_t, var_b_e=v_e_t, D_a=D_a, D_c=D_c, D_e=D_e, lower = c(rep(low_a,n_a),rep(low_c,n_c),low_e), upper = c(rep(upp_a,n_a),rep(upp_c,n_c),upp_e), method = "L-BFGS-B", control=list(maxit = 3000), hessian = TRUE)
  }
  beta_a <- result$par[1:n_a]
  beta_c <- result$par[(1+n_a):(n_a+n_c)]
  beta_e <- result$par[(1+n_a+n_c):(n_a+n_c+n_e)]
  
  lik_a <- lik_c <- lik_e <- 0
  if(mod[1]=='d')
  {
    D_a <- t(delta_a)%*%delta_a
    #beta_a <- ei_a$vectors%*%beta_a
    lik_a <- log(prod(ei_a$values[1:(length(ei_a$values)-2)]))
  }
  if(mod[2]=='d')
  {
    D_c <- t(delta_c)%*%delta_c
    #beta_c <- ei_c$vectors%*%beta_c
    lik_c <- log(prod(ei_c$values[1:(length(ei_c$values)-2)]))
  }
  if(mod[3]=='d')
  {
    D_e <- t(delta_e)%*%delta_e
    #beta_e <- ei_e$vectors%*%beta_e
    lik_e <- log(prod(ei_e$values[1:(length(ei_e$values)-2)]))
  }
  AtCtEtp_model <- list(D_a = D_a, D_c = D_c, D_e=D_e, pheno_m = pheno_m, pheno_d = pheno_d, T_m = T_m, T_d = T_d, knot_a=knots_a, knot_c=knots_c, knot_e=knots_e, beta_a=beta_a, beta_c=beta_c, beta_e=beta_e, con=result$convergence, lik=(min(liks)-lik_a-lik_c-lik_e)/2, iter=(liks)/2, var_b_a=vars[min_i,1], var_b_c=vars[min_i,2], var_b_e=vars[min_i,3], mod=mod, hessian = result$hessian)	
  
}


class(AtCtEtp_model) <- 'AtCtEtp_model'

#print('Estimates of beta_a:')
#print(beta_a)
#print('Estimates of beta_c:')
#print(beta_c)
#print('Estimates of beta_e:')
#print(beta_e)
#print(vars)
#print(betas)
return(invisible(AtCtEtp_model))

}