##########################################################################
## Functions needed for the Robust GaSP
## 
## Robust GaSP Package
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, April 2013.
##
## Copyright (C) 2015-present Mengyang Gu, Jesus Palomo , James O. Berger
##							  
##    
##########################################################################

###############1 dim example function from Higdon 2002 
higdon.1.data <- function(s)
{
 # Implementation of the function in Higdon, D. (2002). 
 #  Space and space-time modeling using process convolutions. 
 #  In Quantitative methods for current environmental issues (pp. 37-56).
 #  Springer London.

 # First term: large scale variation
 # Second term: a fifth of the variation at four times the frequency
 return(sin(2*pi*s/10) + 0.2 * sin(2*pi*s/2.5))
}


#############2 dimension example function from Lim et al. 2002
limetal.2.data <- function(x)
{
  res1 <- 30 + 5*x[1]*sin(5*x[1])
  res2 <- 4 + exp(-5*x[2])
  
  return((res1*res2 - 100) / 6)
}

#############a 3 dimensional example function from Dette and Pepelyshev 2010
dettepepel.3.data <- function(x){

  res1 <- 4 * (x[1] - 2 + 8*x[2] - 8*x[2]^2)^2
  res2 <- (3 - 4*x[2])^2
  res3 <- 16 * sqrt(x[3]+1) * (2*x[3]-1)^2
  
  return(res1 + res2 + res3)
}

# The Borehole Function from Worley (1987)
borehole <- function(x)
{

  rw <- x[1]
  r  <- x[2]
  Tu <- x[3]
  Hu <- x[4]
  Tl <- x[5]
  Hl <- x[6]
  L  <- x[7]
  Kw <- x[8]
  
  res1 <- 2 * pi * Tu * (Hu-Hl)
  res2 <- log(r/rw)
  res3 <- 1 + (2*L*Tu / (res2*rw^2*Kw)) + (Tu / Tl)
  
  return(res1/res2/res3)
}

#############the environmental model with 4 dimensional input from Bliznyuk et al. (2008)
environ.4.data <- function(x, s=c(0.5, 1, 1.5, 2, 2.5), t=seq(from=0.3, to=60, by=0.3)){

  M   <- x[1]
  D   <- x[2]
  L   <- x[3]
  tau <- x[4]
  
  ds <- length(s)
  dt <- length(t)
  dY <- ds * dt
  Y <- matrix(0, ds, dt)
  
  # Create matrix Y, where each row corresponds to si and each column
  # corresponds to tj.
  for (ii in 1:ds) {
    si <- s[ii]
    for (jj in 1:dt) {
      tj <- t[jj]
      
      term1a <- M / sqrt(4*pi*D*tj)
      term1b <- exp(-si^2 / (4*D*tj))
      term1 <- term1a * term1b
      
      term2 <- 0
      if (tau < tj) {
        term2a <- M / sqrt(4*pi*D*(tj-tau))
        term2b <- exp(-(si-L)^2 / (4*D*(tj-tau)))
        term2 <- term2a * term2b
      }
      
      C <- term1 + term2
      Y[ii, jj] <- sqrt(4*pi) * C
    }
  }
  
  # Convert the matrix into a vector (by rows).
  Yrow <- t(Y)
  y <- t(as.vector(Yrow))
  return(y)
}


#############5 dimension example function from Friedman 1991
friedman.5.data <- function(x)
{
  return(10 * sin(pi*x[1]*x[2]) + 20 * (x[3]-0.5)^2 + 10*x[4] + 5*x[5])
}


####################
findInertInputs<-function(object,threshold=0.1){
  P_hat=object@p*object@beta_hat*object@CL/sum(object@beta_hat*object@CL)
  index_inert=which(P_hat<threshold)
  
  cat('The estimated normalized inverse range parameters are :', P_hat,'\n')
  if(length(which(P_hat<0.1))>0){
    cat('The inputs ', index_inert, 'are suspected to be inert inputs','\n')
  }else{
    cat('no input is suspected to be an inert input', '\n')
  }
  P_hat
}

####################
neg_log_marginal_post_approx_ref <- function(param,nugget, nugget.est,R0,X,zero_mean,output,CL,a,b,kernel_type,alpha) {
   #####this has mean X, we should also include the case where X is not zero
   #####
   lml=log_marginal_lik(param,nugget,nugget.est,R0,X,zero_mean,output,kernel_type,alpha);
   lp=log_approx_ref_prior(param,nugget,nugget.est,CL,a,b);
     
  -(lml+lp)

}

neg_log_profile_lik <- function(param,nugget, nugget.est,R0,X,zero_mean,output,CL,a,b,kernel_type,alpha) {
  #####this has mean X, we should also include the case where X is not zero
  #####
  lpl=log_profile_lik(param,nugget,nugget.est,R0,X,zero_mean,output,kernel_type,alpha);

  -lpl
  
}

neg_log_marginal_lik <- function(param,nugget, nugget.est,R0,X,zero_mean,output,CL,a,b,kernel_type,alpha) {
  #####this has mean X, we should also include the case where X is not zero
  #####
  lml=log_marginal_lik(param,nugget,nugget.est,R0,X,zero_mean,output,kernel_type,alpha);
  
  -lml
  
}


####################
neg_log_marginal_post_approx_ref_ppgasp <- function(param,nugget, nugget.est,R0,X,zero_mean,output,CL,a,b,kernel_type,alpha) {
  #####this has mean X, we should also include the case where X is not zero
  #####
  lml=log_marginal_lik_ppgasp(param,nugget,nugget.est,R0,X,zero_mean,output,kernel_type,alpha);
  lp=log_approx_ref_prior(param,nugget,nugget.est,CL,a,b);
  #print(param)
  #print(-(lml+lp))
  
  -(lml+lp)
  
}

neg_log_profile_lik_ppgasp <- function(param,nugget, nugget.est,R0,X,zero_mean,output,CL,a,b,kernel_type,alpha) {
  #####this has mean X, we should also include the case where X is not zero
  #####
  lpl=log_profile_lik_ppgasp(param,nugget,nugget.est,R0,X,zero_mean,output,kernel_type,alpha);
  #print(param)
  #print(-(lml+lp))
  
  -lpl
}

neg_log_marginal_lik_ppgasp <- function(param,nugget, nugget.est,R0,X,zero_mean,output,CL,a,b,kernel_type,alpha) {
  #####this has mean X, we should also include the case where X is not zero
  #####
  lml=log_marginal_lik_ppgasp(param,nugget,nugget.est,R0,X,zero_mean,output,kernel_type,alpha);
  #print(param)
  #print(-(lml+lp))
  
  -lml
}



####################
neg_log_marginal_post_ref<- function(param,nugget, nugget.est,R0,X,zero_mean,output,prior_choice,kernel_type,alpha) {
  
  lmp=log_ref_marginal_post(param,nugget,nugget.est,R0,X,zero_mean,output,kernel_type,alpha);
  
  
  if(prior_choice=='ref_xi'){
     if(nugget.est==T){###note that in this case nu is also needed to be transformed have tail properties
      #-sum(param)-lmp  ###this will let nu be non zero so we need to be careful
       -sum(param[1:(length(param))])-lmp ###this might be fine when intercept is in the trend
     }else{ ###no nugget
      -sum(param)-lmp
     }
  }else if(prior_choice=='ref_gamma'){
    if(nugget.est==T){###note that in this case nu is also needed to be transformed have tail properties
      -2*sum(param[1:(length(param)-1)])-lmp 
    }else{ ###no nugget
      -2*sum(param)-lmp
    }
  }
}

###ppgasp
neg_log_marginal_post_ref_ppgasp<- function(param,nugget, nugget.est,R0,X,zero_mean,output,prior_choice,kernel_type,alpha) {
  
  lmp=log_ref_marginal_post_ppgasp(param,nugget,nugget.est,R0,X,zero_mean,output,kernel_type,alpha);
  
  
  if(prior_choice=='ref_xi'){
    if(nugget.est==T){###note that in this case nu is also needed to be transformed have tail properties
      #-sum(param)-lmp  ###this will let nu be non zero so we need to be careful
      -sum(param[1:(length(param))])-lmp ###this might be fine when intercept is in the trend
    }else{ ###no nugget
      -sum(param)-lmp
    }
  }else if(prior_choice=='ref_gamma'){
    if(nugget.est==T){###note that in this case nu is also needed to be transformed have tail properties
      -2*sum(param[1:(length(param)-1)])-lmp 
    }else{ ###no nugget
      -2*sum(param)-lmp
    }
  }
}


########################
neg_log_marginal_post_approx_ref_deriv<- function(param,nugget,nugget.est,R0,X,zero_mean,output,CL,a,b,kernel_type,alpha) {
  lml_dev=log_marginal_lik_deriv(param,nugget,nugget.est,R0,X,zero_mean,output,kernel_type,alpha)
  lp_dev=log_approx_ref_prior_deriv(param,nugget,nugget.est,CL,a,b)
  
  -(lml_dev+lp_dev)*exp(param)
    
}


neg_log_profile_lik_deriv<- function(param,nugget,nugget.est,R0,X,zero_mean,output,CL,a,b,kernel_type,alpha) {
  lpl_dev=log_profile_lik_deriv(param,nugget,nugget.est,R0,X,zero_mean,output,kernel_type,alpha)

  -(lpl_dev)*exp(param)
  
}


neg_log_marginal_lik_deriv<- function(param,nugget,nugget.est,R0,X,zero_mean,output,CL,a,b,kernel_type,alpha) {
  lml_dev=log_marginal_lik_deriv(param,nugget,nugget.est,R0,X,zero_mean,output,kernel_type,alpha)
  
  -(lml_dev)*exp(param)
  
}


#######################ppgasp
neg_log_marginal_post_approx_ref_deriv_ppgasp<- function(param,nugget,nugget.est,R0,X,zero_mean,output,CL,a,b,kernel_type,alpha) {
  lml_dev=log_marginal_lik_deriv_ppgasp(param,nugget,nugget.est,R0,X,zero_mean,output,kernel_type,alpha)
  lp_dev=log_approx_ref_prior_deriv(param,nugget,nugget.est,CL,a,b)
  
  -(lml_dev+lp_dev)*exp(param)
  
}


neg_log_profile_lik_deriv_ppgasp<- function(param,nugget,nugget.est,R0,X,zero_mean,output,CL,a,b,kernel_type,alpha) {
  lpl_dev=log_profile_lik_deriv_ppgasp(param,nugget,nugget.est,R0,X,zero_mean,output,kernel_type,alpha)

  -(lpl_dev)*exp(param)
  
}

neg_log_marginal_lik_deriv_ppgasp<- function(param,nugget,nugget.est,R0,X,zero_mean,output,CL,a,b,kernel_type,alpha) {
  lml_dev=log_marginal_lik_deriv_ppgasp(param,nugget,nugget.est,R0,X,zero_mean,output,kernel_type,alpha)
  
  -(lml_dev)*exp(param)
  
}


#############this is a function to search the lower bounds for range parameters beta
#####need R0 in the function
search_LB_prob<-function(param, R0,COND_NUM_UB,p,kernel_type,alpha,nugget){
  num_obs=dim(R0[[1]])[1]
  propose_prob=exp(param)/(exp(param)+1)
  LB=NULL
  for( i_LB in 1:p){
    LB=c(LB, log(-log(propose_prob)/(max(R0[[i_LB]]))))    ###LB is log beta
  }
  
  R=separable_multi_kernel(R0,exp(LB),kernel_type,alpha)  
  
  # if(!isotropic){
  #   R=separable_multi_kernel(R0,exp(LB),kernel_type,alpha)  
  # }else{
  #   R=pow_exp_funct(R0[[1]],exp(LB),1)  
  # }
  R=as.matrix(R)
  R=R+nugget*diag(num_obs)
  
  (kappa(R)-COND_NUM_UB)^2
  ##one might change it to
  ##(kappa(R,exact=T)-COND_NUM_UB)^2
  
}

###leave one out 
leave_one_out_rgasp<-function(object){
  R_tilde=object@L%*%t(object@L)+object@nugget
  sigma_2=rep(0,object@num_obs)
  mean=rep(0,object@num_obs)
  if(object@zero_mean=="Yes"){
    for(i in 1:object@num_obs){
      r_sub=R_tilde[-i,i]
      L_sub=t(chol(R_tilde[-i,-i]))
      
      r_sub_t_R_sub_inv=t(backsolve(t(L_sub),forwardsolve(L_sub,r_sub)))
      
      R_sub_inv_y=backsolve(t(L_sub),forwardsolve(L_sub,object@output[-i] ))
      mean[i]=r_sub_t_R_sub_inv%*%object@output[-i]
      sigma_2_hat=object@output[-i]%*%R_sub_inv_y/(object@num_obs-1)
      sigma_2[i]=sigma_2_hat*(R_tilde[i,i]-r_sub_t_R_sub_inv%*%r_sub)
    }
  }else{
    for(i in 1:object@num_obs){
      r_sub=R_tilde[-i,i]
      L_sub=t(chol(R_tilde[-i,-i]))
      r_sub_t_R_sub_inv=t(backsolve(t(L_sub),forwardsolve(L_sub,r_sub)))
      
      R_inv_X=backsolve(t(L_sub),forwardsolve(L_sub,(object@X[-i,]) ))

        
      L_X=t(chol(t(object@X[-i,])%*%R_inv_X))
      theta_hat=backsolve(t(L_X),forwardsolve(L_X,t(R_inv_X)%*%object@output[-i]))

      tilde_output=object@output[-i]-object@X[-i,]%*%theta_hat
      mean[i]=object@X[i,]%*%theta_hat+r_sub_t_R_sub_inv%*%tilde_output
      
      
      if( (object@method=='post_mode') |(object@method=='mmle') ){
        sigma2_hat=t(tilde_output)%*%backsolve(t(L_sub),forwardsolve(L_sub,tilde_output ))/(object@num_obs-1-object@q)
        
        c_star=(R_tilde[i,i]-r_sub_t_R_sub_inv%*%r_sub)
        
  
        h_hat=object@X[i,]-t(object@X[-i,])%*%t(r_sub_t_R_sub_inv)
        c_star_star= c_star+ t(h_hat)%*%backsolve(t(L_X),forwardsolve(L_X,h_hat))
        sigma_2[i]=sigma2_hat*(c_star_star)
      }else if(object@method=='mle'){
        sigma2_hat=t(tilde_output)%*%backsolve(t(L_sub),forwardsolve(L_sub,tilde_output ))/(object@num_obs-1)
        
        c_star=(R_tilde[i,i]-r_sub_t_R_sub_inv%*%r_sub)
        
        sigma_2[i]=sigma2_hat*(c_star)
      }
      #sigma_2[i]=object@sigma2_hat*(c_star)
      
    }
  }
  
  return(list(mean=mean, sd=sqrt(sigma_2)))
  
  #plot((output-mean)/sqrt(sigma_2))

}