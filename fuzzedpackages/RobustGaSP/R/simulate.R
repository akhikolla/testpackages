##########################################################################
## simulate function
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
  
  
simulate.rgasp <- function (object, testing_input, num_sample=1,
                          testing_trend= matrix(1,dim(testing_input)[1],1),
                          r0=NA, rr0=NA,
                          sample_data=T,...){
  
  if(object@zero_mean=="Yes"){
    testing_trend=rep(0,dim(testing_input)[1]);
  }else{
    if( dim(testing_trend)[2]!=dim(object@X)[2]){
      stop("The dimensions of the design trend matrix and testing trend matrix do not match. \n")
    }
  }
  
  if( dim(testing_input)[2]!=dim(object@input)[2]){
    stop("The dimensions of the design matrix and testing inputs matrix do not match. \n")
  }
  
  
  num_testing_input <- dim(testing_input)[1]
  #X_testing = matrix(1,num_testing_input,1) ###testing trend
  
  
  testing_input=as.matrix(testing_input)
  
  p_x <- dim(object@input)[2]
  
  if(class(r0)=='logical'){
    if(!object@isotropic){
      r0 = as.list(1:object@p)
      for(i in 1:object@p){
        r0[[i]] = as.matrix(abs(outer(testing_input[,i], object@input[,i], "-")))
      }
    }else{
      r0 = as.list(1)
      if(p_x<object@num_obs){
        r0_here=0
        for(i in 1:p_x){
          r0_here=r0_here+(as.matrix(abs(outer(testing_input[,i], object@input[,i], "-"))))^2
        }
        r0[[1]]=sqrt(r0_here)
      }else{
        r0[[1]]=euclidean_distance(testing_input,object@input)
      }
    }
  }else if(class(r0)=='matrix'){
    r0_here=r0
    r0 = as.list(1)
    r0[[1]]=r0_here
  }else if(class(r0)!='list'){
    stop("r0 should be either a matrix or a list \n")
  }
  
  if(length(r0)!=object@p){
    stop("the number of R0 matrices should be the same 
         as the number of range parameters in the kernel \n")
  }
  if( (dim(r0[[1]])[1]!=num_testing_input) | (dim(r0[[1]])[2]!=object@num_obs)){
    stop("the dimension of R0 matrices should match the number of observations \n")
  }
  
  ###form matrix for rr0
  
  if(class(rr0)=='logical'){
    if(!object@isotropic){
      rr0 = as.list(1:object@p)
      for(i in 1:object@p){
        rr0[[i]] = as.matrix(abs(outer(testing_input[,i], testing_input[,i], "-")))
      }
    }else{
      rr0 = as.list(1)
      if(p_x<object@num_obs){
        r0_here=0
        for(i in 1:p_x){
          r0_here=r0_here+(as.matrix(abs(outer(testing_input[,i], testing_input[,i], "-"))))^2
        }
        rr0[[1]]=sqrt(r0_here)
      }else{
        rr0[[1]]=euclidean_distance(testing_input,testing_input)
      }
    }
  }else if(class(rr0)=='matrix'){
    rr0_here=rr0
    rr0 = as.list(1)
    rr0[[1]]=rr0_here
  }else if(class(rr0)!='list'){
    stop("r0 should be either a matrix or a list \n")
  }
  
  if(length(rr0)!=object@p){
    stop("the number of R0 matrices should be the same 
         as the number of range parameters in the kernel \n")
  }
  if( (dim(rr0[[1]])[1]!=num_testing_input) | (dim(rr0[[1]])[2]!=num_testing_input)){
    stop("the dimension of R0 matrices should match the number of observations \n")
  } 
  
  
  # rr0 = as.list(1:object@p)
  # for(i in 1:object@p){
  #   rr0[[i]] = as.matrix(abs(outer(testing_input[,i], testing_input[,i], "-")))
  # }
  
  ##change kernel type to integer to pass to C++ code
  kernel_type_num=rep(0,  object@p)
  for(i_p in 1:  object@p){
    if(object@kernel_type[i_p]=="matern_5_2"){
      kernel_type_num[i_p]=as.integer(3)
    }else if (object@kernel_type[i_p]=="matern_3_2"){
      kernel_type_num[i_p]=as.integer(2)
    }else if (object@kernel_type[i_p]=="pow_exp"){
      kernel_type_num[i_p]=as.integer(1)
    }else if (object@kernel_type[i_p]=="periodic_gauss"){  ##this is periodic folding on Gaussian kernel
      kernel_type_num[i_p]=as.integer(4)
    }else if (object@kernel_type[i_p]=="periodic_exp"){   ##this is periodic folding on Exponential kernel
      kernel_type_num[i_p]=as.integer(5)
    }
  }
  
  #####the following the posterior mean and cholesky decomposition of the sigma^2C_Star_star
  mean_cov_list=generate_predictive_mean_cov(object@beta_hat,object@nugget,object@input,object@X,object@zero_mean,object@output,
                       testing_input,testing_trend,object@L,object@LX,object@theta_hat,
                       object@sigma2_hat,rr0,r0,kernel_type_num,object@alpha,object@method,sample_data)
  
  predictive_sample=matrix(0,num_testing_input,num_sample)
  rnorm_sample=matrix(rnorm(num_testing_input*num_sample,0,1),num_testing_input,num_sample)
  
  if( (object@method=='post_mode') | (object@method=='mmle') ){
      df=object@num_obs-object@q
      chisq_sample=rchisq(num_sample,df=df)
      
      for(i_sample in 1:num_sample){
        predictive_sample[,i_sample]=mean_cov_list[[1]]+mean_cov_list[[2]]%*%rnorm_sample[,i_sample]*sqrt(df/chisq_sample[i_sample])
      }
  }else if(object@method=='mle'){
    for(i_sample in 1:num_sample){
      predictive_sample[,i_sample]=mean_cov_list[[1]]+mean_cov_list[[2]]%*%rnorm_sample[,i_sample]
    }
  }
    
  return(predictive_sample)
}

