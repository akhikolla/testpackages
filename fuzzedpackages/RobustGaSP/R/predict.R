##########################################################################
## prediction function
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

predict.rgasp <- function (object, testing_input, testing_trend= matrix(1,dim(as.matrix(testing_input))[1],1), 
                           r0=NA,
                           interval_data=T,outasS3 = T, ...){

  if (!outasS3) {
    output.pred <- new("predrgasp")
    output.pred@call <- match.call()
  } 
  #else warning('In order to have a more informative predic.rgasp output, set the parameter outasS3 = F.\n')
  testing_input=as.matrix(testing_input)
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
  ####let users specify the trend
    num_testing_input <- dim(testing_input)[1]
    #X_testing = matrix(1,num_testing_input,1) ###testing trend
      

      
    testing_input=as.matrix(testing_input)

    ##form the r0 matrix
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
    

    
    if( (object@method=='post_mode')|(object@method=='mmle') ){
      qt_025=qt(0.025,df=(object@num_obs-object@q))
      qt_975=qt(0.975,df=(object@num_obs-object@q))
      pred_list=pred_rgasp(object@beta_hat,object@nugget,object@input,object@X,object@zero_mean,object@output,
                           testing_input,testing_trend,object@L,object@LX,object@theta_hat,
                           object@sigma2_hat,qt_025,qt_975,r0,kernel_type_num,object@alpha,object@method,interval_data)
    }else if(object@method=='mle'){
      qn_025=qnorm(0.025,mean=0,sd=1)
      qn_975=qnorm(0.975,mean=0,sd=1)
      pred_list=pred_rgasp(object@beta_hat,object@nugget,object@input,object@X,object@zero_mean,object@output,
                           testing_input,testing_trend,object@L,object@LX,object@theta_hat,
                           object@sigma2_hat,qn_025,qn_975,r0,kernel_type_num,object@alpha,object@method,interval_data)
      # if(model@q>0){
      #    pred_list[[4]]=pred_list[[4]]/(model@num_obs-model@q)*(model@num_obs-model@q-2);
      # }else{
      #   pred_list[[4]]=pred_list[[4]]/(model@num_obs)*(model@num_obs-2);
      #   
      # }
    }
    output.list <- list()
    output.list$mean=pred_list[[1]]   #####can we all use @ or S? It will be more user friendly in that way 
    output.list$lower95=pred_list[[2]]
    output.list$upper95=pred_list[[3]]
    output.list$sd=sqrt(pred_list[[4]]) 

    if (!outasS3) {
      auxcall <- output.pred@call
      output.pred <- as.S4prediction.predict(output.list)
      output.pred@call <-auxcall
      output.list <- output.pred
    }
    
  return (output.list)
}

as.S4prediction.predict <- function(object){
  auxres <- new("predrgasp")
  auxres@call = match.call()
  auxres@mean = object$mean
  auxres@lower95 = object$lower95
  auxres@upper95 = object$upper95
  auxres@sd = object$sd
  return(auxres)
}



predict.ppgasp <- function (object, testing_input, testing_trend= matrix(1,dim(testing_input)[1],1), 
                            r0=NA,
                            interval_data=T,outasS3 = T, ...){

  if (!outasS3) {
    output.pred <- new("predppgasp")
    output.pred@call <- match.call()
  }
  
  testing_input=as.matrix(testing_input)
  
  #output.pred <- new("predppgasp")
  
  #else warning('In order to have a more informative predic.rgasp output, set the parameter outasS3 = F.\n')

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
  ####let users specify the trend
  num_testing_input <- dim(testing_input)[1]
  #X_testing = matrix(1,num_testing_input,1) ###testing trend

  qt_025=qt(0.025,df=(object@num_obs-object@q))
  qt_975=qt(0.975,df=(object@num_obs-object@q))


  testing_input=as.matrix(testing_input)
 
  ##form the r0 matrix
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
  

  #pred_list=pred_ppgasp(object@beta_hat,object@nugget,object@input,object@X,object@zero_mean,object@output,
  #                     testing_input,testing_trend,object@L,object@LX,object@theta_hat,
  #                     object@sigma2_hat,qt_025,qt_975,r0,kernel_type_num,object@alpha)

  if((object@method=='post_mode')|(object@method=='mmle')){
    qt_025=qt(0.025,df=(object@num_obs-object@q))
    qt_975=qt(0.975,df=(object@num_obs-object@q))
    pred_list=pred_ppgasp(object@beta_hat,object@nugget,object@input,object@X,object@zero_mean,object@output,
                          testing_input,testing_trend,object@L,object@LX,object@theta_hat,
                          object@sigma2_hat,qt_025,qt_975,r0,kernel_type_num,object@alpha,object@method,interval_data)
  }else if(object@method=='mle'){
    qn_025=qnorm(0.025,mean=0,sd=1)
    qn_975=qnorm(0.975,mean=0,sd=1)
    pred_list=pred_ppgasp(object@beta_hat,object@nugget,object@input,object@X,object@zero_mean,object@output,
                          testing_input,testing_trend,object@L,object@LX,object@theta_hat,
                          object@sigma2_hat,qn_025,qn_975,r0,kernel_type_num,object@alpha,object@method,interval_data)
  }
  
  output.list <- list()
  output.list$mean=pred_list[[1]]   #####can we all use @ or S? It will be more user friendly in that way
  output.list$lower95=pred_list[[2]]
  output.list$upper95=pred_list[[3]]
  output.list$sd=sqrt(pred_list[[4]])

  if (!outasS3) {
    auxcall <- output.pred@call
    output.pred <- as.S4prediction.predict_ppgasp(output.list)
    output.pred@call <-auxcall
    output.list <- output.pred
  }

  return (output.list)
}


as.S4prediction.predict_ppgasp <- function(object){
  auxres <- new("predpprgasp")
  auxres@call = match.call()
  auxres@mean = object$mean
  auxres@lower95 = object$lower95
  auxres@upper95 = object$upper95
  auxres@sd = object$sd
  return(auxres)
}

