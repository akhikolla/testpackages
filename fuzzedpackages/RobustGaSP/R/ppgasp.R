##########################################################################
## rgasp fit function
## 
## Robust GaSP Package
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, April 2013.
##
## Copyright (C) 2015-present Mengyang Gu, Jesus Palomo, James O. Berger
##							  
##    
##########################################################################

##function for parallel partial Gaussian stochastic process
ppgasp <- function(design, response,trend=matrix(1,dim(response)[1],1),zero.mean="No",nugget=0,
                  nugget.est=F,range.par=NA,method='post_mode',prior_choice='ref_approx',a=0.2,
                  b=1/(length(response))^{1/dim(as.matrix(design))[2]}*(a+dim(as.matrix(design))[2]),
                  kernel_type='matern_5_2',isotropic=F,R0=NA,optimization='lbfgs',
                  alpha=rep(1.9,dim(as.matrix(design))[2]),lower_bound=T,max_eval=max(30,20+5*dim(design)[2]),
                  initial_values=NA,num_initial_values=2){
  

  if (!is.logical(nugget.est) && length(nugget.est) != 1){
    stop("nugget.est should be boolean (either T or F) \n")
  }
  
  if(nugget!=0 & nugget.est==T){
    stop("one cannot fix and estimate the nugget at the same time \n")  
  }
  
  if(!is.na(range.par)){  
    if(length(range.par)!=dim(as.matrix(design))[2]){
      stop("range.par should either be fixed or estimated.")    
    }
    if(nugget.est){
      stop("We do not support fixing range parameters while estimating the nugget.")      
    }
  }
  
  if(!is.numeric(nugget)){
    stop("nugget should be a numerical value \n")  
  }
  
  if( (optimization!='lbfgs') & (optimization!='nelder-mead')& (optimization!='brent')){
    stop("optimization should be 'lbfgs' or 'nelder-mead' or 'brent' \n")  
    
  }
  
  #model <- new("rgasp")
  model <- new("ppgasp")
  
  model@call <- match.call()
  
  
  #cat('Kernel type is ',model@kernel_type,'\n')
  
  #cat('Multiple starting points: ',multiple_starts, '\n')

  #####Only numeric inputs are allowed
  design=as.matrix(design)
  model@input <- matrix(as.numeric(design), dim(design)[1],dim(design)[2])
  
  # print(model@input)
  model@output <- as.matrix(response)
  
  
  if(dim(model@output)[1]!=dim(model@input)[1]){
    stop("the number of rows in the input should be the same as the number of rows 
         in the output \n")  
    
  }
  
  ##response is n x k matrix, where n is the number of obs and k is the number of loc
  if (zero.mean=="Yes"){
    trend=matrix(0,dim(model@output)[1],1)
  }
  
  # print(model@output)
  p_x <- dim(model@input)[2]
  kernel_type=c(kernel_type) 
  
  if(!is.logical(isotropic)){
    stop("isotropic should be either true or false \n")  
  }
  model@isotropic=isotropic
  
  
  if(!model@isotropic){
    model@alpha <-alpha
  }else{
    model@alpha=alpha[1]
  }
  
  if(model@isotropic){
    model@p =as.integer(1)
  }else{
    model@p=p_x
  }
  
  if(optimization=='brent' & model@p!=1){
    stop('The Brent optimization method can only work for optimizing one parameter \n')
  }
  
  model@num_obs <- dim(model@output)[1]
  model@k <- dim(model@output)[2]
  
  ## Checking the dimensions between the different inputs
  if (model@num_obs != dim(model@output)[1]){
    stop("The dimensions of the design matrix and the response do not match. \n")
  }
  #  stop("The dimension of the training points and simulator values are different. \n")
  
  const_column=test_const_column(model@output);
  
  if(const_column){
    stop("Please delete the column of the response matrix that has the same value across all rows. \n")
  }
  ###add method 
  if( (method!='post_mode') & (method!='mle') & (method!='mmle') ){
    stop("The method should be post_mode or mle or mmle. \n")
  }
  model@method=method
  
  if(optimization=='lbfgs'& ((model@num_obs*model@k)>=10000) ){
    warning("Please consider to select optimization= 'nelder-mead' or 'brent' 
              as the derivative-based optimization can be inaccurate 
             when the number of observations is large.")
  }
    
  ##kernel type
  if(!model@isotropic){
    if(length(kernel_type)==1){
      model@kernel_type=rep(kernel_type,model@p)
      # if(kernel_type!='isotropic'){
      #   kernel_type=rep(kernel_type,model@p)
      # }else{
      #   kernel_type=c('isotropic')
      # }
    }else if(length(kernel_type)!=model@p){
      stop("Please specify the correct number of kernels. \n")
    }else{
      model@kernel_type=kernel_type
    }
  }else{
    model@kernel_type=kernel_type
  }

  
  
  ##change kernel type to integer to pass to C++ code
  ##1 is power exponent, 2 is matern with roughness 3/2, and 3 is matern with roughenss parameter 5/2
  kernel_type_num=rep(0,  model@p)
  for(i_p in 1:  model@p){
    if(model@kernel_type[i_p]=="matern_5_2"){
      kernel_type_num[i_p]=as.integer(3)
    }else if (model@kernel_type[i_p]=="matern_3_2"){
      kernel_type_num[i_p]=as.integer(2)
    }else if (model@kernel_type[i_p]=="pow_exp"){
      kernel_type_num[i_p]=as.integer(1)
    }else if (model@kernel_type[i_p]=="periodic_gauss"){  ##this is periodic folding on Gaussian kernel
      kernel_type_num[i_p]=as.integer(4)
    }else if (model@kernel_type[i_p]=="periodic_exp"){   ##this is periodic folding on Exponential kernel
      kernel_type_num[i_p]=as.integer(5)
    }
  }
  
  

  #####I now build the gasp emulator here
  ##default setting constant mean basis. Of course We should let user to specify it
  #model@X = matrix(1,model@num_obs,1)   ####constant mean
  model@X=trend;               ###If the trend is specified, use it. If not, use a constant mean. 
  model@zero_mean=zero.mean
  #######number of regressor
  if(model@zero_mean=="Yes"){
    model@q=as.integer(0)
  }else{
    model@q = dim(model@X)[2]; # NOTE THIS IS ALWAYS 1 SINCE YOU DEFINE IT THAT WAY ABOVE
  }
  ####################correlation matrix
  
  model@nugget.est <- nugget.est
  if(class(R0)=='logical'){
    if(!model@isotropic){
      model@R0 = as.list(1:model@p)
      for(i in 1:model@p){
        model@R0[[i]] = as.matrix(abs(outer(model@input[,i], model@input[,i], "-")))
      }
    }else{
      model@R0 = as.list(1)
      if(p_x<model@num_obs){
        R0_here=0
        for(i in 1:p_x){
          R0_here=R0_here+(as.matrix(abs(outer(model@input[,i], model@input[,i], "-"))))^2
        }
        model@R0[[1]]=sqrt(R0_here)
      }else{
        model@R0[[1]]=euclidean_distance(model@input,model@input)
      }
    }
  }else if(class(R0)=='matrix'){
    model@R0 = as.list(1)
    model@R0[[1]]=R0
  }else if(class(R0)=='list'){
    model@R0=R0
  }else{
    stop("R0 should be either a matrix or a list \n")
  }
  
  ##check for R0
  if(length(model@R0)!=model@p){
    stop("the number of R0 matrices should be the same 
         as the number of range parameters in the kernel \n")
  }
  if( (dim(model@R0[[1]])[1]!=model@num_obs) | (dim(model@R0[[1]])[2]!=model@num_obs)){
    stop("the dimension of R0 matrices should match the number of observations \n")
  }
  
  
  ###########calculating lower bound for beta
  model@CL = rep(0,model@p)    ###CL is also used in the prior so I make it a model parameter
  
  if(!model@isotropic){
    for(i_cl in 1:model@p){
      #  model@CL[i_cl] = mean(model@R0[[i_cl]][which(model@R0[[i_cl]]>0)])
      
      model@CL[i_cl] = (max(model@input[,i_cl])-min(model@input[,i_cl]))/model@num_obs^{1/model@p}
    }
  }else{
    model@CL[1]=max(model@R0[[1]])/model@num_obs
  }
  
  
  if(is.na(range.par)){
    ########this also depends on the kernel
    ##get some initial values of beta that we can start from a good set of range parameters to optimize
    
   
      COND_NUM_UB = 10^{16}  ###maximum condition number, this might be a little too large
      
      if(lower_bound==T){
      LB_all = optimize(search_LB_prob, interval=c(-5,12), maximum = FALSE, R0=model@R0,COND_NUM_UB= COND_NUM_UB,
                        p=model@p,kernel_type=kernel_type_num,alpha=model@alpha,nugget=nugget) ###find a lower bound for parameter beta
      
      LB_prob = exp(LB_all$minimum)/(exp(LB_all$minimum)+1)
      
      LB = NULL
      
      for( i_LB in 1:model@p){
        LB = c(LB, log(-log(LB_prob)/(max(model@R0[[i_LB]]))))    ###LB is lower bound for log beta, may consider to have it related to p
      }
     }else{
        ##give some empirical bound that passes the initial values if lower bound is F
        ##could can change the first initial value if not search the bound
        LB = NULL
        
        for( i_LB in 1:model@p){
          LB = c(LB, -log(0.1)/((max(model@input[,i_LB])-min(model@input[,i_LB]))*model@p))   
        }
        
      }
   
    if(lower_bound==T){
      if(model@nugget.est){
        model@LB=c(LB,-Inf)   
      }else{
        model@LB=LB
      }
    }else{
      if(model@nugget.est){
        model@LB=rep(-Inf,model@p+1)
      }else{
        model@LB=rep(-Inf,model@p)
      }
    }
    
    cat('The upper bounds of the range parameters are',1/exp(model@LB),'\n')
    
    ############################the lower bound might be needed to discuss
    
    #beta_initial=(a+model@p)/(model@p*model@CL*b)/2  ###half the prior value
    
    if(is.na(initial_values)[1]==T){
      beta_initial=matrix(0,num_initial_values,model@p)
      eta_initial=rep(0,num_initial_values)
      beta_initial[1,]=50*exp(LB) #####one start
      eta_initial[1]=0.0001
      if(num_initial_values>1){
        beta_initial[2,]=(a+model@p)/(model@p*model@CL*b)/2  ###half the prior value
        eta_initial[2]=0.0002

      }
      if(num_initial_values>2){
        for(i_ini in 3:num_initial_values){
          set.seed(i_ini)
          
          beta_initial[i_ini,]=10^3*runif(model@p)/model@CL
          eta_initial[i_ini]=10^(-3)*runif(1)
        }
      }
      initial_values=cbind(log(beta_initial),log(eta_initial))
    }
      
    if(method=='post_mode'){
      object_funct_name='marginal posterior'
    }else if(method=='mle'){
      object_funct_name='profile likelihood'
    }else{
      object_funct_name='marginal likelihood'
    }
    
    model@log_post=-Inf;
    if(optimization=='lbfgs'){
      for(i_ini in 1:num_initial_values){
          if(model@nugget.est){
            ini_value=initial_values[i_ini,]
          }else{
            ini_value=initial_values[i_ini,1:model@p]
            ###without the nugget
          }
          cat('The initial values of range parameters are', 1/exp(ini_value[1:model@p]),'\n')
          cat('Start of the optimization ', i_ini,' : \n')
          if(method=='post_mode'){
            if(prior_choice=='ref_approx'){####this one can be with nugget or without the nugget
              #  if (requireNamespace("lbfgs", quietly = TRUE)) {
              tt_all <- try(nloptr::lbfgs(ini_value, neg_log_marginal_post_approx_ref_ppgasp, 
                                          neg_log_marginal_post_approx_ref_deriv_ppgasp,nugget=nugget, nugget.est=model@nugget.est, 
                                          R0=model@R0,X=model@X, zero_mean=model@zero_mean,output=model@output, CL=model@CL, a=a,b=b,
                                          kernel_type=kernel_type_num,alpha=model@alpha,lower=model@LB,
                                          nl.info = FALSE, control = list(maxeval=max_eval)),TRUE)
              #   }
            }else if(prior_choice=='ref_xi'|prior_choice=='ref_gamma'){####this needs to be revised
              #  if (requireNamespace("lbfgs", quietly = TRUE)) {
              tt_all <- try(nloptr::lbfgs(ini_value, neg_log_marginal_post_ref_ppgasp, 
                                          nugget=nugget, nugget.est=nugget.est, R0=model@R0,
                                          X=model@X, zero_mean=model@zero_mean,output=model@output, prior_choice=prior_choice, kernel_type=kernel_type_num,
                                          alpha=model@alpha,lower=model@LB,nl.info = FALSE, control = list(maxeval=max_eval)),TRUE)
              # }
            }
          }else if(method=='mle'){
            tt_all <- try(nloptr::lbfgs(ini_value, neg_log_profile_lik_ppgasp, 
                                        neg_log_profile_lik_deriv_ppgasp,nugget=nugget, nugget.est=model@nugget.est, 
                                        R0=model@R0,X=model@X, zero_mean=model@zero_mean,output=model@output, CL=model@CL, a=a,b=b,
                                        kernel_type=kernel_type_num,alpha=model@alpha,lower=model@LB,
                                        nl.info = FALSE, control = list(maxeval=max_eval)),TRUE)
            
          }else if(method=='mmle'){
            tt_all <- try(nloptr::lbfgs(ini_value, neg_log_marginal_lik_ppgasp, 
                                        neg_log_marginal_lik_deriv_ppgasp,nugget=nugget, nugget.est=model@nugget.est, 
                                        R0=model@R0,X=model@X, zero_mean=model@zero_mean,output=model@output, CL=model@CL, a=a,b=b,
                                        kernel_type=kernel_type_num,alpha=model@alpha,lower=model@LB,
                                        nl.info = FALSE, control = list(maxeval=max_eval)),TRUE)
            
          }
          #if(class(tt_all)=="try-error"){
            #sink()
          #  stop(tt_all)
          #}
  
          if(class(tt_all)!="try-error"){
            if(model@nugget.est==F){
              nugget_par=nugget
            }else{
              nugget_par=exp(tt_all$par)[model@p+1]
            }
            
            if (tt_all$convergence>0){convergence=T}else{convergence=F}
            cat('The number of iterations is ', tt_all$iter,'\n',
                paste('The value of the ', object_funct_name, ' function is '), -tt_all$value,'\n',
                'Optimized range parameters are', 1/exp(tt_all$par)[1:model@p],'\n',
                'Optimized nugget parameter is', nugget_par,'\n',
                'Convergence: ', convergence,'\n' )
          
  
            if( (-tt_all$value)>model@log_post){
              #log_lik=-tt_all$value
              model@log_post=-tt_all$value
              model@nugget=nugget;
              if(nugget.est){
                model@beta_hat = exp(tt_all$par)[1:model@p];
                model@nugget=exp(tt_all$par)[model@p+1];
              }else{
                model@beta_hat = exp(tt_all$par);
                #  model@nugget=0;
              }
            }
          }
      }
    }else if(optimization=='nelder-mead'){
      for(i_ini in 1:num_initial_values){
        if(model@nugget.est){
          ini_value=initial_values[i_ini,]
        }else{
          ini_value=initial_values[i_ini,1:model@p]
          ###without the nugget
        }
        cat('The initial values of range parameters are', 1/exp(ini_value[1:model@p]),'\n')
        cat('Start of the optimization ', i_ini,' : \n')
        if(method=='post_mode'){
          if(prior_choice=='ref_approx'){####this one can be with nugget or without the nugget
        
            tt_all <- try(optim(ini_value, neg_log_marginal_post_approx_ref_ppgasp,nugget=nugget, nugget.est=model@nugget.est, 
                                        R0=model@R0,X=model@X, zero_mean=model@zero_mean,output=model@output, CL=model@CL, a=a,b=b,
                                        kernel_type=kernel_type_num,alpha=model@alpha),TRUE)
      
          }else if(prior_choice=='ref_xi'|prior_choice=='ref_gamma'){####this needs to be revised
           

            tt_all <- try(optim(ini_value, neg_log_marginal_post_ref_ppgasp, 
                                        nugget=nugget, nugget.est=nugget.est, R0=model@R0,
                                        X=model@X, zero_mean=model@zero_mean,output=model@output, prior_choice=prior_choice, kernel_type=kernel_type_num,
                                        alpha=model@alpha),TRUE)
            
          }
        }else if(method=='mle'){
          tt_all <- try(optim(ini_value, neg_log_profile_lik_ppgasp,nugget=nugget, nugget.est=model@nugget.est, 
                                      R0=model@R0,X=model@X, zero_mean=model@zero_mean,output=model@output, CL=model@CL, a=a,b=b,
                                      kernel_type=kernel_type_num,alpha=model@alpha),TRUE)
          
        }else if(method=='mmle'){
          tt_all <- try(optim(ini_value, neg_log_marginal_lik_ppgasp,nugget=nugget, nugget.est=model@nugget.est, 
                                      R0=model@R0,X=model@X, zero_mean=model@zero_mean,output=model@output, CL=model@CL, a=a,b=b,
                                      kernel_type=kernel_type_num,alpha=model@alpha),TRUE)
          
        }

        if(class(tt_all)!="try-error"){
          if(model@nugget.est==F){
            nugget_par=nugget
          }else{
            nugget_par=exp(tt_all$par)[model@p+1]
          }
          
          if (tt_all$convergence>0){convergence=T}else{convergence=F}
          cat('The number of iterations is ', tt_all$iter,'\n',
              paste('The value of the ', object_funct_name, ' function is '), -tt_all$value,'\n',
              'Optimized range parameters are', 1/exp(tt_all$par)[1:model@p],'\n',
              'Optimized nugget parameter is', nugget_par,'\n',
              'Convergence: ', convergence,'\n' )
          
          
          if( (-tt_all$value)>model@log_post){
            #log_lik=-tt_all$value
            model@log_post=-tt_all$value
            model@nugget=nugget;
            if(nugget.est){
              model@beta_hat = exp(tt_all$par)[1:model@p];
              model@nugget=exp(tt_all$par)[model@p+1];
            }else{
              model@beta_hat = exp(tt_all$par);
              #  model@nugget=0;
            }
          }
        }
      }
      
    }else if(optimization=='brent'){
      #for(i_ini in 1:num_initial_values){
        # if(model@nugget.est){
        #   ini_value=initial_values[i_ini,]
        # }else{
        #   ini_value=initial_values[i_ini,1:model@p]
        #   ###without the nugget
        # }
        #cat('The initial values of range parameters are', 1/exp(ini_value[1:model@p]),'\n')
        #cat('Start of the optimization ', i_ini,' : \n')
        if(method=='post_mode'){
          UB=max(20, 1/log(max(model@R0[[1]])))
          if(prior_choice=='ref_approx'){####this one can be with nugget or without the nugget
            
            tt_all <- try(optimize(neg_log_marginal_post_approx_ref_ppgasp,nugget=nugget, nugget.est=model@nugget.est, 
                                R0=model@R0,X=model@X, zero_mean=model@zero_mean,output=model@output, CL=model@CL, a=a,b=b,
                                kernel_type=kernel_type_num,alpha=model@alpha,lower=LB,upper=UB),TRUE)
            
          }else if(prior_choice=='ref_xi'|prior_choice=='ref_gamma'){####this needs to be revised
            
            
            tt_all <- try(optimize(neg_log_marginal_post_ref_ppgasp, 
                                nugget=nugget, nugget.est=nugget.est, R0=model@R0,
                                X=model@X, zero_mean=model@zero_mean,output=model@output, prior_choice=prior_choice, kernel_type=kernel_type_num,
                                alpha=model@alpha,lower=LB,upper=UB),TRUE)
            
          }
        }else if(method=='mle'){
          tt_all <- try(optimize(neg_log_profile_lik_ppgasp,nugget=nugget, nugget.est=model@nugget.est, 
                              R0=model@R0,X=model@X, zero_mean=model@zero_mean,output=model@output, CL=model@CL, a=a,b=b,
                              kernel_type=kernel_type_num,alpha=model@alpha,lower=LB,upper=UB),TRUE)
          
        }else if(method=='mmle'){
          tt_all <- try(optimize(neg_log_marginal_lik_ppgasp,nugget=nugget, nugget.est=model@nugget.est, 
                              R0=model@R0,X=model@X, zero_mean=model@zero_mean,output=model@output, CL=model@CL, a=a,b=b,
                              kernel_type=kernel_type_num,alpha=model@alpha,lower=LB,upper=UB),TRUE)
          
        }
        
        if(class(tt_all)!="try-error"){
          # if(model@nugget.est==F){
          #   nugget_par=nugget
          # }else{
          #   nugget_par=exp(tt_all$par)[model@p+1]
          # }
          
          #if (tt_all$convergence>0){convergence=T}else{convergence=F}
          cat( paste('The value of the ', object_funct_name, ' function is '), -tt_all$objective,'\n',
              'Optimized range parameters are', 1/exp(tt_all$minimum),'\n')
          

          
            #log_lik=-tt_all$objective
            model@log_post=-tt_all$objective
            model@beta_hat=exp(tt_all$minimum)
            model@nugget=nugget;
            # if(nugget.est){
            #   model@beta_hat = exp(tt_all$par)[1:model@p];
            #   model@nugget=exp(tt_all$par)[model@p+1];
            # }else{
            #   model@beta_hat = exp(tt_all$par);
            #   #  model@nugget=0;
            # }
          }
        
      #}
      
    }
    
  }else{###this is the case where the range parameters and the nugget are all fixed
    model@LB=rep(-Inf,model@p)
    model@beta_hat=1/range.par
    model@nugget=nugget
  }
  
  list_return=construct_ppgasp(model@beta_hat, model@nugget, model@R0, model@X, zero_mean=model@zero_mean,
                              model@output,kernel_type_num,model@alpha); 
  model@L=list_return[[1]];
  model@LX=list_return[[2]];
  model@theta_hat=list_return[[3]];
  #model@sigma2_hat=list_return[[4]];
  if( (method=='post_mode') | (method=='mmle') ){
    model@sigma2_hat=list_return[[4]];
  }else if (method=='mle'){
    if(model@q>0){
      model@sigma2_hat=list_return[[4]]*(model@num_obs-model@q)/model@num_obs;
    }
  }
  
  return(model)
}

show.ppgasp <- function(object) {	
  cat("\n")
  cat("Call:\n")
  print(object@call)
  

  cat('The dimension  of the design is: ',dim(object@input),'\n')
  cat('The dimension  of the output is: ',dim(object@output),'\n')


  cat('Range parameters: ', 1/object@beta_hat,'\n')
  cat('Nugget parameter: ', object@nugget,'\n')
  
  
  # cat('The slots of this object are: ',slotNames(object),'\n')
  
}