
Generalized_Testing_coefficient_of_mediator=function(G,mediator,outcome,covariates=NULL,outcome_type,approxi=TRUE,verbose=FALSE){
  num_outcome=length(outcome_type)
  
  Risk_set_VecSum=function(i,Vec){
    
    return(sum(Vec[1:i]))
  }
  
  
  compute_log_likelihood=function(theta_hat,gamma_tilde_mat,phi_vec,fixed_intercept_covariate,sigma2_gamma_vec){
    
    W_U_list=list()
    log_LL_1=rep(NA,num_outcome)
    
    d=1
    if (outcome_type[d]=="continuous"){
      if (is.null(theta_hat)){
        fixedeffectsdat=cbind(intercept_covariate)
        fixedeffects=c(fixed_intercept_covariate[,d])
      }else{
        fixedeffectsdat=cbind(intercept_covariate,mediator)
        fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
      }
      
      gamma_tilde=gamma_tilde_mat[,d]
      sigma2_gamma=sigma2_gamma_vec[d]
      phi=phi_vec[d]
      
      eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
      mu=eta
      w=as.matrix(rep(1/phi,nrow(G)))
      W_inv=diag(1/w[,1])
      W=diag(w[,1])
      W_U_list[[d]]=W
      Q_inv=eigenMapMatMult(eigenMapMatMult(t(G),W),G)+sigma2_gamma^(-1)*diag(ncol(G))
      #chol_Q_inv=chol(Q_inv)
      log_det_Q_inv=as.numeric(determinant(Q_inv, logarithm = TRUE)$modulus)
      #log_det_Q_inv=2*sum(log(diag(chol_Q_inv)))
      
      
      log_LL_1[d]=-0.5*ncol(G)*log(sigma2_gamma)-0.5*log_det_Q_inv+as.numeric( -0.5*(1/phi)*sum( (outcome-mu)^2 ) -0.5*nrow(G)*log(phi) -0.5*(sigma2_gamma)^(-1)*t(gamma_tilde)%*%gamma_tilde )
      
    }
    
    if (outcome_type[d]=="binary"){
      if (is.null(theta_hat)){
        fixedeffectsdat=cbind(intercept_covariate)
        fixedeffects=c(fixed_intercept_covariate[,d])
      }else{
        fixedeffectsdat=cbind(intercept_covariate,mediator)
        fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
      }
      
      gamma_tilde=gamma_tilde_mat[,d]
      sigma2_gamma=sigma2_gamma_vec[d]
      phi=phi_vec[d]
      
      eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
      mu=1/(1+exp(-eta))
      w=mu*(1-mu)
      W_inv=diag(1/w[,1])
      W=diag(w[,1])
      W_U_list[[d]]=W
      Q_inv=eigenMapMatMult(eigenMapMatMult(t(G),W),G)+sigma2_gamma^(-1)*diag(ncol(G))
      chol_Q_inv=chol(Q_inv)
      log_det_Q_inv=2*sum(log(diag(chol_Q_inv)))
      log_LL_1[d]=-0.5*ncol(G)*log(sigma2_gamma)-0.5*log_det_Q_inv+as.numeric(sum(outcome*eta-log(1+exp(eta)))-0.5*(sigma2_gamma)^(-1)*t(gamma_tilde)%*%gamma_tilde )  
      
    }
    
    if (outcome_type[d]=="count"){
      if (is.null(theta_hat)){
        fixedeffectsdat=cbind(intercept_covariate)
        fixedeffects=c(fixed_intercept_covariate[,d])
      }else{
        fixedeffectsdat=cbind(intercept_covariate,mediator)
        fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
      }
      
      gamma_tilde=gamma_tilde_mat[,d]
      sigma2_gamma=sigma2_gamma_vec[d]
      phi=phi_vec[d]
      
      eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
      mu=exp(eta)
      xi=exp(eta)/(phi+exp(eta))
      u_matrix=diag( as.numeric((outcome+phi)*xi*(1-xi)) )
      W_U_list[[d]]=u_matrix
      Q_inv=eigenMapMatMult(eigenMapMatMult(t(G),u_matrix),G)+sigma2_gamma^(-1)*diag(ncol(G))
      
      chol_Q_inv=chol(Q_inv)
      
      log_det_Q_inv=2*sum(log(diag(chol_Q_inv)))
      
      #log_LL=-0.5*ncol(G)*log(sigma2_gamma)-0.5*log_det_Q_inv+as.numeric( sum( outcome*eta-(outcome+phi)*log(phi+exp(eta))+phi*log(phi)+log(gamma(outcome+phi)/(gamma(phi)*factorial(outcome)) )) -0.5*(sigma2_gamma)^(-1)*t(gamma_tilde)%*%gamma_tilde )       
      log_LL_1[d]=-0.5*ncol(G)*log(sigma2_gamma)-0.5*log_det_Q_inv+as.numeric( sum( outcome*eta-(outcome+phi)*log(phi+exp(eta))+phi*log(phi)+lgamma(outcome+phi)-lgamma(phi)-lfactorial(outcome)  ) -0.5*(sigma2_gamma)^(-1)*t(gamma_tilde)%*%gamma_tilde )  
      
    }
    
    if (outcome_type[d]=="survival"){
      if (is.null(theta_hat)){
        fixedeffectsdat=cbind(intercept_covariate_reorder)
        fixedeffects=c(fixed_intercept_covariate[,d])
      }else{
        fixedeffectsdat=cbind(intercept_covariate_reorder,mediator_reorder)
        fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
      }
      if (is.null(fixedeffectsdat)){
        fixedeffectsdat=matrix(0,nrow=nrow(G)); fixedeffects=0
      }
      gamma_tilde=gamma_tilde_mat[,d]
      sigma2_gamma=sigma2_gamma_vec[d]
      
      eta=fixedeffectsdat%*%fixedeffects+G_reorder%*%gamma_tilde
      exp_eta=exp(eta)
      exp_eta_risk=sapply(1:nrow(G_reorder), function(i) Risk_set_VecSum(i,exp_eta))
      
      mat=lapply(2:nrow(G_reorder), function(i) status_reorder[i]*cov.wt(G_reorder[1:i,],exp_eta[1:i])$cov )
      weighted_covariance=Reduce(`+`,mat)
      Q_inv=weighted_covariance+sigma2_gamma^(-1)*diag(ncol(G_reorder))
      
      chol_Q_inv=chol(Q_inv)
      
      log_det_Q_inv=2*sum(log(diag(chol_Q_inv)))
      
      
      
      #log partial likelihood
      log_PL=sum(log ( (exp_eta/exp_eta_risk)^(status_reorder) ) )
      
      #SecondDeri=diag(apply(G,2,function(g) sum(status*sapply(1:nrow(G), function(i) -Risk_set_VecSum(i,(g^2)*exp_eta)+Risk_set_VecSum(i,g*exp_eta)^2 )/(exp_eta_risk^2)  ) ) )
      #SecondDeri=diag(-colSums(apply(G2_exp_eta_risk, 2, function(x) cumsum(x)*status_reorder/exp_eta_risk^2 ))+colSums(apply(G_exp_eta_risk, 2, function(x) cumsum(x)^2*status_reorder/exp_eta_risk^2 )))
      
      
      log_LL_1[d]=-0.5*ncol(G_reorder)*log(sigma2_gamma)+log_PL -0.5*(sigma2_gamma)^(-1)*t(gamma_tilde)%*%gamma_tilde -0.5*log_det_Q_inv
      
    }
    
    
    return(log_LL_1)
  }
  
  weighted.var = function(x, w, na.rm = FALSE) {
    if (na.rm) {
      w <- w[i <- !is.na(x)]
      x <- x[i]
    }
    sum.w <- sum(w)
    sum.w2 <- sum(w^2)
    mean.w <- sum(x * w) / sum(w)
    (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
                                         na.rm)
  }
  
  
  num_consecutive=10
  num_gap=10
  thres=1e-6
  
  
  # scale #
  G = scale(G,scale = F)
  mediator=scale(mediator,scale = F)
  d=1
  if (outcome_type[d]=="survival"){
    if (!is.null(covariates)){
      covariates=scale(covariates,scale = F)
      intercept_covariate=covariates
    }else{
      intercept_covariate=NULL
    }
  }else{
    if (!is.null(covariates)){
      covariates=scale(covariates,scale = F)
      intercept_covariate=cbind(1,covariates)
    }else{
      intercept_covariate=matrix(1,nrow=nrow(G),ncol=1)
    }
  }
  
  
  d=1
  if (outcome_type[d]=="survival"){
    time=outcome[,1,drop=F]
    status=outcome[,2]
    
    
    # reorder the data based on the time (order is from max to min)
    reorder=order(-time)
    time_reorder=time[reorder,,drop=F]
    status_reorder=status[reorder]
    mediator_reorder=mediator[reorder,,drop=F]
    G_reorder=G[reorder,,drop=F]
    intercept_covariate_reorder=intercept_covariate[reorder,,drop=F]
    
  }
  
  # under the alternative #
  if (verbose){
    print("Under the alternative: theta is nonzero")
  }
  # initial values
  d=1
  if (outcome_type[d]=="survival"){
    if (is.null(covariates)){
      fixed_intercept_covariate=NULL
    }else{
      fixed_intercept_covariate=matrix(0,nrow=ncol(covariates),ncol=num_outcome)
    }
  }else{
    if (is.null(covariates)){
      fixed_intercept_covariate=matrix(0,nrow=1,ncol=num_outcome)
    }else{
      fixed_intercept_covariate=matrix(0,nrow=ncol(covariates)+1,ncol=num_outcome)
    }
  }
  
  
  theta_hat=rep(0,num_outcome)
  gamma_tilde_mat=matrix(0,nrow=ncol(G),ncol=num_outcome)
  phi_vec=rep(1,num_outcome)
  sigma2_gamma_vec=rep(0.001,num_outcome) 
  
  
  
  difference=1
  iteration=0
  check_difference=separate_iteration=rep(0,5)
  names(check_difference)=names(separate_iteration)=c( "fixed_intercept_covariate", "theta_hat", "gamma_tilde_mat", "phi_vec","sigma2_gamma_vec")   
  
  while( (difference>thres)&(iteration<100) ){
    
    iteration=iteration+1
    if (verbose){
      print(iteration)
    }
    old=list("theta_hat"=theta_hat,"gamma_tilde_mat"=gamma_tilde_mat,"sigma2_gamma_vec"=sigma2_gamma_vec,
             "phi_vec"=phi_vec,"fixed_intercept_covariate"=fixed_intercept_covariate)
    
    
    # update fixed_intercept_covariate (intercept and coefficient of covariates) #
    update_fixed_intercept_covariate=function(j){
      W_U_list=list()
      W_first=list()
      W_second=list()
      gradient_list=list()
      hess_list=list()
      dif=1
      iter=0
      step=1
      
      while(dif>thres){
        before=fixed_intercept_covariate[j]
        d=1
        iter=iter+1
        if (iter>100){
          step=step*0.5
        }
        if (outcome_type[d]=="continuous"){
          deri_eta_wrt_zeta=intercept_covariate[,j]
          
          fixedeffectsdat=cbind(intercept_covariate,mediator)
          fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
          gamma_tilde=gamma_tilde_mat[,d]
          phi=phi_vec[d]
          sigma2_gamma=sigma2_gamma_vec[d]
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=eta
          w=as.matrix(rep(1/phi,nrow(G)))
          W_inv=diag(1/w[,1])
          W=diag(w[,1])
          W_U_list[[d]]=W
          
          gradient_list[[d]]=as.numeric( t(deri_eta_wrt_zeta)%*%(outcome-mu))/phi  # L2 approximation
          #gradient=grad(LL_theta,fixedeffects[4])
          hess_list[[d]]=as.numeric(-t(deri_eta_wrt_zeta)%*%W%*%deri_eta_wrt_zeta)
          #hess=hessian(LL_theta,fixedeffects[4])
          
        }
        
        if (outcome_type[d]=="binary"){
          deri_eta_wrt_zeta=intercept_covariate[,j]
          
          fixedeffectsdat=cbind(intercept_covariate,mediator)
          fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
          gamma_tilde=gamma_tilde_mat[,d]
          phi=phi_vec[d]
          sigma2_gamma=sigma2_gamma_vec[d]
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=1/(1+exp(-eta))
          w=mu*(1-mu)
          W_inv=diag(1/w[,1])
          W=diag(w[,1])
          W_U_list[[d]]=W
          if (approxi){
            gradient_list[[d]]=as.numeric( t(deri_eta_wrt_zeta)%*%(outcome-mu))
            hess_list[[d]]=as.numeric(-t(deri_eta_wrt_zeta)%*%W%*%deri_eta_wrt_zeta)
          }else{
            A_1=sigma2_gamma^(-1)*diag(ncol(G))+eigenMapMatMult(eigenMapMatMult(t(G),W),G)
            #A_1_inv=MASS::ginv(A_1)
            A_1_inv=chol2inv(chol(A_1))
            W_first[[d]]=diag( as.numeric(w*(1-2*mu)*deri_eta_wrt_zeta) )
            W_second[[d]]=diag( as.numeric(w*(1-6*w)*deri_eta_wrt_zeta^2) )
            deri_logdet_A1=eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_first[[d]]),G) )
            gradient_list[[d]]=as.numeric( t(deri_eta_wrt_zeta)%*%(outcome-mu))-0.5*sum(diag( deri_logdet_A1 ))  
            hess_list[[d]]=as.numeric(-t(deri_eta_wrt_zeta)%*%W%*%deri_eta_wrt_zeta)-0.5*( sum(diag( -deri_logdet_A1^2+eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_second[[d]]),G) ) ))  ) 
          }
          
        }
        
        if (outcome_type[d]=="count"){
          deri_eta_wrt_zeta=intercept_covariate[,j]
          
          fixedeffectsdat=cbind(intercept_covariate,mediator)
          fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
          gamma_tilde=gamma_tilde_mat[,d]
          phi=phi_vec[d]
          sigma2_gamma=sigma2_gamma_vec[d]
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=exp(eta)
          xi=mu/(phi+mu)
          W=diag( as.numeric((outcome+phi)*xi*(1-xi)) )
          W_U_list[[d]]=W
          if (approxi){
            gradient_list[[d]]=as.numeric(t(deri_eta_wrt_zeta)%*%outcome-t(deri_eta_wrt_zeta*xi)%*%(outcome+phi))
            hess_list[[d]]=as.numeric(-t(deri_eta_wrt_zeta)%*%W%*%deri_eta_wrt_zeta)
          }else{
            A_1=sigma2_gamma^(-1)*diag(ncol(G))+eigenMapMatMult(eigenMapMatMult(t(G),W),G)
            #A_1_inv=MASS::ginv(A_1)
            A_1_inv=chol2inv(chol(A_1))
            W_first[[d]]=diag(as.numeric((outcome+phi)*phi*mu*(phi-mu)/((phi+mu)^3)*deri_eta_wrt_zeta))
            W_second[[d]]=diag(as.numeric((outcome+phi)*phi*mu*(phi^2-4*phi*mu+mu^2)/((phi+mu)^4)*deri_eta_wrt_zeta^2))
            deri_logdet_A1=eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_first[[d]]),G) )
            gradient_list[[d]]=as.numeric(t(deri_eta_wrt_zeta)%*%outcome-t(deri_eta_wrt_zeta*xi)%*%(outcome+phi))-0.5*sum(diag( deri_logdet_A1 ))  
            hess_list[[d]]=as.numeric(-t(deri_eta_wrt_zeta)%*%W%*%deri_eta_wrt_zeta)-0.5*( sum(diag( -deri_logdet_A1^2+eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_second[[d]]),G) ) ))  ) 
          }
          
        }
        
        if (outcome_type[d]=="survival"){
          deri_eta_wrt_zeta=intercept_covariate_reorder[,j]
          
          fixedeffectsdat=cbind(intercept_covariate_reorder,mediator_reorder)
          fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
          gamma_tilde=gamma_tilde_mat[,d]
          sigma2_gamma=sigma2_gamma_vec[d]
          
          eta=fixedeffectsdat%*%fixedeffects+G_reorder%*%gamma_tilde
          exp_eta=exp(eta)
          
          
          if (approxi) {
            weighted_mean=sapply(2:nrow(G_reorder),  function(i) deri_eta_wrt_zeta[i]-weighted.mean(deri_eta_wrt_zeta[1:i] , exp_eta[1:i]) )
            weighted_variance=suppressWarnings(sapply(2:nrow(G_reorder),  function(i) weighted.var( deri_eta_wrt_zeta[1:i] , exp_eta[1:i]) ))
            gradient_list[[d]]=sum(weighted_mean*status_reorder[-1])
            hess_list[[d]]=-sum(weighted_variance*status_reorder[-1])
          }else{
            print("Sorry, we don't have the exact version for survival data; please set approxi=TRUE to use the approximated version.")
          }
          
        }
        
        
        d=1
        if (outcome_type[d]=="continuous"){
          fixed_intercept_covariate[j]=fixed_intercept_covariate[j]+step*as.numeric(gradient_list[[d]]/abs(hess_list[[d]]) )
        }
        
        if (outcome_type[d]=="binary"){
          fixed_intercept_covariate[j]=fixed_intercept_covariate[j]+step*as.numeric(gradient_list[[d]]/abs(hess_list[[d]]) )
        }
        
        if (outcome_type[d]=="count"){
          fixed_intercept_covariate[j]=fixed_intercept_covariate[j]+step*as.numeric(gradient_list[[d]]/abs(hess_list[[d]]) )
        }
        
        if (outcome_type[d]=="survival"){
          fixed_intercept_covariate[j]=fixed_intercept_covariate[j]+step*as.numeric(gradient_list[[d]]/abs(hess_list[[d]]) )
        }
        
        after=fixed_intercept_covariate[j]
        dif=sum((before-after)^2)
        
        
        
      }
      return(fixed_intercept_covariate[j])
    }
    
    
    
    param="fixed_intercept_covariate"
    if (verbose){
      print(param)
    }
    
    if ( !is.null(intercept_covariate) ){
      if (check_difference[param]<num_consecutive){
        
        for (j in 1:ncol(intercept_covariate) ){
          fixed_intercept_covariate[j]=update_fixed_intercept_covariate(j)
        }
        
      }else{
        separate_iteration[param]=separate_iteration[param]+1
        if (separate_iteration[param] %% num_gap == 0){
          for (j in 1:ncol(intercept_covariate) ){
            fixed_intercept_covariate[j]=update_fixed_intercept_covariate(j)
          }
        }
      }
    }
    
    
    if(verbose){
      cat("approximated_log_likelihood =",compute_log_likelihood(theta_hat,gamma_tilde_mat,phi_vec,fixed_intercept_covariate,sigma2_gamma_vec),"\n")     
    }
    
    
    
    # update gamma_tilde_mat
    update_gamma_tilde_mat=function(){
      
      d=1
      iter=0
      step=1
      
      if (outcome_type[d]=="continuous"){
        dif=1
        
        fixedeffectsdat=cbind(intercept_covariate,mediator)
        fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
        gamma_tilde=gamma_tilde_mat[,d]
        phi=phi_vec[d]
        sigma2_gamma=sigma2_gamma_vec[d]
        
        while(dif>thres){
          iter=iter+1
          if (iter>100){
            step=step*0.5
          }
          before=gamma_tilde
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=eta
          w=as.matrix(rep(1/phi,nrow(G)))
          W_inv=diag(1/w[,1])
          W=diag(w[,1])
          Q_inv=eigenMapMatMult(eigenMapMatMult(t(G),W),G)+sigma2_gamma^(-1)*diag(ncol(G))
          #Q=MASS::ginv(Q_inv)
          Q=chol2inv(chol(Q_inv))
          S_prime=t(G)%*%(outcome-mu)/phi-sigma2_gamma^(-1)*gamma_tilde
          gamma_tilde=as.numeric(gamma_tilde+step*Q%*%S_prime)
          
          after=gamma_tilde
          dif=sum((before-after)^2)
        }
        gamma_tilde_mat[,d]=gamma_tilde
      }
      
      if (outcome_type[d]=="binary"){
        dif=1
        
        fixedeffectsdat=cbind(intercept_covariate,mediator)
        fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
        gamma_tilde=gamma_tilde_mat[,d]
        phi=phi_vec[d]
        sigma2_gamma=sigma2_gamma_vec[d]
        while(dif>thres){
          iter=iter+1
          if (iter>100){
            step=step*0.5
          }
          before=gamma_tilde
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=1/(1+exp(-eta))
          w=mu*(1-mu)
          W_inv=diag(1/w[,1])
          W=diag(w[,1])
          Q_inv=eigenMapMatMult(eigenMapMatMult(t(G),W),G)+sigma2_gamma^(-1)*diag(ncol(G))
          #Q=MASS::ginv(Q_inv)
          Q=chol2inv(chol(Q_inv))
          S_prime=t(G)%*%(outcome-mu)-sigma2_gamma^(-1)*gamma_tilde
          gamma_tilde=as.numeric(gamma_tilde+step*Q%*%S_prime)
          after=gamma_tilde
          dif=sum((before-after)^2)
        }
        gamma_tilde_mat[,d]=gamma_tilde
      }
      
      if (outcome_type[d]=="count"){
        dif=1
        
        fixedeffectsdat=cbind(intercept_covariate,mediator)
        fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
        gamma_tilde=gamma_tilde_mat[,d]
        phi=phi_vec[d]
        sigma2_gamma=sigma2_gamma_vec[d]
        while(dif>thres){
          iter=iter+1
          if (iter>100){
            step=step*0.5
          }
          before=gamma_tilde
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=exp(eta)
          xi=mu/(phi+mu)
          u_matrix=diag( as.numeric((outcome+phi)*xi*(1-xi)) )
          Q_inv=eigenMapMatMult(eigenMapMatMult(t(G),u_matrix),G)+sigma2_gamma^(-1)*diag(ncol(G))
          #Q=MASS::ginv(Q_inv)
          Q=chol2inv(chol(Q_inv))
          S_prime=t(G)%*%(outcome- (outcome+phi)*xi)-sigma2_gamma^(-1)*gamma_tilde
          gamma_tilde=as.numeric(gamma_tilde+step*Q%*%S_prime)
          after=gamma_tilde
          dif=sum((before-after)^2)
        }
        gamma_tilde_mat[,d]=gamma_tilde
      }
      
      if (outcome_type[d]=="survival"){
        dif=1
        
        fixedeffectsdat=cbind(intercept_covariate_reorder,mediator_reorder)
        fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
        gamma_tilde=gamma_tilde_mat[,d]
        sigma2_gamma=sigma2_gamma_vec[d]
        
        while(dif>thres){
          iter=iter+1
          if (iter>100){
            step=step*0.5
          }
          before=gamma_tilde
          eta=fixedeffectsdat%*%fixedeffects+G_reorder%*%gamma_tilde
          exp_eta=exp(eta)
          
          if(F){
            gra=function(j){
              deri_eta_wrt_zeta=G_reorder[,j]
              weighted_mean=sapply(2:nrow(G_reorder),  function(i) deri_eta_wrt_zeta[i]-weighted.mean(deri_eta_wrt_zeta[1:i] , exp_eta[1:i]) )
              sum(weighted_mean*status_reorder[-1])
            }
            a=sapply(1:ncol(G_reorder), gra)
          }
          
          
          weighted_mean=sapply(1:ncol(G_reorder),function(j) 
            {
              sum(sapply(2:nrow(G_reorder),  function(i) 
                status_reorder[i]*(G_reorder[i,j]-weighted.mean(G_reorder[1:i,j] , exp_eta[1:i]))  ))
            }
          )
          
          mat=lapply(2:nrow(G_reorder), function(i) status_reorder[i]*cov.wt(G_reorder[1:i,],exp_eta[1:i])$cov )
          weighted_covariance=Reduce(`+`,mat)
          
          gradient=weighted_mean -sigma2_gamma^(-1)*gamma_tilde

          hess=-weighted_covariance-sigma2_gamma^(-1)*diag(ncol(G_reorder))
          gamma_tilde=as.numeric(gamma_tilde- step*MASS::ginv(hess)%*%gradient )
         
          after=gamma_tilde
          dif=sum((before-after)^2)
          
          
        }
        gamma_tilde_mat[,d]=gamma_tilde
      }
      
      return(gamma_tilde_mat)
    }
    
    param="gamma_tilde_mat"
    if (verbose){
      print(param)
    }
    if (check_difference[param]<num_consecutive){
      gamma_tilde_mat=update_gamma_tilde_mat()
    }else{
      separate_iteration[param]=separate_iteration[param]+1
      if (separate_iteration[param] %% num_gap == 0){
        gamma_tilde_mat=update_gamma_tilde_mat()
      }
    }
    
    if(verbose){
      cat("approximated_log_likelihood =",compute_log_likelihood(theta_hat,gamma_tilde_mat,phi_vec,fixed_intercept_covariate,sigma2_gamma_vec),"\n")  
      }
    
    # update phi_vec
    update_phi_vec=function(){
      W_U_list=list()
      W_first=list()
      W_second=list()
      gradient_list=list()
      hess_list=list()
      dif=1
      iter=0
      step=1
      
      while(dif>thres){
        before=phi_vec
        d=1
        iter=iter+1
        if (iter>100){
          step=step*0.5
        }
        
        if (outcome_type[d]=="continuous"){
          
          fixedeffectsdat=cbind(intercept_covariate,mediator)
          fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
          gamma_tilde=gamma_tilde_mat[,d]
          phi=phi_vec[d]
          sigma2_gamma=sigma2_gamma_vec[d]
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=eta
          w=as.matrix(rep(1/phi,nrow(G)))
          W_inv=diag(1/w[,1])
          W=diag(w[,1])
          W_U_list[[d]]=W
          if (approxi){
            gradient_list[[d]]=sum((outcome-mu)^2)/(2*phi^2)-nrow(G)/(2*phi)
            hess_list[[d]]=-sum((outcome-mu)^2)/(phi^3)+nrow(G)/(2*phi^2)
          }else{
            A_1=sigma2_gamma^(-1)*diag(ncol(G))+eigenMapMatMult(eigenMapMatMult(t(G),W),G)
            A_1_inv=chol2inv(chol(A_1))
            W_first[[d]]=diag( rep(-1/(phi^2),nrow(G)) )
            W_second[[d]]=diag( rep(2/(phi^3),nrow(G)) )
            deri_logdet_A1=eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_first[[d]]),G) )
            
            gradient_list[[d]]=sum((outcome-mu)^2)/(2*phi^2)-nrow(G)/(2*phi)-0.5*sum(diag( deri_logdet_A1 ))  
            hess_list[[d]]=-sum((outcome-mu)^2)/(phi^3)+nrow(G)/(2*phi^2) -0.5*( sum(diag( -deri_logdet_A1^2+eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_second[[d]]),G) ) ))  ) 
            
          }
          
        }
        
        if (outcome_type[d]=="binary"){
          # dispersion parameter is equal to 1
        }
        
        if (outcome_type[d]=="count"){
          
          fixedeffectsdat=cbind(intercept_covariate,mediator)
          fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
          gamma_tilde=gamma_tilde_mat[,d]
          phi=phi_vec[d]
          sigma2_gamma=sigma2_gamma_vec[d]
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=exp(eta)
          xi=mu/(phi+mu)
          W=diag( as.numeric((outcome+phi)*xi*(1-xi)) )
          W_U_list[[d]]=W
          if (approxi){
            gradient_list[[d]]=sum( -log(phi+mu)-(outcome+phi)/(phi+mu)+1+log(phi)-(digamma(phi)-digamma(phi+outcome)) )
            hess_list[[d]]=sum( -1/(phi+mu) + (outcome-mu)/((phi+mu)^2) + 1/phi - (psigamma(phi,deriv=1)-psigamma(phi+outcome,deriv=1)) )
          }else{
            A_1=sigma2_gamma^(-1)*diag(ncol(G))+eigenMapMatMult(eigenMapMatMult(t(G),W),G)
            #A_1_inv=MASS::ginv(A_1)
            A_1_inv=chol2inv(chol(A_1))
            W_first[[d]]= diag( as.numeric(xi*(1-xi)+(outcome+phi)*( (mu^2-phi*mu)/((phi+mu)^3) )) )
            W_second[[d]]=diag(as.numeric( 2*(mu^2-phi*mu)/((phi+mu)^3) + (outcome+phi)*(2*phi*mu-4*mu^2)/((phi+mu)^4) )) 
            deri_logdet_A1=eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_first[[d]]),G) )
            gradient_list[[d]]=sum( -log(phi+mu)-(outcome+phi)/(phi+mu)+1+log(phi)-(digamma(phi)-digamma(phi+outcome)) )-0.5*sum(diag( deri_logdet_A1 ))  
            hess_list[[d]]=sum( -1/(phi+mu) + (outcome-mu)/((phi+mu)^2) + 1/phi - (psigamma(phi,deriv=1)-psigamma(phi+outcome,deriv=1)) ) -0.5*( sum(diag( -deri_logdet_A1^2+eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_second[[d]]),G) ) ))  ) 
          }
          
        }
        
        if (outcome_type[d]=="survival"){
          # no dispersion parameter for survival data
        }
        
        
        
        
        
        d=1
        if (outcome_type[d]=="continuous"){
          phi_vec[d]=phi_vec[d]+step*as.numeric(gradient_list[[d]]/abs(hess_list[[d]]) )
        }
        
        if (outcome_type[d]=="binary"){
          phi_vec[d]=1
        }
        
        if (outcome_type[d]=="count"){
          phi_vec[d]=max(min(phi_vec[d]+step*as.numeric(gradient_list[[d]]/abs(hess_list[[d]]) ),1000),0.001)
        }
        
        if (outcome_type[d]=="survival"){
          # no dispersion parameter
        }
       
        after=phi_vec
        dif=sum((before-after)^2)
      }
      
      
      
      return(phi_vec)
    }
    
    param="phi_vec"
    if (verbose){
      print(param)
    }
    if (check_difference[param]<num_consecutive){
      phi_vec=update_phi_vec()
    }else{
      separate_iteration[param]=separate_iteration[param]+1
      if (separate_iteration[param] %% num_gap == 0){
        phi_vec=update_phi_vec()
      }
    }
    
    if(verbose){
      cat("approximated_log_likelihood =",compute_log_likelihood(theta_hat,gamma_tilde_mat,phi_vec,fixed_intercept_covariate,sigma2_gamma_vec),"\n")  
    }
    
    # update sigma2_gamma_vec
    update_sigma2_gamma_vec=function(){
      d=1
      iter=0
      step=1
      if (outcome_type[d]=="continuous"){
        
        fixedeffectsdat=cbind(intercept_covariate,mediator)
        fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
        gamma_tilde=gamma_tilde_mat[,d]
        phi=phi_vec[d]
        sigma2_gamma=sigma2_gamma_vec[d]
        dif=1
        while(dif>thres){
          iter=iter+1
          if (iter>100){
            step=step*0.5
          }
          before=sigma2_gamma
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=eta
          w=as.matrix(rep(1/phi,nrow(G)))
          W_inv=diag(1/w[,1])
          W=diag(w[,1])
          #tGWG=eigenMapMatMult(eigenMapMatMult(t(G),W),G)
          Iq=diag(ncol(G))
          Q_inv=eigenMapMatMult(eigenMapMatMult(t(G),W),G)+sigma2_gamma^(-1)*Iq
          Q=MASS::ginv(Q_inv)
          #Q=chol2inv(chol(Q_inv))
          
          A=sigma2_gamma^(-1)*Iq-sigma2_gamma^(-2)*Q
          A2=eigenMapMatMult(A,A)
          #gradient=0.5*( -sum(diag( MASS::ginv(Iq+sigma2_gamma*tGWG)%*%tGWG ))+t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-2) )
          gradient=0.5*( -sum(diag( A )) + t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-2) )  
          #grad(LL_sigma2_gamma,sigma2_gamma)
          #hess=0.5*(sum(diag( (MASS::ginv(Iq+sigma2_gamma*tGWG)%*%tGWG)^2 )) -2*t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-3) )
          hess=0.5*( sum(diag(A2))  -2*t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-3) )
          #hessian(LL_sigma2_gamma,sigma2_gamma)
          
          sigma2_gamma=min(max(sigma2_gamma + step*as.numeric(gradient/ abs(hess) ),0.001),10)
          
          after=sigma2_gamma
          dif=((before-after)^2)
        }
        sigma2_gamma_vec[d]=sigma2_gamma
      }
      
      
      if (outcome_type[d]=="binary"){
        
        fixedeffectsdat=cbind(intercept_covariate,mediator)
        fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
        gamma_tilde=gamma_tilde_mat[,d]
        phi=phi_vec[d]
        sigma2_gamma=sigma2_gamma_vec[d]
        dif=1
        while(dif>thres){
          iter=iter+1
          if (iter>100){
            step=step*0.5
          }
          before=sigma2_gamma
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=1/(1+exp(-eta))
          w=mu*(1-mu)
          W_inv=diag(1/w[,1])
          W=diag(w[,1])
          #tGWG=eigenMapMatMult(eigenMapMatMult(t(G),W),G)
          Iq=diag(ncol(G))
          Q_inv=eigenMapMatMult(eigenMapMatMult(t(G),W),G)+sigma2_gamma^(-1)*Iq
          #Q=MASS::ginv(Q_inv)
          Q=chol2inv(chol(Q_inv))
          
          A=sigma2_gamma^(-1)*Iq-sigma2_gamma^(-2)*Q
          A2=eigenMapMatMult(A,A)
          #gradient=0.5*( -sum(diag( MASS::ginv(Iq+sigma2_gamma*tGWG)%*%tGWG ))+t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-2) )
          gradient=0.5*( -sum(diag( A )) + t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-2) )  
          #grad(LL_sigma2_gamma,sigma2_gamma)
          #hess=0.5*(sum(diag( (MASS::ginv(Iq+sigma2_gamma*tGWG)%*%tGWG)^2 )) -2*t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-3) )
          hess=0.5*( sum(diag(A2))  -2*t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-3) )
          #hessian(LL_sigma2_gamma,sigma2_gamma)
          sigma2_gamma=min(max(sigma2_gamma + step*as.numeric(gradient/ abs(hess) ),0.001),10)
          after=sigma2_gamma
          dif=((before-after)^2)
          
        }
        sigma2_gamma_vec[d]=sigma2_gamma
      }
      
      if (outcome_type[d]=="count"){
        
        fixedeffectsdat=cbind(intercept_covariate,mediator)
        fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
        gamma_tilde=gamma_tilde_mat[,d]
        phi=phi_vec[d]
        sigma2_gamma=sigma2_gamma_vec[d]
        dif=1
        while(dif>thres){
          iter=iter+1
          if (iter>100){
            step=step*0.5
          }
          before=sigma2_gamma
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=exp(eta)
          xi=mu/(phi+mu)
          u_matrix=diag( as.numeric((outcome+phi)*xi*(1-xi)) )
          
          Q_inv=eigenMapMatMult(eigenMapMatMult(t(G),u_matrix),G)+sigma2_gamma^(-1)*diag(ncol(G))
          #Q=MASS::ginv(Q_inv)
          Q=chol2inv(chol(Q_inv))
          gradient=0.5*sum(diag(sigma2_gamma^(-2)*(gamma_tilde%*%t(gamma_tilde)+Q-sigma2_gamma*diag(ncol(G)))))
          # gradient=numDeriv::grad(LL_sigma2_gamma,sigma2_gamma)
          hess=-sum(diag(gamma_tilde%*%t(gamma_tilde)))*sigma2_gamma^(-3)-0.5*ncol(G)
          # hess=numDeriv::hessian(LL_sigma2_gamma,sigma2_gamma)
          
          sigma2_gamma=min(max(sigma2_gamma + step*as.numeric(gradient/ abs(hess) ),0.001),10)
          after=sigma2_gamma
          dif=((before-after)^2)
        }
        sigma2_gamma_vec[d]=sigma2_gamma
      }
      
      
      if (outcome_type[d]=="survival"){  
        dif=1
        fixedeffectsdat=cbind(intercept_covariate_reorder,mediator_reorder)
        fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
        gamma_tilde=gamma_tilde_mat[,d]
        sigma2_gamma=sigma2_gamma_vec[d]
        while(dif>thres){
          iter=iter+1
          if (iter>100){
            step=step*0.5
          }
          before=sigma2_gamma
          gradient=0.5*(sigma2_gamma^(-2)*t(gamma_tilde)%*%gamma_tilde-sigma2_gamma^(-1)*ncol(G_reorder) )
          hess=-t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-3)+0.5*ncol(G_reorder)*sigma2_gamma^(-2)
          sigma2_gamma=min(max(sigma2_gamma + step*as.numeric(gradient/ abs(hess) ),0.001),10)
          after=sigma2_gamma
          dif=((before-after)^2)
        }
        sigma2_gamma_vec[d]=sigma2_gamma
      }
      return(sigma2_gamma_vec)
    }
    
    param="sigma2_gamma_vec"
    if (verbose){
      print(param)
    }
    if (check_difference[param]<num_consecutive){
      sigma2_gamma_vec=update_sigma2_gamma_vec()
    }else{
      separate_iteration[param]=separate_iteration[param]+1
      if (separate_iteration[param] %% num_gap == 0){
        sigma2_gamma_vec=update_sigma2_gamma_vec()
      }
    }
    
    if(verbose){
      cat("approximated_log_likelihood =",compute_log_likelihood(theta_hat,gamma_tilde_mat,phi_vec,fixed_intercept_covariate,sigma2_gamma_vec),"\n")  
    }
    
    # update theta_hat 
    update_theta_hat=function(){
      W_U_list=list()
      W_first=list()
      W_second=list()
      gradient_list=list()
      hess_list=list()
      dif=1
      iter=0
      step=1
      
      while(dif>thres){
        iter=iter+1
        if (iter>100){
          step=step*0.5
        }
        
        before=theta_hat
        
        d=1
        if (outcome_type[d]=="continuous"){
          deri_eta_wrt_zeta=mediator
          
          fixedeffectsdat=cbind(intercept_covariate,mediator)
          fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
          gamma_tilde=gamma_tilde_mat[,d]
          phi=phi_vec[d]
          sigma2_gamma=sigma2_gamma_vec[d]
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=eta
          w=as.matrix(rep(1/phi,nrow(G)))
          W_inv=diag(1/w[,1])
          W=diag(w[,1])
          W_U_list[[d]]=W
          
          gradient_list[[d]]=as.numeric( t(deri_eta_wrt_zeta)%*%(outcome-mu))/phi  # L2 approximation
          #gradient=grad(LL_theta,fixedeffects[4])
          hess_list[[d]]=as.numeric(-t(deri_eta_wrt_zeta)%*%W%*%deri_eta_wrt_zeta)
          #hess=hessian(LL_theta,fixedeffects[4])
          
        }
        
        if (outcome_type[d]=="binary"){
          deri_eta_wrt_zeta=mediator
          
          fixedeffectsdat=cbind(intercept_covariate,mediator)
          fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
          gamma_tilde=gamma_tilde_mat[,d]
          phi=phi_vec[d]
          sigma2_gamma=sigma2_gamma_vec[d]
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=1/(1+exp(-eta))
          w=mu*(1-mu)
          W_inv=diag(1/w[,1])
          W=diag(w[,1])
          W_U_list[[d]]=W
          if (approxi){
            gradient_list[[d]]=as.numeric( t(deri_eta_wrt_zeta)%*%(outcome-mu))
            hess_list[[d]]=as.numeric(-t(deri_eta_wrt_zeta)%*%W%*%deri_eta_wrt_zeta)
          }else{
            A_1=sigma2_gamma^(-1)*diag(ncol(G))+eigenMapMatMult(eigenMapMatMult(t(G),W),G)
            #A_1_inv=MASS::ginv(A_1)
            A_1_inv=chol2inv(chol(A_1))
            W_first[[d]]=diag( as.numeric(w*(1-2*mu)*deri_eta_wrt_zeta) )
            W_second[[d]]=diag( as.numeric(w*(1-6*w)*deri_eta_wrt_zeta^2) )
            deri_logdet_A1=eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_first[[d]]),G) )
            gradient_list[[d]]=as.numeric( t(deri_eta_wrt_zeta)%*%(outcome-mu))-0.5*sum(diag( deri_logdet_A1 ))  
            hess_list[[d]]=as.numeric(-t(deri_eta_wrt_zeta)%*%W%*%deri_eta_wrt_zeta)-0.5*( sum(diag( -deri_logdet_A1^2+eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_second[[d]]),G) ) ))  ) 
          }
          
        }
        
        if (outcome_type[d]=="count"){
          deri_eta_wrt_zeta=mediator
          
          fixedeffectsdat=cbind(intercept_covariate,mediator)
          fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
          gamma_tilde=gamma_tilde_mat[,d]
          phi=phi_vec[d]
          sigma2_gamma=sigma2_gamma_vec[d]
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=exp(eta)
          xi=mu/(phi+mu)
          W=diag( as.numeric((outcome+phi)*xi*(1-xi)) )
          W_U_list[[d]]=W
          if (approxi){
            gradient_list[[d]]=as.numeric(t(deri_eta_wrt_zeta)%*%outcome-t(deri_eta_wrt_zeta*xi)%*%(outcome+phi))
            hess_list[[d]]=as.numeric(-t(deri_eta_wrt_zeta)%*%W%*%deri_eta_wrt_zeta)
          }else{
            A_1=sigma2_gamma^(-1)*diag(ncol(G))+eigenMapMatMult(eigenMapMatMult(t(G),W),G)
            #A_1_inv=MASS::ginv(A_1)
            A_1_inv=chol2inv(chol(A_1))
            W_first[[d]]=diag(as.numeric((outcome+phi)*phi*mu*(phi-mu)/((phi+mu)^3)*deri_eta_wrt_zeta))
            W_second[[d]]=diag(as.numeric((outcome+phi)*phi*mu*(phi^2-4*phi*mu+mu^2)/((phi+mu)^4)*deri_eta_wrt_zeta^2))
            deri_logdet_A1=eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_first[[d]]),G) )
            gradient_list[[d]]=as.numeric(t(deri_eta_wrt_zeta)%*%outcome-t(deri_eta_wrt_zeta*xi)%*%(outcome+phi))-0.5*sum(diag( deri_logdet_A1 ))  
            hess_list[[d]]=as.numeric(-t(deri_eta_wrt_zeta)%*%W%*%deri_eta_wrt_zeta)-0.5*( sum(diag( -deri_logdet_A1^2+eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_second[[d]]),G) ) ))  ) 
          }
          
        }
        
        if (outcome_type[d]=="survival"){
          deri_eta_wrt_zeta=mediator_reorder
          
          fixedeffectsdat=cbind(intercept_covariate_reorder,mediator_reorder)
          fixedeffects=c(fixed_intercept_covariate[,d],theta_hat[d])
          gamma_tilde=gamma_tilde_mat[,d]
          sigma2_gamma=sigma2_gamma_vec[d]
          
          eta=fixedeffectsdat%*%fixedeffects+G_reorder%*%gamma_tilde
          exp_eta=exp(eta)
          
          if (approxi) {
            weighted_mean=sapply(2:nrow(G_reorder),  function(i) deri_eta_wrt_zeta[i]-weighted.mean(deri_eta_wrt_zeta[1:i] , exp_eta[1:i]) )
            weighted_variance=suppressWarnings(sapply(2:nrow(G_reorder),  function(i) weighted.var( deri_eta_wrt_zeta[1:i] , exp_eta[1:i]) ))
            gradient_list[[d]]=sum(weighted_mean*status_reorder[-1])
            hess_list[[d]]=-sum(weighted_variance*status_reorder[-1])
          }else{
            print("Sorry, we don't have the exact version for survival data; please set approxi=TRUE to use the approximated version.")
          }
        }
        
        
        
        
        
        d=1
        if (outcome_type[d]=="continuous"){
          theta_hat[d]=theta_hat[d]+step*as.numeric(gradient_list[[d]]/abs(hess_list[[d]]) )
        }
        
        if (outcome_type[d]=="binary"){
          theta_hat[d]=theta_hat[d]+step*as.numeric(gradient_list[[d]]/abs(hess_list[[d]]) )
        }
        
        if (outcome_type[d]=="count"){
          theta_hat[d]=theta_hat[d]+step*as.numeric(gradient_list[[d]]/abs(hess_list[[d]]) )
        }
        
        if (outcome_type[d]=="survival"){
          theta_hat[d]=theta_hat[d]+step*as.numeric(gradient_list[[d]]/abs(hess_list[[d]]) )
        }
        
        
        
        after=theta_hat
        dif=sum((before-after)^2)
        
      }
      return(theta_hat)
    }
    
    param="theta_hat"
    if (verbose){
      print(param)
    }
    if (check_difference[param]<num_consecutive){
      theta_hat=update_theta_hat()
    }else{
      separate_iteration[param]=separate_iteration[param]+1
      if (separate_iteration[param] %% num_gap == 0){
        theta_hat=update_theta_hat()
      }
    }
    
    if(verbose){
      cat("approximated_log_likelihood =",compute_log_likelihood(theta_hat,gamma_tilde_mat,phi_vec,fixed_intercept_covariate,sigma2_gamma_vec),"\n")   
    }
    
    new=list("theta_hat"=theta_hat,"gamma_tilde_mat"=gamma_tilde_mat,"sigma2_gamma_vec"=sigma2_gamma_vec,
             "phi_vec"=phi_vec,"fixed_intercept_covariate"=fixed_intercept_covariate)
    
    for (param in c( "fixed_intercept_covariate", "theta_hat", "gamma_tilde_mat", "phi_vec", "sigma2_gamma_vec") ){
      if (verbose){
        cat(param,"_diff=",sum((old[[param]]-new[[param]])^2),"\n")
      }
      if ( sum((old[[param]]-new[[param]])^2) < thres ){
        check_difference[param]=check_difference[param]+1
      }
    }
    
    
    
    difference=sum((as.numeric( unlist(new) )-as.numeric(unlist(old)) )^2)
    if(verbose){
      cat("check_difference=",check_difference,"\n")
      cat("separate_iteration=",separate_iteration,"\n")
      cat("fixed_intercept_covariate=",fixed_intercept_covariate,"\n")
      cat("phi_vec=",phi_vec,"\n")
      cat("theta_hat=",theta_hat,"\n")
      cat( "difference=",difference ,"\n")
    }
    
  }
  
  LL_alt=compute_log_likelihood(theta_hat,gamma_tilde_mat,phi_vec,fixed_intercept_covariate,sigma2_gamma_vec)
  
  
  # under the null #
  if (verbose){
    print("Under the null: theta is zero")
  }
  # initial values
  
  d=1
  if (outcome_type[d]=="survival"){
    if (is.null(covariates)){
      fixed_intercept_covariate=NULL
    }else{
      fixed_intercept_covariate=matrix(0,nrow=ncol(covariates),ncol=num_outcome)
    }
  }else{
    if (is.null(covariates)){
      fixed_intercept_covariate=matrix(0,nrow=1,ncol=num_outcome)
    }else{
      fixed_intercept_covariate=matrix(0,nrow=ncol(covariates)+1,ncol=num_outcome)
    }
  }
  
  gamma_tilde_mat=matrix(0,nrow=ncol(G),ncol=num_outcome)
  phi_vec=rep(1,num_outcome)
  sigma2_gamma_vec=rep(0.001,num_outcome) 
  
  
  
  difference=1
  iteration=0
  check_difference=separate_iteration=rep(0,5)
  names(check_difference)=names(separate_iteration)=c( "fixed_intercept_covariate", "theta_hat", "gamma_tilde_mat", "phi_vec","sigma2_gamma_vec")   
  
  
  while( (difference>thres)&(iteration<100) ){
    
    iteration=iteration+1
    if (verbose){
      print(iteration)
    }
    
    old=list("gamma_tilde_mat"=gamma_tilde_mat,"sigma2_gamma_vec"=sigma2_gamma_vec,
             "phi_vec"=phi_vec,"fixed_intercept_covariate"=fixed_intercept_covariate)
    
    
    # update fixed_intercept_covariate (intercept and coefficient of covariates) #
    update_fixed_intercept_covariate=function(j){
      W_U_list=list()
      W_first=list()
      W_second=list()
      gradient_list=list()
      hess_list=list()
      dif=1
      iter=0
      step=1
      while(dif>thres){
        before=fixed_intercept_covariate[j]
        iter=iter+1
        if (iter>100){
          step=step*0.5
        }
        d=1
        if (outcome_type[d]=="continuous"){
          deri_eta_wrt_zeta=intercept_covariate[,j]
          
          fixedeffectsdat=cbind(intercept_covariate)
          fixedeffects=c(fixed_intercept_covariate[,d])
          gamma_tilde=gamma_tilde_mat[,d]
          phi=phi_vec[d]
          sigma2_gamma=sigma2_gamma_vec[d]
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=eta
          w=as.matrix(rep(1/phi,nrow(G)))
          W_inv=diag(1/w[,1])
          W=diag(w[,1])
          W_U_list[[d]]=W
          
          gradient_list[[d]]=as.numeric( t(deri_eta_wrt_zeta)%*%(outcome-mu))/phi  # L2 approximation
          #gradient=grad(LL_theta,fixedeffects[4])
          hess_list[[d]]=as.numeric(-t(deri_eta_wrt_zeta)%*%W%*%deri_eta_wrt_zeta)
          #hess=hessian(LL_theta,fixedeffects[4])
          
        }
        
        if (outcome_type[d]=="binary"){
          deri_eta_wrt_zeta=intercept_covariate[,j]
          
          fixedeffectsdat=cbind(intercept_covariate)
          fixedeffects=c(fixed_intercept_covariate[,d])
          gamma_tilde=gamma_tilde_mat[,d]
          phi=phi_vec[d]
          sigma2_gamma=sigma2_gamma_vec[d]
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=1/(1+exp(-eta))
          w=mu*(1-mu)
          W_inv=diag(1/w[,1])
          W=diag(w[,1])
          W_U_list[[d]]=W
          if (approxi){
            gradient_list[[d]]=as.numeric( t(deri_eta_wrt_zeta)%*%(outcome-mu))
            hess_list[[d]]=as.numeric(-t(deri_eta_wrt_zeta)%*%W%*%deri_eta_wrt_zeta)
          }else{
            A_1=sigma2_gamma^(-1)*diag(ncol(G))+eigenMapMatMult(eigenMapMatMult(t(G),W),G)
            #A_1_inv=MASS::ginv(A_1)
            A_1_inv=chol2inv(chol(A_1))
            W_first[[d]]=diag( as.numeric(w*(1-2*mu)*deri_eta_wrt_zeta) )
            W_second[[d]]=diag( as.numeric(w*(1-6*w)*deri_eta_wrt_zeta^2) )
            deri_logdet_A1=eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_first[[d]]),G) )
            gradient_list[[d]]=as.numeric( t(deri_eta_wrt_zeta)%*%(outcome-mu))-0.5*sum(diag( deri_logdet_A1 ))  
            hess_list[[d]]=as.numeric(-t(deri_eta_wrt_zeta)%*%W%*%deri_eta_wrt_zeta)-0.5*( sum(diag( -deri_logdet_A1^2+eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_second[[d]]),G) ) ))  ) 
          }
          
        }
        
        if (outcome_type[d]=="count"){
          deri_eta_wrt_zeta=intercept_covariate[,j]
          
          fixedeffectsdat=cbind(intercept_covariate)
          fixedeffects=c(fixed_intercept_covariate[,d])
          gamma_tilde=gamma_tilde_mat[,d]
          phi=phi_vec[d]
          sigma2_gamma=sigma2_gamma_vec[d]
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=exp(eta)
          xi=mu/(phi+mu)
          W=diag( as.numeric((outcome+phi)*xi*(1-xi)) )
          W_U_list[[d]]=W
          if (approxi){
            gradient_list[[d]]=as.numeric(t(deri_eta_wrt_zeta)%*%outcome-t(deri_eta_wrt_zeta*xi)%*%(outcome+phi))
            hess_list[[d]]=as.numeric(-t(deri_eta_wrt_zeta)%*%W%*%deri_eta_wrt_zeta)
          }else{
            A_1=sigma2_gamma^(-1)*diag(ncol(G))+eigenMapMatMult(eigenMapMatMult(t(G),W),G)
            #A_1_inv=MASS::ginv(A_1)
            A_1_inv=chol2inv(chol(A_1))
            W_first[[d]]=diag(as.numeric((outcome+phi)*phi*mu*(phi-mu)/((phi+mu)^3)*deri_eta_wrt_zeta))
            W_second[[d]]=diag(as.numeric((outcome+phi)*phi*mu*(phi^2-4*phi*mu+mu^2)/((phi+mu)^4)*deri_eta_wrt_zeta^2))
            deri_logdet_A1=eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_first[[d]]),G) )
            gradient_list[[d]]=as.numeric(t(deri_eta_wrt_zeta)%*%outcome-t(deri_eta_wrt_zeta*xi)%*%(outcome+phi))-0.5*sum(diag( deri_logdet_A1 ))  
            hess_list[[d]]=as.numeric(-t(deri_eta_wrt_zeta)%*%W%*%deri_eta_wrt_zeta)-0.5*( sum(diag( -deri_logdet_A1^2+eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_second[[d]]),G) ) ))  ) 
          }
          
        }
        
        if (outcome_type[d]=="survival"){
          deri_eta_wrt_zeta=intercept_covariate_reorder[,j]
          
          fixedeffectsdat=cbind(intercept_covariate_reorder)
          fixedeffects=c(fixed_intercept_covariate[,d])
          gamma_tilde=gamma_tilde_mat[,d]
          sigma2_gamma=sigma2_gamma_vec[d]
          
          eta=fixedeffectsdat%*%fixedeffects+G_reorder%*%gamma_tilde
          exp_eta=exp(eta)
          
          if (approxi) {
            weighted_mean=sapply(2:nrow(G_reorder),  function(i) deri_eta_wrt_zeta[i]-weighted.mean(deri_eta_wrt_zeta[1:i] , exp_eta[1:i]) )
            weighted_variance=suppressWarnings(sapply(2:nrow(G_reorder),  function(i) weighted.var( deri_eta_wrt_zeta[1:i] , exp_eta[1:i]) ))
            gradient_list[[d]]=sum(weighted_mean*status_reorder[-1])
            hess_list[[d]]=-sum(weighted_variance*status_reorder[-1])
          }else{
            print("Sorry, we don't have the exact version for survival data; please set approxi=TRUE to use the approximated version.")
          }
        }
        
        
        
        
        
        d=1
        if (outcome_type[d]=="continuous"){
          fixed_intercept_covariate[j]=fixed_intercept_covariate[j]+step*as.numeric(gradient_list[[d]]/abs(hess_list[[d]]) )
        }
        
        if (outcome_type[d]=="binary"){
          fixed_intercept_covariate[j]=fixed_intercept_covariate[j]+step*as.numeric(gradient_list[[d]]/abs(hess_list[[d]]) )
        }
        
        if (outcome_type[d]=="count"){
          fixed_intercept_covariate[j]=fixed_intercept_covariate[j]+step*as.numeric(gradient_list[[d]]/abs(hess_list[[d]]) )
        }
        
        if (outcome_type[d]=="survival"){
          fixed_intercept_covariate[j]=fixed_intercept_covariate[j]+step*as.numeric(gradient_list[[d]]/abs(hess_list[[d]]) )
        }
        
        
        
        after=fixed_intercept_covariate[j]
        dif=sum((before-after)^2)
        
       
        
      }
      return(fixed_intercept_covariate[j])
    }
    
    
    
    param="fixed_intercept_covariate"
    if (verbose){
      print(param)
    }
    if ( !is.null(intercept_covariate) ){
      if (check_difference[param]<num_consecutive){
        
        for (j in 1:ncol(intercept_covariate) ){
          fixed_intercept_covariate[j]=update_fixed_intercept_covariate(j)
        }
        
      }else{
        separate_iteration[param]=separate_iteration[param]+1
        if (separate_iteration[param] %% num_gap == 0){
          for (j in 1:ncol(intercept_covariate) ){
            fixed_intercept_covariate[j]=update_fixed_intercept_covariate(j)
          }
        }
      }
    }
    
    
    if(verbose){
      cat("approximated_log_likelihood =",compute_log_likelihood(theta_hat=NULL,gamma_tilde_mat,phi_vec,fixed_intercept_covariate,sigma2_gamma_vec),"\n")   
    }
    
    
    
    # update gamma_tilde_mat
    update_gamma_tilde_mat=function(){
      iter=0
      step=1
      d=1
      if (outcome_type[d]=="continuous"){
        dif=1
        
        fixedeffectsdat=cbind(intercept_covariate)
        fixedeffects=c(fixed_intercept_covariate[,d])
        gamma_tilde=gamma_tilde_mat[,d]
        phi=phi_vec[d]
        sigma2_gamma=sigma2_gamma_vec[d]
        
        while(dif>thres){
          iter=iter+1
          if (iter>100){
            step=step*0.5
          }
          before=gamma_tilde
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=eta
          w=as.matrix(rep(1/phi,nrow(G)))
          W_inv=diag(1/w[,1])
          W=diag(w[,1])
          Q_inv=eigenMapMatMult(eigenMapMatMult(t(G),W),G)+sigma2_gamma^(-1)*diag(ncol(G))
          #Q=MASS::ginv(Q_inv)
          Q=chol2inv(chol(Q_inv))
          S_prime=t(G)%*%(outcome-mu)/phi-sigma2_gamma^(-1)*gamma_tilde
          gamma_tilde=as.numeric(gamma_tilde+step*Q%*%S_prime)
          
          after=gamma_tilde
          dif=sum((before-after)^2)
        }
        gamma_tilde_mat[,d]=gamma_tilde
      }
      
      if (outcome_type[d]=="binary"){
        dif=1
        
        fixedeffectsdat=cbind(intercept_covariate)
        fixedeffects=c(fixed_intercept_covariate[,d])
        gamma_tilde=gamma_tilde_mat[,d]
        phi=phi_vec[d]
        sigma2_gamma=sigma2_gamma_vec[d]
        while(dif>thres){
          iter=iter+1
          if (iter>100){
            step=step*0.5
          }
          before=gamma_tilde
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=1/(1+exp(-eta))
          w=mu*(1-mu)
          W_inv=diag(1/w[,1])
          W=diag(w[,1])
          Q_inv=eigenMapMatMult(eigenMapMatMult(t(G),W),G)+sigma2_gamma^(-1)*diag(ncol(G))
          #Q=MASS::ginv(Q_inv)
          Q=chol2inv(chol(Q_inv))
          S_prime=t(G)%*%(outcome-mu)-sigma2_gamma^(-1)*gamma_tilde
          gamma_tilde=as.numeric(gamma_tilde+step*Q%*%S_prime)
          after=gamma_tilde
          dif=sum((before-after)^2)
        }
        gamma_tilde_mat[,d]=gamma_tilde
      }
      
      if (outcome_type[d]=="count"){
        dif=1
        
        fixedeffectsdat=cbind(intercept_covariate)
        fixedeffects=c(fixed_intercept_covariate[,d])
        gamma_tilde=gamma_tilde_mat[,d]
        phi=phi_vec[d]
        sigma2_gamma=sigma2_gamma_vec[d]
        while(dif>thres){
          iter=iter+1
          if (iter>100){
            step=step*0.5
          }
          before=gamma_tilde
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=exp(eta)
          xi=mu/(phi+mu)
          u_matrix=diag( as.numeric((outcome+phi)*xi*(1-xi)) )
          Q_inv=eigenMapMatMult(eigenMapMatMult(t(G),u_matrix),G)+sigma2_gamma^(-1)*diag(ncol(G))
          #Q=MASS::ginv(Q_inv)
          Q=chol2inv(chol(Q_inv))
          S_prime=t(G)%*%(outcome- (outcome+phi)*xi)-sigma2_gamma^(-1)*gamma_tilde
          gamma_tilde=as.numeric(gamma_tilde+step*Q%*%S_prime)
          after=gamma_tilde
          dif=sum((before-after)^2)
        }
        gamma_tilde_mat[,d]=gamma_tilde
      }
      
      
      if (outcome_type[d]=="survival"){
        dif=1
        
        fixedeffectsdat=cbind(intercept_covariate_reorder)
        fixedeffects=c(fixed_intercept_covariate[,d])
        if (is.null(fixedeffectsdat)){
          fixedeffectsdat=matrix(0,nrow=nrow(G)); fixedeffects=0
        }
        gamma_tilde=gamma_tilde_mat[,d]
        sigma2_gamma=sigma2_gamma_vec[d]
        
        while(dif>thres){
          iter=iter+1
          if (iter>100){
            step=step*0.5
          }
          before=gamma_tilde
          eta=fixedeffectsdat%*%fixedeffects+G_reorder%*%gamma_tilde
          exp_eta=exp(eta)
          
          
          weighted_mean=sapply(1:ncol(G_reorder),function(j) 
          {
            sum(sapply(2:nrow(G_reorder),  function(i) 
              status_reorder[i]*(G_reorder[i,j]-weighted.mean(G_reorder[1:i,j] , exp_eta[1:i]))  ))
          }
          )
          
          mat=lapply(2:nrow(G_reorder), function(i) status_reorder[i]*cov.wt(G_reorder[1:i,],exp_eta[1:i])$cov )
          weighted_covariance=Reduce(`+`,mat)
          
          gradient=weighted_mean -sigma2_gamma^(-1)*gamma_tilde
          
          hess=-weighted_covariance-sigma2_gamma^(-1)*diag(ncol(G_reorder))
          gamma_tilde=as.numeric(gamma_tilde- step*MASS::ginv(hess)%*%gradient )
          
          after=gamma_tilde
          dif=sum((before-after)^2)
        }
        gamma_tilde_mat[,d]=gamma_tilde
      }
      
      return(gamma_tilde_mat)
    }
    
    param="gamma_tilde_mat"
    if (verbose){
      print(param)
    }
    if (check_difference[param]<num_consecutive){
      gamma_tilde_mat=update_gamma_tilde_mat()
    }else{
      separate_iteration[param]=separate_iteration[param]+1
      if (separate_iteration[param] %% num_gap == 0){
        gamma_tilde_mat=update_gamma_tilde_mat()
      }
    }
    
    if(verbose){
      cat("approximated_log_likelihood =",compute_log_likelihood(theta_hat=NULL,gamma_tilde_mat,phi_vec,fixed_intercept_covariate,sigma2_gamma_vec),"\n") 
    }
    
    # update phi_vec
    update_phi_vec=function(){
      W_U_list=list()
      W_first=list()
      W_second=list()
      gradient_list=list()
      hess_list=list()
      dif=1
      iter=0
      step=1
      
      while(dif>thres){
        before=phi_vec
        d=1
        iter=iter+1
        if (iter>100){
          step=step*0.5
        }
        if (outcome_type[d]=="continuous"){
          
          fixedeffectsdat=cbind(intercept_covariate)
          fixedeffects=c(fixed_intercept_covariate[,d])
          gamma_tilde=gamma_tilde_mat[,d]
          phi=phi_vec[d]
          sigma2_gamma=sigma2_gamma_vec[d]
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=eta
          w=as.matrix(rep(1/phi,nrow(G)))
          W_inv=diag(1/w[,1])
          W=diag(w[,1])
          W_U_list[[d]]=W
          if (approxi){
            gradient_list[[d]]=sum((outcome-mu)^2)/(2*phi^2)-nrow(G)/(2*phi)
            hess_list[[d]]=-sum((outcome-mu)^2)/(phi^3)+nrow(G)/(2*phi^2)
          }else{
            A_1=sigma2_gamma^(-1)*diag(ncol(G))+eigenMapMatMult(eigenMapMatMult(t(G),W),G)
            A_1_inv=chol2inv(chol(A_1))
            W_first[[d]]=diag( rep(-1/(phi^2),nrow(G)) )
            W_second[[d]]=diag( rep(2/(phi^3),nrow(G)) )
            deri_logdet_A1=eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_first[[d]]),G) )
            
            gradient_list[[d]]=sum((outcome-mu)^2)/(2*phi^2)-nrow(G)/(2*phi)-0.5*sum(diag( deri_logdet_A1 ))  
            hess_list[[d]]=-sum((outcome-mu)^2)/(phi^3)+nrow(G)/(2*phi^2) -0.5*( sum(diag( -deri_logdet_A1^2+eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_second[[d]]),G) ) ))  ) 
          }
          
          
        }
        
        if (outcome_type[d]=="binary"){
          # dispersion parameter is equal to 1
        }
        
        if (outcome_type[d]=="count"){
          
          fixedeffectsdat=cbind(intercept_covariate)
          fixedeffects=c(fixed_intercept_covariate[,d])
          gamma_tilde=gamma_tilde_mat[,d]
          phi=phi_vec[d]
          sigma2_gamma=sigma2_gamma_vec[d]
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=exp(eta)
          xi=mu/(phi+mu)
          W=diag( as.numeric((outcome+phi)*xi*(1-xi)) )
          W_U_list[[d]]=W
          if (approxi){
            gradient_list[[d]]=sum( -log(phi+mu)-(outcome+phi)/(phi+mu)+1+log(phi)-(digamma(phi)-digamma(phi+outcome)) )
            hess_list[[d]]=sum( -1/(phi+mu) + (outcome-mu)/((phi+mu)^2) + 1/phi - (psigamma(phi,deriv=1)-psigamma(phi+outcome,deriv=1)) )
          }else{
            A_1=sigma2_gamma^(-1)*diag(ncol(G))+eigenMapMatMult(eigenMapMatMult(t(G),W),G)
            #A_1_inv=MASS::ginv(A_1)
            A_1_inv=chol2inv(chol(A_1))
            W_first[[d]]= diag( as.numeric(xi*(1-xi)+(outcome+phi)*( (mu^2-phi*mu)/((phi+mu)^3) )) )
            W_second[[d]]=diag(as.numeric( 2*(mu^2-phi*mu)/((phi+mu)^3) + (outcome+phi)*(2*phi*mu-4*mu^2)/((phi+mu)^4) )) 
            deri_logdet_A1=eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_first[[d]]),G) )
            gradient_list[[d]]=sum( -log(phi+mu)-(outcome+phi)/(phi+mu)+1+log(phi)-(digamma(phi)-digamma(phi+outcome)) )-0.5*sum(diag( deri_logdet_A1 ))  
            hess_list[[d]]=sum( -1/(phi+mu) + (outcome-mu)/((phi+mu)^2) + 1/phi - (psigamma(phi,deriv=1)-psigamma(phi+outcome,deriv=1)) ) -0.5*( sum(diag( -deri_logdet_A1^2+eigenMapMatMult(A_1_inv,eigenMapMatMult(eigenMapMatMult(t(G),W_second[[d]]),G) ) ))  ) 
          }
          
        }
        
        if (outcome_type[d]=="survival"){
          # no dispersion parameter for survival data
        }
        
        
        
        
        
        d=1
        if (outcome_type[d]=="continuous"){
          phi_vec[d]=phi_vec[d]+step*as.numeric(gradient_list[[d]]/abs(hess_list[[d]]) )
        }
        
        if (outcome_type[d]=="binary"){
          phi_vec[d]=1
        }
        
        if (outcome_type[d]=="count"){
          phi_vec[d]=max(min(phi_vec[d]+step*as.numeric(gradient_list[[d]]/abs(hess_list[[d]]) ),1000),0.001)
        }
        
        if (outcome_type[d]=="survival"){
          # no dispersion parameter
        }
        
        
        
        after=phi_vec
        dif=sum((before-after)^2)
      }
      
      
      
      return(phi_vec)
    }
    
    param="phi_vec"
    if (verbose){
      print(param)
    }
    if (check_difference[param]<num_consecutive){
      phi_vec=update_phi_vec()
    }else{
      separate_iteration[param]=separate_iteration[param]+1
      if (separate_iteration[param] %% num_gap == 0){
        phi_vec=update_phi_vec()
      }
    }
    if(verbose){
      cat("approximated_log_likelihood =",compute_log_likelihood(theta_hat=NULL,gamma_tilde_mat,phi_vec,fixed_intercept_covariate,sigma2_gamma_vec),"\n")  
    }
    
    # update sigma2_gamma_vec
    update_sigma2_gamma_vec=function(){
      d=1
      if (outcome_type[d]=="continuous"){
        
        fixedeffectsdat=cbind(intercept_covariate)
        fixedeffects=c(fixed_intercept_covariate[,d])
        gamma_tilde=gamma_tilde_mat[,d]
        phi=phi_vec[d]
        sigma2_gamma=sigma2_gamma_vec[d]
        dif=1
        while(dif>thres){
          before=sigma2_gamma
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=eta
          w=as.matrix(rep(1/phi,nrow(G)))
          W_inv=diag(1/w[,1])
          W=diag(w[,1])
          #tGWG=eigenMapMatMult(eigenMapMatMult(t(G),W),G)
          Iq=diag(ncol(G))
          Q_inv=eigenMapMatMult(eigenMapMatMult(t(G),W),G)+sigma2_gamma^(-1)*Iq
          Q=MASS::ginv(Q_inv)
          #Q=chol2inv(chol(Q_inv))
          
          A=sigma2_gamma^(-1)*Iq-sigma2_gamma^(-2)*Q
          A2=eigenMapMatMult(A,A)
          #gradient=0.5*( -sum(diag( MASS::ginv(Iq+sigma2_gamma*tGWG)%*%tGWG ))+t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-2) )
          gradient=0.5*( -sum(diag( A )) + t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-2) )  
          #grad(LL_sigma2_gamma,sigma2_gamma)
          #hess=0.5*(sum(diag( (MASS::ginv(Iq+sigma2_gamma*tGWG)%*%tGWG)^2 )) -2*t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-3) )
          hess=0.5*( sum(diag(A2))  -2*t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-3) )
          #hessian(LL_sigma2_gamma,sigma2_gamma)
          
          sigma2_gamma=min(max(sigma2_gamma + as.numeric(gradient/ abs(hess) ),0.001),10)
          
          after=sigma2_gamma
          dif=((before-after)^2)
        }
        sigma2_gamma_vec[d]=sigma2_gamma
      }
      
      
      if (outcome_type[d]=="binary"){
        
        fixedeffectsdat=cbind(intercept_covariate)
        fixedeffects=c(fixed_intercept_covariate[,d])
        gamma_tilde=gamma_tilde_mat[,d]
        phi=phi_vec[d]
        sigma2_gamma=sigma2_gamma_vec[d]
        dif=1
        while(dif>thres){
          before=sigma2_gamma
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=1/(1+exp(-eta))
          w=mu*(1-mu)
          W_inv=diag(1/w[,1])
          W=diag(w[,1])
          #tGWG=eigenMapMatMult(eigenMapMatMult(t(G),W),G)
          Iq=diag(ncol(G))
          Q_inv=eigenMapMatMult(eigenMapMatMult(t(G),W),G)+sigma2_gamma^(-1)*Iq
          #Q=MASS::ginv(Q_inv)
          Q=chol2inv(chol(Q_inv))
          
          A=sigma2_gamma^(-1)*Iq-sigma2_gamma^(-2)*Q
          A2=eigenMapMatMult(A,A)
          #gradient=0.5*( -sum(diag( MASS::ginv(Iq+sigma2_gamma*tGWG)%*%tGWG ))+t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-2) )
          gradient=0.5*( -sum(diag( A )) + t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-2) )  
          #grad(LL_sigma2_gamma,sigma2_gamma)
          #hess=0.5*(sum(diag( (MASS::ginv(Iq+sigma2_gamma*tGWG)%*%tGWG)^2 )) -2*t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-3) )
          hess=0.5*( sum(diag(A2))  -2*t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-3) )
          #hessian(LL_sigma2_gamma,sigma2_gamma)
          sigma2_gamma=min(max(sigma2_gamma + as.numeric(gradient/ abs(hess) ),0.001),10)
          after=sigma2_gamma
          dif=((before-after)^2)
          
        }
        sigma2_gamma_vec[d]=sigma2_gamma
      }
      
      if (outcome_type[d]=="count"){
        
        fixedeffectsdat=cbind(intercept_covariate)
        fixedeffects=c(fixed_intercept_covariate[,d])
        gamma_tilde=gamma_tilde_mat[,d]
        phi=phi_vec[d]
        sigma2_gamma=sigma2_gamma_vec[d]
        dif=1
        while(dif>thres){
          before=sigma2_gamma
          eta=fixedeffectsdat%*%fixedeffects+G%*%gamma_tilde
          mu=exp(eta)
          xi=mu/(phi+mu)
          u_matrix=diag( as.numeric((outcome+phi)*xi*(1-xi)) )
          
          Q_inv=eigenMapMatMult(eigenMapMatMult(t(G),u_matrix),G)+sigma2_gamma^(-1)*diag(ncol(G))
          #Q=MASS::ginv(Q_inv)
          Q=chol2inv(chol(Q_inv))
          gradient=0.5*sum(diag(sigma2_gamma^(-2)*(gamma_tilde%*%t(gamma_tilde)+Q-sigma2_gamma*diag(ncol(G)))))
          # gradient=numDeriv::grad(LL_sigma2_gamma,sigma2_gamma)
          hess=-sum(diag(gamma_tilde%*%t(gamma_tilde)))*sigma2_gamma^(-3)-0.5*ncol(G)
          # hess=numDeriv::hessian(LL_sigma2_gamma,sigma2_gamma)
          
          sigma2_gamma=min(max(sigma2_gamma + as.numeric(gradient/ abs(hess) ),0.001),10)
          after=sigma2_gamma
          dif=((before-after)^2)
        }
        sigma2_gamma_vec[d]=sigma2_gamma
      }
      
      
      if (outcome_type[d]=="survival"){  
        dif=1
        fixedeffectsdat=cbind(intercept_covariate_reorder)
        fixedeffects=c(fixed_intercept_covariate[,d])
        if (is.null(fixedeffectsdat)){
          fixedeffectsdat=matrix(0,nrow=nrow(G)); fixedeffects=0
        }
        gamma_tilde=gamma_tilde_mat[,d]
        sigma2_gamma=sigma2_gamma_vec[d]
        while(dif>thres){
          before=sigma2_gamma
          gradient=0.5*(sigma2_gamma^(-2)*t(gamma_tilde)%*%gamma_tilde-sigma2_gamma^(-1)*ncol(G_reorder) )
          hess=-t(gamma_tilde)%*%gamma_tilde*sigma2_gamma^(-3)+0.5*ncol(G_reorder)*sigma2_gamma^(-2)
          sigma2_gamma=min(max(sigma2_gamma + as.numeric(gradient/ abs(hess) ),0.001),10)
          after=sigma2_gamma
          dif=((before-after)^2)
        }
        sigma2_gamma_vec[d]=sigma2_gamma
      }
      return(sigma2_gamma_vec)
    }
    
    param="sigma2_gamma_vec"
    if (verbose){
      print(param)
    }
    if (check_difference[param]<num_consecutive){
      sigma2_gamma_vec=update_sigma2_gamma_vec()
    }else{
      separate_iteration[param]=separate_iteration[param]+1
      if (separate_iteration[param] %% num_gap == 0){
        sigma2_gamma_vec=update_sigma2_gamma_vec()
      }
    }
    
    if(verbose){
      cat("approximated_log_likelihood =",compute_log_likelihood(theta_hat=NULL,gamma_tilde_mat,phi_vec,fixed_intercept_covariate,sigma2_gamma_vec),"\n") 
    }
    
    new=list("gamma_tilde_mat"=gamma_tilde_mat,"sigma2_gamma_vec"=sigma2_gamma_vec,
             "phi_vec"=phi_vec,"fixed_intercept_covariate"=fixed_intercept_covariate)
    
    for (param in c( "fixed_intercept_covariate", "theta_hat", "gamma_tilde_mat", "phi_vec", "sigma2_gamma_vec") ){
      if (verbose){
        cat(param,"_diff=",sum((old[[param]]-new[[param]])^2),"\n")
      }
      if ( sum((old[[param]]-new[[param]])^2) < thres ){
        check_difference[param]=check_difference[param]+1
      }
    }
    
    
    
    difference=sum((as.numeric( unlist(new) )-as.numeric(unlist(old)) )^2)
    
    if(verbose){
      cat("check_difference=",check_difference,"\n")
      cat("separate_iteration=",separate_iteration,"\n")
      cat("fixed_intercept_covariate=",fixed_intercept_covariate,"\n")
      cat("phi_vec=",phi_vec,"\n")
      cat( "difference=",difference ,"\n")
    }
    
  }
  
  LL_null=compute_log_likelihood(theta_hat=NULL,gamma_tilde_mat,phi_vec,fixed_intercept_covariate,sigma2_gamma_vec)
  
  LL_stat=2*(LL_alt-LL_null)
  pvalue=as.numeric(1-pchisq(LL_stat,df=num_outcome))
  result=list("pvalue"=pvalue,"theta_hat"=theta_hat)
  return(result)
}

