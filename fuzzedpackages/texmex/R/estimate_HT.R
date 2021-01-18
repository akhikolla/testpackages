estimate_HT <-
function(list,u,pars,params=TRUE)
                                          
  {
    res  <- optim(par=pars,Profile_likelihood_HT_unc,
                  listr=list,x=u,
                  control=list(fnscale=-1,maxit=100000))        
    ifelse(params==TRUE,return(res$par), return(res))
  }
