profile_minmax_joint_posneg_KT <-
function(pars,listdata,u,q1=0,
                                    q2=1,...,sill=-10^(40))
  {
    loglik_min     <- NULL
    loglik_max     <- NULL
    loglik_neg_min <- NULL
    loglik_neg_max <- NULL    
    loglik         <- NULL
    
    loglik_min <- Profile_likelihood_cd_nm_joint_D_KT(par=pars,
                                                      listr=listdata,
                                                      x=u,Zestfun=quantile,
                                                      probs=q1,          
                                                      silly=sill,...)
    
    loglik_max <- Profile_likelihood_cd_nm_joint_D_KT(par=pars,
                                                      listr=listdata,
                                                      x=u,Zestfun=quantile,
                                                      probs=q2,
                                                      silly=sill,...)

    loglik_neg_min <- Profile_likelihood_cd_nm_joint_D_KT_neg(par=pars,
                                                             listr=listdata,
                                                              x=u,Zestfun=quantile,
                                                              probs=q1,
                                                              silly=sill,...)

    loglik_neg_max <- Profile_likelihood_cd_nm_joint_D_KT_neg(par=pars,
                                                              listr=listdata,
                                                              x=u,Zestfun=quantile,
                                                              probs=q2,
                                                              silly=sill,...)
    
    
    if(loglik_min == sill || loglik_max == sill ||
       loglik_neg_min == sill || loglik_neg_max == sill)
      {
        return( sill )
      }
    if(loglik_min != sill  & loglik_max != sill &
       loglik_neg_min != sill & loglik_neg_max != sill)
      {
        return( loglik_max )
      }
  }
