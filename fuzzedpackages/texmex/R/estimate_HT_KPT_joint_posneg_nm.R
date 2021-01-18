estimate_HT_KPT_joint_posneg_nm <-
function(pars,x,listr,params=TRUE,...,k=3)
  {    
    temp <- NULL
    temp <- optim(par=pars,profile_minmax_joint_posneg_KT,
                  listdata=listr,u=x,...,
                  method="Nelder-Mead",
                  control=list(fnscale=-1,maxit=100000))

    tempar <- temp$par
    for(j in 1:k)
      {
        temp <- optim(par=tempar,profile_minmax_joint_posneg_KT,
                      listdata=listr,u=x,...,
                      method="Nelder-Mead",
                      control=list(fnscale=-1,maxit=100000))
        tempar <- temp$par
      }

    ifelse(params==TRUE,return(tempar), return(temp))
  }
