Profile_likelihood_cd_nm_joint_D_KT_neg <-
function (par,listr,x,
                                                     Zestfun,...,v,
                                                     silly=-10^(40))
  {
    n                <- NULL
    sig              <- NULL
    sumX             <- NULL
    no_of_roots      <- NULL
    no_of_roots_star <- NULL
    temp             <- NULL
    Zq               <- NULL
    Zstarq           <- NULL
    xstar            <- NULL
    xdstar           <- NULL
    s                <- NULL
    cond_alphas      <- NULL
    cond_ord_dep     <- NULL
    cond_ord_pairs   <- NULL
    vdep             <- NULL
    z                <- list()            
    Pl               <- silly
    X                <- vector('list',length(listr))
    Y                <- vector('list',length(listr))
    Z                <- vector('list',length(listr))
    Zstar            <- vector('list',length(listr))
    index_alpha      <- seq(1,((2*(length(listr)) ) -1),by=2)
    index_beta       <- seq(2,((2*(length(listr)) )   ),by=2)
    alpha            <- par[index_alpha]
    beta             <- par[index_beta]
    xstar            <- rep(v,(length(listr)-1))
    xdstar           <- rep(v,(length(listr)-1))
    xdepstar         <- rep(vdep,length(listr))
    
    for(i in 1:length(listr))
      {
        cond_alphas[i] <- ((alpha[i]) <= 1)        
        temp       <- as.matrix(listr[[i]])
        X[[i]]     <- temp[,1][temp[,1]>x]
        vdep[i]    <- max(X[[i]])
        n[i]       <- length(X[[i]])
        Y[[i]]     <- temp[,2][temp[,1]>x]
        Z[[i]]     <- (Y[[i]] - alpha[i]*X[[i]])/(X[[i]]^beta[i])
        Zstar[[i]] <- (Y[[i]] + X[[i]])
        Zq[i]      <- Zestfun(Z[[i]],...)        
        Zstarq[i]  <- Zestfun(Zstar[[i]],...)
        sig[i]     <- (1/n[i]) * sum ((Z[[i]]-mean(Z[[i]]))^2)
        sumX[i]    <- sum(beta[i]*log(X[[i]]))
      }
    
    if(all(cond_alphas==TRUE))
      {
        for(i in 1:length(listr))
          {
            temp_roots_star      <- roots(lev=v,a=alpha[i],c=-1,
                                          b=beta[i],d=0,Zj=Zq[i],            
                                          Zk=Zstarq[i])        
            xdepstar[i]          <- temp_roots_star$xstar            
          }
      }
    
    if(all(alpha <= 1) & all(alpha >= -1) & all(beta <= 1) &
       all(cond_alphas==TRUE)) 
      {                
        for(j in 1:length(listr))
          {
            cond_ord_dep[j] <-  ( -1 <=   # (alpha[j])
                                 #mark:change v to vdep[j]
                                 (min(alpha[j],(Dcond(v,alpha[j],beta[j],-1,0,
                                               Zq[j],Zstarq[j])),#-(1e-10)),
                                      (Dcond(xdepstar[j],alpha[j],beta[j],-1,
                                             0,Zq[j],
                                             Zstarq[j])))))#-(1e-10)) )))            
          }
        
        condition <- (all(cond_ord_dep==TRUE))
        
        if( condition == TRUE )
          {
            Pl  <- sum(((-(n/2)*log (2*pi*sig)) - sumX - (n/2)))# note sig is actually sigmaSquared in the normal density
          }    
        if(condition == FALSE)
          {
            Pl  <- silly
          }          
      }
    if((all(alpha <= 1) ==FALSE) ||  (all(alpha >= -1)==FALSE) ||
       (all(beta < 1)==FALSE) || (all(cond_alphas==TRUE)==FALSE))
      {
        Pl <- silly
      }
    
     z$Pl <- Pl
#    z$Zs <- Zstarq
#    z$Zq <- Zq
      
    return(z$Pl)
  }
