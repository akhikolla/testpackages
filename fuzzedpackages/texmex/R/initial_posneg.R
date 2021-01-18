initial_posneg <-
function(D,...)
  {    
    a    <- runif((1000*D),-1,1)
    b    <- runif((1000*D),-3,0.99)
    prop <- matrix(rbind(a,b),nrow=2*D) 
        n <- 1000
        j <- 1
        Pl <- profile_minmax_joint_posneg_KT(pars=prop[,j],...)
    while(Pl<=-10^10)
      {
        j <- j+1
        if(j<=1000)
          {
            Pl <-  profile_minmax_joint_posneg_KT(pars=prop[,j],...)
          }
        if(j>1000)break        
      }
    if(j<=1000)
      {
        return(prop[,j])
      }
  }
