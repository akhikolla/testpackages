Profile_likelihood_HT_unc <-
function (par,listr,x,silly=-10^(40))
{
  n                <- NULL
  sig              <- NULL
  sumX             <- NULL
  temp             <- NULL
  temp2            <- NULL
  z                <- list()            
  Pl               <- silly
  X                <- vector('list',length(listr))
  Y                <- vector('list',length(listr))
  Z                <- vector('list',length(listr))
  index_alpha      <- seq(1,((2*(length(listr)) ) -1),by=2)
  index_beta       <- seq(2,((2*(length(listr)) )   ),by=2)
  alpha            <- par[index_alpha]
  beta             <- par[index_beta]    
  Z                <- vector('list',length(listr))
  Zstar            <- vector('list',length(listr))
  
   
  for(i in 1:length(listr))
    {
      temp           <- as.matrix(listr[[i]])
      X[[i]]         <- temp[,1][temp[,1]>x]
      n[i]           <- length(X[[i]])
      Y[[i]]         <- temp[,2][temp[,1]>x]
      Z[[i]]         <- (Y[[i]]  - alpha[i]*X[[i]])/(X[[i]]^beta[i])
      Zstar[[i]]     <- (Y[[i]]  - X[[i]])
      sig[i]         <- (1/n[i]) * sum ((Z[[i]]-mean(Z[[i]]))^2)
      sumX[i]        <- sum(beta[i]*log(X[[i]]))
    }
  
  if(all(alpha <= 1) & all(alpha >= -1) & all(beta < 1) ) 
    {                                       
      Pl  <- sum(((-(n/2)*log (2*pi*sig)) - sumX - (n/2)))      
    }
  if((all(alpha <= 1) ==FALSE) ||  (all(alpha >= -1)==FALSE) ||
     (all(beta < 1)==FALSE) )
    {
      Pl <- silly
    }  
  z$Pl <- Pl    
  return(z$Pl)
}
