#' Performs posterior sampling for the MABOUST design and determines whether the trial should continue and what treatment(s) are optimal.
#' @param nCat Number of ordinal outcome categories, i.e. J.
#' @param theta Vector of (J-1)*K specific parameters for the MABOUST model. One row of output from MCMC_MABOUST function.
#' @return Estimated treatment-specific outcome probabilities for a given \eqn{\bf{\theta}} vector.
#' @export
GetProbs=function(nCat,theta){

  gamma=rep(0,nCat-1)
  nTreat=length(theta)/length(gamma)

  PROB=as.list(rep(NA,nTreat))


  theta1=theta[1:length(gamma)]

  for(j in 2:nTreat){
    INDICIES = (length(gamma)*(j-1)+1):(length(gamma)*j)
    theta1=rbind(theta1,theta[INDICIES])
  }

  theta=theta1


  HOLD=rep(NA,length(gamma)+1)
  for(j in 1:length(PROB)){
    m1=1
    HOLD[m1]=exp(gamma[m1]+theta[j,m1])/(1+exp(gamma[m1]+theta[j,m1]))
    m1=length(HOLD)
    HOLD[m1]=1-exp(gamma[m1-1]+theta[j,m1-1])/(1+exp(gamma[m1-1]+theta[j,m1-1]))

    for(m1 in 2:(length(HOLD)-1)){
      HOLD[m1]=exp(gamma[m1]+theta[j,m1])/(1+exp(gamma[m1]+theta[j,m1]))-exp(gamma[m1-1]+theta[j,m1-1])/(1+exp(gamma[m1-1]+theta[j,m1-1]))
    }

    PROB[[j]]=HOLD

  }
  return(PROB)

}


