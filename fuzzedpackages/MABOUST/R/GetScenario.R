#' Returns the superiority or futility cutoff during a MABOUST trial.
#' @param nTreat Number of active treatments in consideration, i.e. 1,...,K.
#' @param RANGES J-list containing ranges of plausible marginal treatment outcome probabilities.
#' @param RANGES1 J-list containing ranges of plausible covariate adjusted outcome probabilities.
#' @param XPROB List of matrices containing discrete values various covariates can take, along with their probabilities.
#' @return Randomly generate marginal ordinal outcome probabilities for each treatment and a covariate vector.
#' @references
#' [1] Chapple and Clement (2020), MABOUST: A Multi-Armed Bayesian Ordinal Outcome Utility-Based Sequential Trial. Submitted.
#' @importFrom stats runif
#' @examples
#' ###Trial parameters
#' nTreat = 3
#' nCat=6
#' ###Marginal Probability Ranges
#' RANGES = as.list(rep(NA,nCat))
#' RANGES[[1]]=c(.1,.35)
#' RANGES[[2]]=c(.1,.3)
#' RANGES[[3]]=c(.4,.7)
#' RANGES[[4]]=c(0,.1)
#' RANGES[[5]]=c(.1,.3)
#' RANGES[[6]]=c(.0,.1)
#' ###Covariate Adjusted Probability Ranges
#' RANGES1=RANGES
#' RANGES1[[1]]=c(0,.5)
#' RANGES1[[2]]=c(0,.5)
#' RANGES1[[3]]=c(0,.8)
#' RANGES1[[4]]=c(0,.45)
#' RANGES1[[5]]=c(0,.45)
#' RANGES1[[6]]=c(0,.30)
#' XPROB = as.list(rep(NA,3))
#' XPROB[[1]]=rbind(0:10,round(dpois(0:10,2),2)) ###CCI
#' XPROB[[2]]=rbind(c(-1,0,1),c(.5,.4,.1)) ###O2 Status
#' XPROB[[3]]=rbind(c(-2,-1,0,1),c(.27,.38,.18,.17))
#' GetScenario(nTreat,RANGES,RANGES1, XPROB)
#' @export
GetScenario = function(nTreat,
                       RANGES,
                       RANGES1,
                       XPROB){
  nCat=length(RANGES)

  THETAPROB = matrix(nrow=nTreat,ncol=length(RANGES))

  PROBS = as.list(rep(NA,nTreat))


  for(j in 1:nTreat){

    X=rep(NA,length(RANGES))
    X[1]=runif(1,RANGES[[1]][1],RANGES[[1]][2])
    count = 1
    while(count<nCat){
      count=1
      for(k in 2:(nCat-1)){
        X[k]=runif(1,0,1)*(1-sum(X[1:(k-1)]))
        if((X[k]>RANGES[[k]][1])&(X[k]<RANGES[[k]][2])){
          count=count+1
        }
      }

      k=nCat
      X[k]=(1-sum(X[1:(k-1)]))
      if((X[k]>RANGES[[k]][1])&(X[k]<RANGES[[k]][2])){
        count=count+1
      }


    }

    PROBS[[j]]=X
    THETAPROB[j,]=X

  }

  PROBHOLD= THETAPROB




  n=10000

  X = matrix(nrow=n,ncol=length(XPROB))

  for(j in 1:ncol(X)){
    X[,j]=sample(XPROB[[j]][1,],n,prob=XPROB[[j]][2,],replace=TRUE)
  }


  T1=sample(1:nTreat,n,replace=T)


  PROBHOLD = matrix(rep(-10,n*nCat),ncol=nCat)
  count = -10
  m=0
  while(count<nCat){

    STOP=0
    ###Get random Age/CCI/O2 Effects



    OREFF = runif(ncol(X),-2,0)

    for(i in 1:n){
      prob = PROBS[[T1[[i]]]]
      cumprob = cumsum(prob)[-length(prob)]
      ETA = log(cumprob/(1-cumprob))
      ETA = ETA + c(as.vector(X[i,])%*%OREFF)
      cumprob = c(0,exp(ETA)/(1+exp(ETA)),1)
      prob=diff(cumprob)
      PROBHOLD[i,]=prob
      for(j in 1:length(RANGES1)){
        if(prob[j]< RANGES1[[j]][1] | prob[j]>RANGES1[[j]][2]){
          STOP=1
          break   ###One observation violated the conditions
        }
      }

      if(STOP==1){
        break
      }

    }


    ###Made it through...
    if(STOP==0){
      count=10000
    }





  }


  HOLD=as.list(rep(NA,3))

  rownames(THETAPROB)=paste0("Treatment ",1:nTreat )
  colnames(THETAPROB)=paste0("Outcome ",1:ncol(THETAPROB))

  HOLD[[1]]=THETAPROB
  HOLD[[2]]=OREFF
  HOLD[[3]]=XPROB

  names(HOLD)=c("Marginal Outcome Probabilities",
                "Beta",
                "Covariate Distribution")


  return(HOLD)

}
