#'Returns prior dose-specific means.
#'
#'Uses elicited efficacy and toxicity dose-specific parameters along with latent prior variance, dose-specific mean hypervariance and frailty variance to determine dose-specific prior means for efficacy and toxicity and prints the prior effective sample size associated with the specified prior parameters.
#'@param PROBST Elicited prior toxicity probability at each dose.
#'@param PROBSE Elicted prior efficacy probability at each dose.
#'@param Var Latent parameter variance for normal probability of efficacy and toxicity.
#'@param HypVar Hypervariance on dose specific mean efficacy and toxicity parameters.
#'@param tau Frailty variance parameter.
#'@param B Number of prior samples to draw for calculating ESS. Suggested values of ten thousand.
#'@return A list contianing the vector of dose-specific efficacy probability prior mean parameters and the vector of dose-specific toxicity probability prior mean parameters.
#'@importFrom mvtnorm rmvnorm
#'@importFrom stats rnorm runif var
#'@examples
#'library(mvtnorm)
#'PROBST=c(.05,.10,.15,.20,.30)
#'PROBSE=c(.2,.4,.6,.65,.7)
#'Var=1
#'HypVar=36
#'tau=1
#'B=100
#'Z=GetPriors(PROBST,PROBSE,Var,HypVar,tau,B)
#'@export
GetPriors=function(PROBST,PROBSE,Var,HypVar,tau,B){
  
  
  
  MeansT= rep(NA,length(PROBST))
  
  for(m in 1:length(MeansT)){
    ##Start at small value
    x=-5
    prob = 0
    while(prob<PROBST[m]){
      prob=1-pnorm(0,x,sqrt(Var+tau))
      x=x+.001
    }
    
    
    MeansT[m]=x
    
    
    
  }
  
  MeansE=MeansT
  
  
  for(m in 1:length(MeansE)){
    ##Start at small value
    x=-5
    prob = 0
    while(prob<PROBSE[m]){
      prob=1-pnorm(0,x,sqrt(Var+tau))
      x=x+.001
    }
    
    
    MeansE[m]=x
    
    
    
  }
  
  
  ##Ok Now lets generate a bunch of probabilities for each dose

  PIT=matrix(rep(NA,B*length(MeansE)),nrow=B)
  PIE=matrix(rep(NA,B*length(MeansE)),nrow=B)
  
  MU1=c(0,0)
  
  
  
  
  for(b in 1:nrow(PIT)){
    ##Generate prior dose means from the hyperprior
    MUE=rep(-10,length(MeansE))
    MUT=MUE
    ##Generate means of first dose
    MUE[1]=rnorm(1,MeansE[1],sqrt(HypVar))
    MUT[1]=rnorm(1,MeansT[1],sqrt(HypVar))
    for(m in 2:length(MeansE)){
      ##Truncate prior dose means
      MUE[m]=pmax(rnorm(1,MeansE[m],sqrt(HypVar)),MUE[m-1])
      MUT[m]=pmax(rnorm(1,MeansT[m],sqrt(HypVar)),MUE[m-1])
    }
    
    
    rho=runif(1,-1,1)
    
    
    for(m in 1:length(MeansE)){
      
      MU1=c(MUE[m],MUT[m])
      
      ##Now we compute the mean for each individual, randomly generate!!
      MU = rmvnorm(1,MU1,
                   matrix(c(Var+tau,tau*rho,tau*rho,Var+tau),byrow=TRUE,nrow=2))
      
      ##Ok now we have our latent vector, compute the EFF and TOX probs
      PIT[b,m]=1-pnorm(0,MU[2],sqrt(Var+tau))
      PIE[b,m]=1-pnorm(0,MU[1],sqrt(Var+tau))
      
      
    }
  }
  
  PriorESS = mean(colMeans(PIT)*(1-colMeans(PIT))/apply(PIT,2,var)) + 
    mean(colMeans(PIE)*(1-colMeans(PIE))/apply(PIE,2,var))-2
  
  
  
  
  print(PriorESS)
  
  
  Z=as.list(c(0,0))
  
  Z[[1]]=MeansE
  Z[[2]]=MeansT
  
  return(Z)
  
  
}




