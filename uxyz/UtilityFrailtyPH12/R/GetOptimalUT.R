#' Gives the dose to assign the next patient cohort using monotone utility based phase 12 trial.
#'
#' Gives the dose to assign the next patient cohort using the monotone utility based phase 12 trial using either adaptive randomization or fixed dose assignment. If the trial has been completed, AR equals FALSE will give the optimal dose level for the trial.
#' @param YE Vector of patient efficacy status.
#' @param YT Vector of patient toxicity status.
#' @param Doses Vector of patient dose assignment.
#' @param DoseTried Vector containing 1s for doses tried and 0 otherwise.
#' @param AR Logical stating whether or not to adaptively randomize the next cohort of patients. If the trial has been completed, AR equals FALSE will give the optimal dose level for the trial.
#' @param UT Utility Matrix with entries U11, U22 elicited and U12 equals 100, U21 equals 0.
#' @param CutE Cutoff for efficacy acceptability.
#' @param CutT Cutoff for toxicity acceptability.
#' @param AcceptE Probability threshold for efficacy acceptability.
#' @param AcceptT Probability threshold for toxicity acceptability.
#' @param HypermeansE Dose-specific hypermeans for efficacy.
#' @param HypermeansT Dose-specific hypermeans for toxcity.
#' @param Hypervars Length 3 vector of hypervariances. Hypervars[1] contains the Latent parameter variance for normal probability of efficacy and toxicity. Hypervars[2] contains the hypervariance on dose specific mean efficacy and toxicity parameters and Hypervars[3] contains the frailty variance parameter.
#' @return A numerical value of the dose to assign the next patient cohort to. If the trial has been completed, this is the optimal dose. If a value of 0 is returned, no doses are acceptable and the trial should be stopped.
#' @examples
#' library(mvtnorm)
#' ##Data Here
#' YE=rbinom(30,1,.8)
#' YT=rbinom(30,1,.3)
#' Doses=sample(1:3,30,replace=TRUE)
#' DoseTried=c(1,1,1,0,0)
#' ##UTILITIES
#' UT = matrix(c(38.23529,100,0,61.76471),nrow=2,byrow=TRUE)
#' ##Safety Parameters
#' CutE=.3
#' CutT=.4
#' AcceptE=.1
#' AcceptT=.1
#' ##Hyperparameters for Utility
#' HypermeansE=c(-1.189, -0.357,  0.360,  0.546,  0.743)
#' HypermeansT=c(-2.325, -1.811, -1.464, -1.189, -0.740)
#' Hypervars=c(1,36,1)  
#' ##Adaptively randomize or not?
#' AR=FALSE
#' GetOptimalUT( YE,YT, Doses,DoseTried,AR, UT, CutE, CutT,
#' AcceptE, AcceptT, HypermeansE, HypermeansT,Hypervars)
#'@export
GetOptimalUT=function(
  YE, ##Vector of patient efficacy status.
  YT, ##Vector of patient toxicity status.
  Doses, ##Vector of patient dose assignment.
  DoseTried, ##Vector containing 1s for doses tried and 0 otherwise.
  AR, ##Logical value on whether or not to adaptively randomize.
  UT, ##Utility Matrix
  CutE, ##Cutoff For efficacy acceptability
  CutT, ##Cutoff for toxicity acceptability
  AcceptE, ##Probability threshold for eff acceptability
  AcceptT, ##Probability threshold for tox acceptability
  HypermeansE, ##Hypermeans for efficacy
  HypermeansT, ##Hypermeans for Toxcity
  Hypervars ##Hypervariances
){
  NSkip=0  ##Will  count the number of skipped simulations
  
  nDose=length(HypermeansE)
  Dosetried=DoseTried
  
  Hypervars=c(Hypervars[1],Hypervars[1],Hypervars[2],Hypervars[2],Hypervars[3])
  
  
  B=2000

      
      Z=UTEFFTOX(YE,YT, Doses, HypermeansE, HypermeansT,  Hypervars, B )
      

      Sigma=matrix(c(Hypervars[1]+Hypervars[5],Hypervars[5]*mean(Z[[3]]),
                     Hypervars[5]*mean(Z[[3]]), Hypervars[2]+Hypervars[5]), nrow=2, byrow=TRUE)
      PMAT=Sigma
      
      
      MeanUT=rep(NA,nDose)
      
      for(D in 1:nDose){
        MU=c(mean(Z[[1]][,D]),mean(Z[[2]][,D]))
        
        
        ##YE,YT=0,0
        lower <- rep(-Inf, 2)
        upper <- rep(0, 2)
        PMAT[1,1] <- pmvnorm(lower, upper, MU, sigma=Sigma)
        ##YE,YT=0,1
        lower <- c(-Inf,0)
        upper <- c(0,Inf)
        PMAT[2,1] <- pmvnorm(lower, upper, MU, sigma=Sigma)
        ##YE,YT = 1,0
        lower <- c(0,-Inf)
        upper <- c(Inf,0)
        PMAT[1,2] <- pmvnorm(lower, upper, MU, sigma=Sigma)
        ##YE, YT = 1,1
        lower <- c(0,0)
        upper <- c(Inf,Inf)
        PMAT[2,2] <- pmvnorm(lower, upper, MU, sigma=Sigma)
        
        MeanUT[D]=sum(UT*PMAT)
        
        
      }
      
      ACCEPTE=rep(0,nDose)
      ACCEPTT=rep(0,nDose)
      ##What Doses are acceptable?
      for(D in 1:nDose){
        ##Calculate the probability of eff over all the samples
        PEFF = 1-pnorm(0,Z[[1]][,D],sqrt(Hypervars[1]+Hypervars[5]))
        
        ##IS THIS DOSE ACCEPTABLE IN TERMS OF EFFICACY?
        ACCEPTE[D]=mean(PEFF>CutE)>AcceptE
        
        ##Calculate the probability of tox over all the samples
        PTOX = 1-pnorm(0,Z[[2]][,D],sqrt(Hypervars[2]+Hypervars[5]))
        
        ##IS THIS DOSE ACCEPTABLE IN TERM OF TOXICITY?
        ACCEPTT[D]=mean(PTOX<CutT)>AcceptT
        
      }
      
      
      if(sum(ACCEPTE*ACCEPTT)==0){
        OptDose=0
        

        
      }else{
        
        
        ##Multiply MeanUT by accept
        
        MeanUT= MeanUT*ACCEPTE*ACCEPTT
          
          if(AR==FALSE){
            ##Assign Next Dose Deterministically
            if(sum(Dosetried>0)<nDose){
              ##We haven't tried all the doses
              Dosetried1=Dosetried
              Dosetried1[max(which(Dosetried>0))+1]=1
              Dosetried1 = Dosetried1>0
            }
            MeanUT=MeanUT*Dosetried1
            OptDose = which(MeanUT==max(MeanUT))
            
            Dosetried[OptDose]=1
            
            
            
          }else{
            ##Adaptive randomization
            if(sum(Dosetried>0)<nDose){
              ##We haven't tried all the doses
              Dosetried1=Dosetried
              Dosetried1[max(which(Dosetried>0))+1]=1
              Dosetried1 = Dosetried1>0
            }
            MeanUT=MeanUT*Dosetried1
            
            
            Indexes = which(MeanUT>0)
            MeanUT1=MeanUT[MeanUT>0]
            
            if(length(Indexes)>1){
              
              U1= MeanUT1/10
              ##Calculate ranodmization probabilities
              probs=exp(U1)/sum(exp(U1))
              OptDose=sample(Indexes,1,FALSE, prob=probs)
            }else{
              
              OptDose=Indexes
            }
            
            
            Dosetried[OptDose]=1
            
            
            
          }
          
          
        
   

    
  }
  
      
      cat("Posterior Mean Utitliy Scores
          ")
      print(MeanUT)
      
      
      cat("OPTIMAL Dose
          ")
      print(OptDose)
      
      return(OptDose)
}
  