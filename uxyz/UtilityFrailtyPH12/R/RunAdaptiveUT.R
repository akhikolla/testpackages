#' Simulates replications from the monotone utility based phase 12 trial.
#'
#' Simulates trial replications from the monotone utility based phase 12 trial using either adaptive randomization or fixed dose assignment. Prints the true utility scores, dose selection probability, 
#' average number of patients treated at each dose, average number of responses, average number of toxicities, and Delta value. Returns trial outcomes.
#' @param NSims Number of simulations.
#' @param PE True Efficacy Probability for each dose.
#' @param PT True toxicity probaiblity for each dose.
#' @param corET Correlation parameter between Efficacy and Toxicity status.
#' @param Nmax Maximum Sample size.
#' @param cohort Cohort Size.
#' @param NF Number of fixed assignment patients until adaptive randomization. If NF equals Nmax, the trial is conducted without adaptive randomization.
#' @param UT Utility Matrix with entries U11, U22 elicited and U12 equals 100, U21 equals 0.
#' @param CutE Cutoff for efficacy acceptability.
#' @param CutT Cutoff for toxicity acceptability.
#' @param AcceptE Probability threshold for efficacy acceptability.
#' @param AcceptT Probability threshold for toxicity acceptability.
#' @param HypermeansE Dose-specific hypermeans for efficacy.
#' @param HypermeansT Dose-specific hypermeans for toxcity.
#' @param Hypervars Length 3 vector of hypervariances. Hypervars[1] contains the Latent parameter variance for normal probability of efficacy and toxicity. Hypervars[2] contains the hypervariance on dose specific mean efficacy and toxicity parameters and Hypervars[3] contains the frailty variance parameter.
#' @return A list of size NSims with results from each simulated trial. Each entry contains a list with (1) the optimal dose selected, (2) the posterior mean utility for each dose level, (3) a matrix containing the dose given, the efficacy outcome and the toxicity outcome for each patient.
#' @importFrom stats pnorm sd
#' @importFrom bindata rmvbin
#' @importFrom mvtnorm pmvnorm
#' @examples
#' library(bindata)
#' library(mvtnorm)
#' ##Trial PArameters here
#' Nmax=30 ##Number of patients to enroll
#' NF=30 ##Number until AR if NF=Nmax, there's no AR.
#' cohort=3
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
#' #True Efficacy and Toxicity probabilities
#' PE=c(.2,.4,.6,.7,.7)
#' PT=c(.2,.2,.2,.3,.5)
#' corET=0
#' ##Number of simulations
#' NSims=2
#' RESULTS=RunAdaptiveUT(NSims, PE, PT, corET,  Nmax, 
#' cohort,NF, UT, CutE, CutT, AcceptE, 
#' AcceptT, HypermeansE, HypermeansT, Hypervars)
#'@export
RunAdaptiveUT=function(
  NSims, ##Number of simulations
  PE, ##True Efficacy Probability for each dose
  PT, ##True toxicity probaiblity for each dose
  corET, ##Correlation parameter between Eff and Tox
  Nmax, ##MAximum Sample size
  cohort, ##Cohort Size
  NF, ##Number of fixed assignment patients until adaptive randomization
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
  
Hypervars=c(Hypervars[1],Hypervars[1],Hypervars[2],Hypervars[2],Hypervars[3])
  
  TRUEUT=rep(NA,nDose)
  
  
  for(k in 1:nDose){
    Z1 =matrix(c( (1-PT[k])*(1-PE[k]), PE[k]*(1-PT[k]),
                  PT[k]*(1-PE[k]), PT[k]*PE[k]), byrow=TRUE, nrow=2)
    
    TRUEUT[k]=sum(Z1*UT)
    
  }
  
  
  
  NumTrt=matrix(rep(NA,nDose*NSims),nrow=NSims)
  
  CORMAT <- matrix(c(1,corET,corET,1), ncol=2)   
  
  ##Setup Simulation parameters
  SIMSTORE = as.list(rep(0,NSims))
  
  DoseOpt=rep(NA,NSims)
  NTox = rep(NA,NSims)
  NEff=rep(NA,NSims)
  
  DoseStore=rep(NA,Nmax)
  YEStore = rep(NA,Nmax)
  YTStore =rep(NA,Nmax)
  
  B=2000
  
  
  ##No correlation
  for(m1 in 1:NSims){
    
    if(m1%%1000==0){
      
      cat(paste(m1,"Simulations
                ",sep="  "))
      
    }
    
    
    Dosetried=c(1,rep(0,nDose-1))
    
    ##Start at lowest dose
    Doses=rep(1,cohort)
    
    
    OUT <- rmvbin(cohort, margprob = c(PE[1], PT[1]), bincorr = CORMAT) 
    
    YE=OUT[,1]
    YT=OUT[,2]
    
    
    ##Ok now we have the starting values
    
    
    
    
    
    for(i in 2:(Nmax/cohort)){
      
      
      
      
      Z=UTEFFTOX(YE,YT, Doses, HypermeansE, HypermeansT,  Hypervars, B )
      
      
      ##   plot(1:nrow(Z[[5]]),Z[[5]],type="l")
      
      ##     print(length(unique(Z[[5]]))/length(Z[[5]]))
      
      ## print((Z[[5]]))
      
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
      
      
      if(sum(ACCEPTE*ACCEPTT)==0 && i>2){
        OptDose=0
        
        break
        
        
      }else{
        
        
        ##Multiply MeanUT by accept
        
        MeanUT= MeanUT*ACCEPTE*ACCEPTT
        
        if(sum(MeanUT)==0){
          OptDose=1
        }else{
          
          if((i*cohort)<=NF){
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
        
        
        ##Now enroll the next three patients
        Doses=c(Doses,rep(OptDose,cohort))
        #What are their binary outcomes
        OUT <- rmvbin(cohort, margprob = c(PE[OptDose], PT[OptDose]), bincorr = CORMAT) 
        
        YE=c(YE,OUT[,1])
        YT=c(YT,OUT[,2])
        
        
        
        
      }
      
      
      
    }
    
    if(OptDose==0){
      NSkip=NSkip+1
      ##This simulation doesn't count now.
      
      X1=as.list(c(0,0,0))
      X1[[1]]=OptDose
      X1[[2]]=MeanUT
      X1[[3]]=cbind(Doses,YE,YT)
      
      SIMSTORE[[m1]]=X1
      NTox[m1]=sum(YT)
      NEff[m1]=sum(YE)
      DoseOpt[m1]=OptDose
      for(j in 1:nDose){
        NumTrt[m1,j]=sum(Doses==j)
      }
      
      
    }else{
      
      ##Ok now weve done the whole trial, whats the optimal dose
      
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
      
      
      ##Multiply MeanUT by accept
      
      if(sum(ACCEPTE*ACCEPTT)==0){
        OptDose=0
        
        
      }else{
        MeanUT= MeanUT*ACCEPTE*ACCEPTT
        
        ##Assign Next Dose Deterministically
        OptDose = which(MeanUT==max(MeanUT))
      }
      
      ##Fill this list with our values
      X1=as.list(c(0,0,0))
      X1[[1]]=OptDose
      X1[[2]]=MeanUT
      X1[[3]]=cbind(Doses,YE,YT)
      
      SIMSTORE[[m1]]=X1
      
      NTox[m1]=sum(YT)
      NEff[m1]=sum(YE)
      
      for(j in 1:nDose){
        NumTrt[m1,j]=sum(Doses==j)
      }
      
      DoseOpt[m1]=OptDose
      
      
    }
    
    
    
  }
  
  
  ##Calculate True Accept
  TRUECEPT = (PE>=CutE)*(PT<=CutT)
  
  
  
  cat("True Mean Utility Scores
      
      ")  
  
  print(TRUEUT)
  
  cat("Dose Selection Probability
      
      ")
  print(table(DoseOpt)/NSims)
  
  
  cat("Average Number of patients treated at each dose
      
      ")
  print(colMeans(NumTrt))
  
  cat("Mean Number of Responses
      
      ")
  print(mean(NEff))
  
  
  cat("Mean Number of Toxicities
      
      ")
  print(mean(NTox))
  
  
  prob1=rep(0,nDose+1)
  
  for(j in 0:nDose){
    prob1[j]=mean(DoseOpt==j)
  }
  
  U1=TRUEUT
  
 # U1=(TRUEUT-mean(TRUEUT))/sd(TRUEUT)
  U1=c(0,U1*TRUECEPT)
  cat("Delta U
      
      ")
  
  print(sum(prob1*(max(U1)-U1)))
  
  
  return(SIMSTORE)
  
  
  
  
  
  }