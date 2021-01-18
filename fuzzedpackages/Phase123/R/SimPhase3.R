#' Performs one replication of phase 3 for the phase 123 design, given phase 12 data.
#'
#' This function simulates the phase 3 potion of the phase 123 trial, given phase 12 outcomes.
#' @param Dose Vector of standardized doses considered in the trial.
#' @param PE True efficacy dose-toxicity vector.
#' @param PT True toxicity dose-toxicity vector.
#' @param Phase12 Matrix Consisting of patient data from a phase 12 trial. The columns are in order: Doses given, YE, YT, Accrual Times
#' @param Hypermeans Prior Means for the Eff-Tox design of length 6.
#' @param Hypervars Prior Variances for the Eff-Tox design of length 6.
#' @param betaA True linear term for the rate or mean parameter (beta_1,exp(beta_E),-exp(beta_T),beta_2,beta_0) for agent A.
#' @param ProbC Probability of efficacy and toxicity for the control therapy.
#' @param betaC Linear term for efficacy, toxicity and beta_0 for the control group.
#' @param Family Time to event distribution. Options include: Exponential, Gamma, Weibull, Lognormal.
#' @param alpha Shape parameter or standard deviation of a lognormal distribution.
#' @param Nmax Maximum number of patients to enroll in phase 3.
#' @param Opt Dose used for A to begin randomization in phase 3.
#' @param Accrue Accrual rate for patients in phase 3.
#' @param Time12 Time window for phase 12.
#' @param Twait Waiting time in between phase 12 and phase 3.
#' @param NLookSwitch Number of patient events to determine if we re-optimize doses for A.
#' @param NLook Vector of information criteria for making interim looks.
#' @param Sup Vector of superiority boundaries.
#' @param Fut Vector of futility boundaries.
#' @importFrom  stats sd rbinom rexp rgamma rlnorm rweibull
#' @importFrom survival coxph
#' @importFrom Rcpp evalCpp
#' @importFrom stats quantile
#' @useDynLib Phase123
#' @references
#' [1] Chapple and Thall (2018).A Hybrid Phase I-II/III Clinical Trial Design Allowing Dose Re-Optimization in Phase III. Biometrics. In Press,
#' @examples
#' library(survival)
#'##True Efficacy and Toxicity Probabilities
#'PT = c(.1,.15,.25,.35,.5)
#'PE=c(.2,.4,.6,.65,.7)
#'##Dose Levels considered
#'Dose = c(1,2,3,3.5,5)
#'Dose=(Dose-mean(Dose))/sd(Dose)
#'##Average accrual rate for phase III
#'Accrue = 10
#'#'##Hypermeans for Eff-Tox
#'Hypermeans = c(.022,3.45,0,-4.23,3.1,0)
#'Hypervars = c(2.6761, 2.6852, .2, 3.1304, 3.1165, 1)
#'Hypervars=Hypervars^2
#'Contour = c(.35, .75,.7,.4)
#'PiLim = c(.3,.4)
#'ProbLim=c(.1,.1)
#'###Family of Distributions
#'Family="Exponential"
#'###Shape parameter ## Doesn't matter for exponential distribution
#'alpha=1
#'###True Beta vector
#'betaA = c(.75,-.5, .3, -.25,2.143)
#'##True beta vector for efficacy, toxicity and intercept of the control treatment
#'betaC=c(.3,-.25,2.389)
#'##True efficacy and toxicity probability for control group
#'ProbC = c(.4,.15)
#'##Waiting time in between
#'Twait=1
#'###How long is the time window in phase 12?
#'Time12=1
#'##Dose to start phase 3 with
#'Opt=3
#'##Make matrix with old phase 12 data
#'Doses= c(1,1,1,2,2,2,1,1,1,3,3,3,1,1,1,2,2,2)
#'YE = c(0,0,1,1,1,0,0,0,0,1,1,1,0,0,1,1,1,0)
#'YT=c(0,0,0,1,1,0,1,0,0,1,1,1,0,0,0,1,0,0)
#'##Accrual Times for old data
#'Accrue12=2
#'##Size of phase 12 cohort
#'cohort=3
#'ACC1=cumsum(rexp(length(YT),Accrue12))
#'##Accrual times are the same for each cohort in phase 12
#'Grab = rep(NA,length(YT)/cohort)
#'for(m in 1:length(Grab)){Grab[m]=ACC1[m*3]}
#'for(m in 1:length(Grab)){ACC1[((m-1)*cohort+1):((m-1)*cohort+cohort)]=rep(Grab[m],cohort)}
#'Phase12 = cbind(Doses,YE,YT,ACC1)
#'betaC=c(.3,-.25,2.389)
#'##True efficacy and toxicity probability for control group
#'ProbC = c(.4,.15)
#'##Max Sample Size
#'Nmax=500
#'###Number of patient events to Re-optimize doses
#'NLookSwitch = 50
#'##Number of patient events for interim looks
#'NLook = c(200,300,400)
#'##Superiority Boundaries
#'Sup = c(2.96, 2.53,1.99)
#'##Futility Boundaries (0 means no futility decision)
#'Fut = c(0,1.001,0)
#'##Starting Dose, hat(x)_ET
#'Opt=3
#'##Number of simulations to run
#'nSims=10
#'SimPhase3(Dose,Phase12,PE,PT,Hypermeans,Hypervars,betaA,
#'ProbC,betaC,Family,alpha,Nmax,Opt,Accrue,
#'Time12,Twait,NLookSwitch,NLook,Sup,Fut)
#' @export
SimPhase3=function(Dose,Phase12,PE,PT,Hypermeans,Hypervars,betaA,ProbC,betaC,Family,alpha,Nmax,Opt,Accrue,Time12,Twait,NLookSwitch,NLook,Sup,Fut){


  B=2e3

  ProbE=PE
  ProbT=PT
  probE=PE
  probT=PT

  YE=Phase12[,2]
  YT=Phase12[,3]

  Nmax1=Nmax

##Generate Survival Times for patients from phase 12
  NET=nrow(Phase12)
  TIMES=rep(NA,NET)

  Doses=Dose[Phase12[,1]]


    if(Family=="Gamma"){

      for(h in 1:NET){
    TIMES[h]=rgamma(1,alpha,1/exp(betaA[1]*Doses[h]+betaA[2]*YE[h]+betaA[3]*YT[h]+betaA[4]*Doses[h]^2+betaA[5]))
    }
    }


    if(Family=="Weibull"){
      for(h in 1:NET){
      TIMES[h]=rgamma(1,alpha,exp(betaA[1]*Doses[h]+betaA[2]*YE[h]+betaA[3]*YT[h]+betaA[4]*Doses[h]^2+betaA[5]))
      }
    }

    if(Family=="Lognormal"){
      for(h in 1:NET){
      TIMES[h]=rlnorm(1,betaA[1]*Doses[h]+betaA[2]*YE[h]+betaA[3]*YT[h]+betaA[4]*Doses[h]^2+betaA[5],alpha)
      }
    }

  if(Family=="Exponential"){
    for(h in 1:NET){
      TIMES[h]=rexp(1,1/exp(betaA[1]*Doses[h]+betaA[2]*YE[h]+betaA[3]*YT[h]+betaA[4]*Doses[h]^2+betaA[5]))
    }
  }








  ##Ok Now we have the data corresponding to the phase 12 portion of the trial.
  ###Get Accrual in Phase III trial
  ACC = cumsum(rexp(Nmax*2,Accrue))
  ##This Will store our overall treament in the trial
  Trt = rep(NA,Nmax*2)
  ##Generate the outcomes for the control, we generate Nmax here to have extra if we switch doses
  TimeCont=rep(NA,Nmax)
  YECont=TimeCont
  YTCont=TimeCont


    if(Family=="Gamma"){
      for(b in 1:Nmax){
        YE=rbinom(1,1,ProbC[1])
        YT=rbinom(1,1,ProbC[2])


    TimeCont[b]=rgamma(1,alpha,1/exp(YE*betaC[1]+YT*betaC[2]+betaC[3]))
      }

    }

  if(Family=="Exponential"){
    for(b in 1:Nmax){
      YE=rbinom(1,1,ProbC[1])
      YT=rbinom(1,1,ProbC[2])


      TimeCont[b]=rgamma(1,1,1/exp(YE*betaC[1]+YT*betaC[2]+betaC[3]))
    }

  }


  if(Family=="Weibull"){
    for(b in 1:Nmax){
      YE=rbinom(1,1,ProbC[1])
      YT=rbinom(1,1,ProbC[2])


      TimeCont[b]=rweibull(1,alpha,exp(YE*betaC[1]+YT*betaC[2]+betaC[3]))
    }

  }


  if(Family=="Lognormal"){
    for(b in 1:Nmax){
      YE=rbinom(1,1,ProbC[1])
      YT=rbinom(1,1,ProbC[2])


      TimeCont[b]=rlnorm(1,YE*betaC[1]+YT*betaC[2]+betaC[3],alpha)
    }

  }



  ##Now we have the control outcomes

  INDSTOP=0

  ##Set the trial time to 0
  trial.time=0
  trial.time1=0
  ##Setup Storage matrices for this simulation rep
  Times=rep(NA,Nmax*2)
  TimesREG=Times
  YE=Times
  YT=Times
  Doses=Times
  Treat=Times

























  ##Phase 12 Info
  XOLD = Phase12[,1]
  YEOLD=Phase12[,2]
  YTOLD=Phase12[,3]
  ACCOLD =Phase12[,4]
    TIMEOLD=TIMES
OptDose=Opt


##Start Phase 3
  trial.time=ACC[1]
  ##Setup Storage matrices
  Times=rep(NA,Nmax1*2)
  YE=Times
  YT=Times
  Doses=Times
  Treat=Times
  ##Stores Best treatment
  Best=NA

  ###Enroll patients in the trial between two doses until the first interim look for dose switching
  NDeath=0
  i=1
  Nmax=Nmax1

  trial.time1=trial.time
Time2=Time12




  ###Simulate Trial until decision to switch doses
  while(NDeath<NLookSwitch){


    trial.time = trial.time +(ACC[i+1]-ACC[i])

    Trt[i]=rbinom(1,1,.5)


    ##If Group=0 it's in the control
    if(Trt[i]==0){
      Times[i]=TimeCont[i]
      Treat[i]=0
      Doses[i]=0
    }else{
      ##Generate Binary indicators for first part of trial
      YE[i]  = rbinom(1,1,probE[OptDose])
      YT[i] =  rbinom(1,1,probT[OptDose])

      if(Family=="Gamma"){
      Times[i]= rgamma(1,alpha,1/exp(betaA[1]*Dose[OptDose]+betaA[4]*Dose[OptDose]^2+YE[i]*betaA[2]+YT[i]*betaA[3]+betaA[5]))
      }

      if(Family=="Exponential"){
        Times[i]= rgamma(1,1,1/exp(betaA[1]*Dose[OptDose]+betaA[4]*Dose[OptDose]^2+YE[i]*betaA[2]+YT[i]*betaA[3]+betaA[5]))
      }

      if(Family=="Weibull"){
        Times[i]= rweibull(1,alpha,exp(betaA[1]*Dose[OptDose]+betaA[4]*Dose[OptDose]^2+YE[i]*betaA[2]+YT[i]*betaA[3]+betaA[5]))
      }


      if(Family=="Lognormal"){
        Times[i]= rlnorm(1,(betaA[1]*Dose[OptDose]+betaA[4]*Dose[OptDose]^2+YE[i]*betaA[2]+YT[i]*betaA[3]+betaA[5]),alpha)
      }


      Doses[i]=OptDose

    }





    ##Calculate the numnber of toxicities in the trial thusfar
    NDeath = sum((ACC+Times)<trial.time,na.rm=TRUE)


    ##Add to i
    i=i+1


  }
  ##Now we have NLook patient events. Let's see what dose we should use for the remainder of the trial
  TimesREG = Times
  YEREG=YE
  YTREG=YT
  DosesREG=Doses
  iREG=i
  trial.reg = trial.time
  TrtREG=Trt
  Y = Times[Trt>0]
  ACC1=ACC[Trt>0]
  ACC1=ACC1[!is.na(Y)]
  Y=Y[!is.na(Y)]
  Y = pmin(Y,trial.time-ACC1)
  I = (Y+ACC1)<trial.time
  YE1 = subset(YE,!(is.na(YE)))
  YT1 = subset(YT, !(is.na(YT)))
  YE1=c(YE1,YEOLD)
  YT1=c(YT1,YTOLD)
  ##Package up the phase I/II data too
  Doses1= subset(Doses, !(is.na(Doses)))
  Doses1=Doses1[Doses1>0]
  ##Time accrued into Phase I/II portion of trial.
  ACCOLD1=max(ACCOLD)-ACCOLD - Time2 + Twait
  TIMEOLD1 = pmin(TIMEOLD,ACCOLD1+trial.time)
  IOLD1=(TIMEOLD==TIMEOLD1)
  Y=c(Y,TIMEOLD1)
  I=c(I,IOLD1)
  Doses2=c(Doses1,XOLD)


OptDose3=OptDose
OptDose=  Reoptimize1(Y,I,YE1,YT1, Doses2, Dose, Hypermeans,  Hypervars, B )
OptDose1=OptDose



NumLooks = length(NLook)

##Now we have the newly reoptimized dose in terms of mean survival. Let's finish the trial, and report the results of both
##This and the conventional approach.
DECISION=NA
DECISIONREG=NA

Decision=NA
DecisionREG=NA
DECISON=NA
DECISIONREG=NA
if(OptDose==OptDose3){
  ##We just need to run one trial, all results are the same
  NDeath=0
  i=iREG

  INDSTOP=0

  for(k in 1:NumLooks){

    if(INDSTOP==0){

      NDeath=0

      while(NDeath<NLook[k]){

        if(i==Nmax){
          INDSTOP=1
          break
        }
        trial.time = trial.time +(ACC[i+1]-ACC[i])

        Trt[i]=rbinom(1,1,.5)


        ##If Group=0 it's in the control
        if(Trt[i]==0){
          Times[i]=TimeCont[i]
          Treat[i]=0
          Doses[i]=0
        }else{
          ##Generate Binary indicators for first part of trial
          YE[i]  = rbinom(1,1,probE[OptDose])
          YT[i] =  rbinom(1,1,probT[OptDose])


          if(Family=="Gamma"){
            Times[i]= rgamma(1,alpha,1/exp(betaA[1]*Dose[OptDose]+betaA[4]*Dose[OptDose]^2+YE[i]*betaA[2]+YT[i]*betaA[3]+betaA[5]))
          }

          if(Family=="Exponential"){
            Times[i]= rgamma(1,1,1/exp(betaA[1]*Dose[OptDose]+betaA[4]*Dose[OptDose]^2+YE[i]*betaA[2]+YT[i]*betaA[3]+betaA[5]))
          }

          if(Family=="Weibull"){
            Times[i]= rweibull(1,alpha,exp(betaA[1]*Dose[OptDose]+betaA[4]*Dose[OptDose]^2+YE[i]*betaA[2]+YT[i]*betaA[3]+betaA[5]))
          }


          if(Family=="Lognormal"){
            Times[i]= rlnorm(1,(betaA[1]*Dose[OptDose]+betaA[4]*Dose[OptDose]^2+YE[i]*betaA[2]+YT[i]*betaA[3]+betaA[5]),alpha)
          }





          Doses[i]=OptDose

        }





        ##Calculate the numnber of toxicities in the trial thusfar
        NDeath = sum((ACC+Times)<trial.time,na.rm=TRUE)


        ##Add to i
        i=i+1


      }


      if(INDSTOP==1){


          trial.time = quantile((ACC+Times), NLook[k]/Nmax,na.rm=TRUE)



          ##Calculate the numnber of toxicities in the trial thusfar
          NDeath = sum((ACC+Times)<trial.time,na.rm=TRUE)






      }

    }else{



        trial.time = quantile((ACC+Times), NLook[k]/Nmax,na.rm=TRUE)



        ##Calculate the numnber of toxicities in the trial thusfar
        NDeath = sum((ACC+Times)<trial.time,na.rm=TRUE)







    }




    ##Lets remove all Data that's NOT part of the pruned dose
    Y = Times[!is.na(Times)]
    Trt1 = Trt[!is.na(Times)]
    ACC1=ACC[!is.na(Times)]
    Y = pmin(Y,trial.time-ACC1)
    I = (Y+ACC1)<trial.time





    LR= coxph(Surv(Y,I)~Trt1)


    if(sqrt(LR$wald.test)<Fut[k]){
      ##Stop Due to Nutility
      DECISION=0
      Npats=i
      TRIALTIMES=trial.time
      LOOK=k
      LOOKREG=k



      DECISIONREG=0
      NpatsREG=i
      TRIALTIMESREG=trial.time

      break

    }else{

      if(sqrt(LR$wald.test)>Sup[k]){


        TRIALTIMESREG=trial.time
        NpatsREG=i


        TRIALTIMES=trial.time
        Npats=i


        if(LR$coefficients[1]<0){
          DECISIONREG=1
          DECISION=1

        }else{
          DECISION=-1
          DECISIONREG=-1

        }


        break

      }

    }


  }



  if(is.na(DECISION)){

    ##Stop Due to Nutility
    DECISION=0
    Npats=i
    TRIALTIMES=trial.time

    DECISIONREG=0
    NpatsREG=i
    TRIALTIMESREG=trial.time




  }








  }else{

  ##We need two separate trials, one for each control and regular
  ##Do Regular Trial First

  NDeath=0
  i=iREG

  INDSTOP=0
  DECISION=NA
  DECISIONREG=NA

  for(k in 1:NumLooks){

    OptDose=OptDose3

    if(INDSTOP==0){

      while(NDeath<NLook[k]){
        if(i==Nmax){
          INDSTOP=1
          break
        }

        trial.reg = trial.reg +(ACC[i+1]-ACC[i])

        TrtREG[i]=rbinom(1,1,.5)



        ##If Group=0 it's in the control
        if(TrtREG[i]==0){
          TimesREG[i]=TimeCont[i]
          Treat[i]=0
          Doses[i]=0
        }else{

          ##Generate Binary indicators for first part of trial
          YE[i]  = rbinom(1,1,probE[OptDose])
          YT[i] =  rbinom(1,1,probT[OptDose])


          if(Family=="Gamma"){
            TimesREG[i]= rgamma(1,alpha,1/exp(betaA[1]*Dose[OptDose]+betaA[4]*Dose[OptDose]^2+YE[i]*betaA[2]+YT[i]*betaA[3]+betaA[5]))
          }

          if(Family=="Exponential"){
            TimesREG[i]= rgamma(1,1,1/exp(betaA[1]*Dose[OptDose]+betaA[4]*Dose[OptDose]^2+YE[i]*betaA[2]+YT[i]*betaA[3]+betaA[5]))
          }

          if(Family=="Weibull"){
            TimesREG[i]= rweibull(1,alpha,exp(betaA[1]*Dose[OptDose]+betaA[4]*Dose[OptDose]^2+YE[i]*betaA[2]+YT[i]*betaA[3]+betaA[5]))
          }


          if(Family=="Lognormal"){
            TimesREG[i]= rlnorm(1,(betaA[1]*Dose[OptDose]+betaA[4]*Dose[OptDose]^2+YE[i]*betaA[2]+YT[i]*betaA[3]+betaA[5]),alpha)
          }







        }





        ##Calculate the numnber of toxicities in the trial thusfar
        NDeath = sum((ACC+TimesREG)<trial.reg,na.rm=TRUE)


        ##Add to i
        i=i+1


      }


      if(INDSTOP==1){






          trial.time = quantile(ACC+TimesREG,NLook[k]/Nmax,na.rm=TRUE)


          ##Calculate the numnber of toxicities in the trial thusfar
          NDeath = sum((ACC+TimesREG)<trial.time,na.rm=TRUE)





      }



    }else{



        trial.reg = quantile((ACC+TimesREG),NLook[k]/Nmax,na.rm=TRUE)








        ##Calculate the numnber of toxicities in the trial thusfar
        NDeath = sum((ACC+TimesREG)<trial.reg,na.rm=TRUE)








    }

    ##Lets remove all Data that's NOT part of the pruned dose
    Y = TimesREG[1:(i-1)]
    Trt1 = TrtREG[1:(i-1)]
    ACC1=ACC[1:(i-1)]
    Y = pmin(Y,trial.reg-ACC1)
    I = (Y+ACC1)<trial.time





    LR= coxph(Surv(Y,I)~Trt1)


    if(sqrt(LR$wald.test)<Fut[k]){
      ##Stop Due to Nutility

      DECISIONREG=0
      NpatsREG=i
      TRIALTIMESREG=trial.reg




      break
    }else{

      if(sqrt(LR$wald.test)>Sup[k]){


        TRIALTIMESREG=trial.reg
        NpatsREG=i



        if(LR$coefficients[1]<0){
          DECISIONREG=1

        }else{
          DECISIONREG=-1

        }


        break

      }

    }


  }



  if(is.na(DECISIONREG)){

    DECISIONREG=0
    NpatsREG=i
    TRIALTIMESREG=trial.reg

  }




  i=iREG


  OptDose=OptDose1

  INDSTOP=0


  for(k in 1:NumLooks){


    NDeath=0
    if(INDSTOP==0){
      while(NDeath<NLook[k]){
        if(i==Nmax){
          INDSTOP=1
          break
        }

        trial.time = trial.time +(ACC[i+1]-ACC[i])

        Trt[i]=rbinom(1,1,.5)


        ##If Group=0 it's in the control
        if(Trt[i]==0){
          Times[i]=TimeCont[i]
          Treat[i]=0
          Doses[i]=0
        }else{
          ##Generate Binary indicators for first part of trial
          YE[i]  = rbinom(1,1,probE[OptDose])
          YT[i] =  rbinom(1,1,probT[OptDose])


          if(Family=="Gamma"){
            Times[i]= rgamma(1,alpha,1/exp(betaA[1]*Dose[OptDose]+betaA[4]*Dose[OptDose]^2+YE[i]*betaA[2]+YT[i]*betaA[3]+betaA[5]))
          }

          if(Family=="Exponential"){
            Times[i]= rgamma(1,1,1/exp(betaA[1]*Dose[OptDose]+betaA[4]*Dose[OptDose]^2+YE[i]*betaA[2]+YT[i]*betaA[3]+betaA[5]))
          }

          if(Family=="Weibull"){
            Times[i]= rweibull(1,alpha,exp(betaA[1]*Dose[OptDose]+betaA[4]*Dose[OptDose]^2+YE[i]*betaA[2]+YT[i]*betaA[3]+betaA[5]))
          }


          if(Family=="Lognormal"){
            Times[i]= rlnorm(1,(betaA[1]*Dose[OptDose]+betaA[4]*Dose[OptDose]^2+YE[i]*betaA[2]+YT[i]*betaA[3]+betaA[5]),alpha)
          }





          Doses[i]=OptDose

        }





        ##Calculate the numnber of toxicities in the trial thusfar
        NDeath = sum((ACC+Times)<trial.time,na.rm=TRUE)


        ##Add to i
        i=i+1


      }


      if(INDSTOP==1){



          trial.time = quantile((ACC+Times),NLook[k]/Nmax,na.rm=TRUE)






          ##Calculate the numnber of toxicities in the trial thusfar
          NDeath = sum((ACC+Times)<trial.time,na.rm=TRUE)







      }

    }else{

        trial.time = quantile((ACC+Times),NLook[k]/Nmax,na.rm=TRUE)






        ##Calculate the numnber of toxicities in the trial thusfar
        NDeath = sum((ACC+Times)<trial.time,na.rm=TRUE)






    }

    ##Lets remove all Data that's NOT part of the pruned dose

    ##Lets remove all Data that's NOT part of the pruned dose
    Y = Times[1:(i-1)]
    Trt1 = Trt[1:(i-1)]
    Doses1 = Doses[1:(i-1)]

    ACC1=ACC[1:(i-1)]
    Y = pmin(Y,trial.time-ACC1)
    I = (Y+ACC1)<trial.time


    Trt2=subset(Trt1, Trt1==0 | Doses1==OptDose)
    Y=subset(Y, Trt1==0 | Doses1==OptDose)
    I=subset(I, Trt1==0 | Doses1==OptDose)





    LR= coxph(Surv(Y,I)~Trt2)


    if(sqrt(LR$wald.test)<Fut[k]){
      ##Stop Due to Nutility

      DECISION=0
      Npats=i
      TRIALTIMES=trial.time



      break

    }else{

      if(sqrt(LR$wald.test)>Sup[k]){


        TRIALTIMES=trial.time
        Npats=i


        TRIALTIMESSET=trial.time
        NpatsSET=i


        if(LR$coefficients[1]<0){
          DECISION=1

        }else{
          DECISION=-1
        }


        break

      }

    }


  }

  if(is.na(DECISION)){

    DECISION=0
    Npats=i
    TRIALTIMES=trial.time



  }


}



Phase123 = rep(0,4)
Phase123[1]=OptDose1
Phase123[2]=DECISION
Phase123[3]=Npats
Phase123[4]=TRIALTIMES


G = rep(0,4)
G[1]=OptDose3
G[2]=DECISIONREG
G[3]=NpatsREG
G[4]=TRIALTIMESREG




Z=as.list(c(0,0))
Z[[1]]=Phase123
Z[[2]]=G

return(Z)


}





