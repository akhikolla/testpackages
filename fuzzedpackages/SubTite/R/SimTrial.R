#' Simulates a Sub-TITE trial design
#'
#' Simulates replicates from a Sub-TITE trial with user specified true toxicity time distributions for different doses and subgroups and returns average summary statistics of the trial.
#' @importFrom stats pexp pweibull pgamma plnorm rexp rbinom runif rweibull rgamma rlnorm rnorm var nls
#' @importFrom Rcpp evalCpp
#' @param nSims Number of Trials to Simulate.
#' @param Nmax Maximum Number of Patients to enroll in the trial.
#' @param T1 Reference time for toxicity.
#' @param Target Target cumulative toxicity probability (or subgroup specific vector) at time T1.
#' @param Dose Standardized vector of doses to try.
#' @param DoseStart Dose (or vector of Doses) to enroll the first patient in each subgroup at.
#' @param Accrue Expected montly patient accrual rate.
#' @param groupprob Probability vector of subgroup assignment.
#' @param Upper Cutoff values used to determine if accrual in a subgroup should be suspended.
#' @param meanmu Prior mean of the baseline intercept parameter.
#' @param meanslope Prior mean of the baseline slope parameter.
#' @param MeanInts G-1 length vector of subgroup specific prior intercept means.
#' @param MeanSlopes G-1 length vector of subgroup specific prior slope means.
#' @param Family What distribution Family to simulate from. Options include: Exponential,Gamma, Lognormal, Uniform, Weibull.
#' @param SimTruth List of 2 #Groups by #Doses matrices containing the true parameter values needed for simulating from different true time to toxicity distributions.
#' @param VarInt Prior Variance of Intercept Parameters.
#' @param VarSlope Prior Variance of Slope Parameters.
#' @param phetero Prior prob of clustering
#' @param NSep Number of patients to assign based on no borrowing.
#' @param NBorrow Number of patients to assign based on no clustering
#' @param cohort Number of patients to enroll before escalating.
#' @param FULL Do we have to fully evaluate a cohort before escalating?
#' @return A list with first entry corresponding to summaries of the operating characteristics of the design including
#' @examples
#' ##Note: nSims  should be set larger than the example below.
#' nSims=1
#' ###TRIAL PARAMETERS###
#' ##Specify reference toxicity time and target
#' T1=6
#' Target=.3
#' ##Number of Groups
#' ##Specify upper bound for determining if the lowest dose is too toxic in a subgroup
#' Upper=c(.95,.95)
#' #' ##Standardized Dose Values and starting dose index
#' Dose=sort(rnorm(5))
#' DoseStart=1
#' ##Maximum Sample Size
#' Nmax=25
#' ##Number of patients to run separately
#' NSep=0
#' ##Number of patients to borrow, but NOT cluster
#' NBorrow=0
#' ##Number of patients to fully evaluate or TREAT before ESCALATING
#' cohort=3
#' ##Do we fully evaluate a cohort before escalating?
#' FULL=0
#' #HYPERPARAMETERS#
#' ##Hypermeans for baseline terms
#' meanmu=2.21
#' meanslope=-.57
#' ##Hypervectors for subgroup specific terms
#' MeanInts = c(.46)
#' MeanSlopes = c(.04)
#' ##Hypervariances
#' VarInt=5
#' VarSlope=1
#' ######SIMULATION TRUTH####
#' ##True Accrual Rate
#' Accrue=2
#' ##True Distribution of subgroups
#' groupprob=c(.5,.5)
#' ##True Group Toxicity probabilities at each dose level
#' GroupProb =matrix(c(.05,.3,.6,.7,.8,.01,.02,.13,.27,.5),nrow=2,byrow=TRUE)
#' ##True Simulation distribution
#' Family="Uniform"
#' SimTruth = as.list(c(0,0))
#' SimTruth[[1]]=GroupProb
#' SimTruth[[2]]=GroupProb
#' phetero=.9
#' RESULTS=SimTrial(nSims,Nmax,T1,Target,Dose,DoseStart,
#'               Upper,Accrue,groupprob,meanmu,meanslope,
#'               MeanInts,MeanSlopes,VarInt,VarSlope,phetero,
#'               Family,SimTruth,NSep,NBorrow,cohort,FULL)
#'               RESULTS[[1]]
#' @export
SimTrial = function(nSims,Nmax,T1,Target,Dose, DoseStart,Upper,Accrue,groupprob,meanmu,meanslope,MeanInts,MeanSlopes,
                    VarInt,VarSlope,phetero,Family,SimTruth,NSep,NBorrow,cohort,FULL){

pmono=phetero

  Param1=SimTruth[[1]]
  Param2=SimTruth[[2]]


  nGroups=nrow(Param1)
  nGroup=nGroups
  nDose=ncol(Param1)
  upper=Upper

  Dist = matrix(ncol=nGroups,nrow=nSims)


  target=Target
  B=2000




  nDoses=length(Dose)
  target=Target



  ##If length(Target)==1 or length(DoseStart)==1, let's make a vector containing the group specific targets

  if(length(Target)==1){
    Target1=Target
    Target=rep(NA,nrow(Param1))
    for(k in 1:nrow(Param1)){
      Target[k]=Target1
    }

  }


  if(length(DoseStart)==1){
    DoseStart1=DoseStart
    DoseStart=rep(NA,nrow(Param1))
    for(k in 1:nrow(Param1)){
      DoseStart[k]=DoseStart1
    }

  }


  if(length(Upper)==1){
    Upper1=Upper
    Upper=rep(NA,nrow(Param1))
    for(k in 1:nrow(Param1)){
      Upper[k]=Upper1
    }

  }


  ##Now check for errors...


  ERRHOLD=c(length(Target), nrow(Param1), nrow(Param2), length(Upper), length(MeanInts)+1, length(MeanSlopes)+1)

  HOLD=0
  ##Check for errors in dimension specification
  for(k in 1:length(ERRHOLD)){
    for(m in 1:length(ERRHOLD)){
      if(ERRHOLD[k] != ERRHOLD[m]){
        HOLD=1
      }
    }
  }

  if(HOLD==1){
    message("Dose tried matrix,target toxicity vector, toxicity threshold, or subgroup hyperparameter vector has incorrect dimensions")
  }



  ##Repackage MeanInts and MeanSlopes
  MeanInts=c(0,MeanInts)
  MeanSlopes=c(0,MeanSlopes)





  if(Family=="Exponential"){
    GroupProb = pexp(T1,1/Param1)
    Fam=0
  }

  if(Family=="Uniform"){

    GroupProb=Param1
    Fam=3
  }

  if(Family=="Weibull"){
    GroupProb=pweibull(T1,Param1,Param2)
    Fam=4
  }

  if(Family=="Gamma"){
    GroupProb=pgamma(T1,Param1,Param2)
    Fam=1
  }

  if(Family=="Lognormal"){
    GroupProb=plnorm(T1,Param1,Param2)
    Fam=2
  }




  POPT=rep(NA,nGroups)
  WHICHOPT=POPT
  ##What is POPT

  for(m in 1:nGroups){
    min1 = min(abs(GroupProb[m,]-Target[m]))
    which1 =  which(abs(GroupProb[m,]-Target[m])==min1)
    POPT[m]=GroupProb[m,which1]
    WHICHOPT[m]=which1

  }



  ##Feed Everything into the trial
  Results= SimTrial1( nSims, Nmax,  T1, Target,  Dose,  DoseStart,Upper, Accrue, groupprob,
                      Fam, Param1,Param2,meanmu, meanslope,MeanInts, MeanSlopes,VarInt,VarSlope,pmono, NSep,NBorrow, cohort, FULL )

  DoseOpt = Results[[1]]
  NTox=Results[[2]]


  ##Calculate \Delta_1,\Delta_2
  for(b in 1:nSims){

    for(m in 1:nGroups){
      if(DoseOpt[b,m]==0){
        Dist[b,m]=POPT[m]
      }else{
        Dist[b,m]=abs(GroupProb[m,DoseOpt[b,m]]-POPT[m])
      }
    }







  }





  ##Obtain the prob of best

  OptProb=WHICHOPT
  ###GET Psel and PM
  PM=OptProb
  for(m in 1:nGroups){
    OptProb[m]=mean(DoseOpt[,m]==WHICHOPT[m])
    if(WHICHOPT[m]==1){
      PM[m]=mean(DoseOpt[,m] %in% c(1,2))
    }else{
      if(WHICHOPT[m]==ncol(GroupProb)){
        PM[m]=mean(DoseOpt[,m] %in% c(ncol(GroupProb),ncol(GroupProb)-1))

      }else{
        ##Dose in the middle some where
        PM[m]=mean(DoseOpt[,m] %in% c(WHICHOPT[m]-1,WHICHOPT[m],WHICHOPT[m]+1))

      }


    }
    }

  for(m in 1:nGroups){
    if(GroupProb[m,1]>Target[m]){
      OptProb[m]=NA
      Dist[m]=NA
      PM[m]=NA
    }
  }









  Z=as.list(rep(0,8))
  GroupProb=as.data.frame(GroupProb)

  for(k in 1:length(MeanInts)){
    rownames(GroupProb)[k]=paste0("Subgroup ",k)
  }



  for(k in 1:ncol(GroupProb)){
    colnames(GroupProb)[k]=paste0(k)
  }


  Z[[1]]=GroupProb

  X=data.frame(matrix(ncol=length(Dose)+1,nrow=length(MeanInts)))

  for(k in 1:nrow(X)){
    rownames(X)[k]=paste0("Subgroup ",k)
  }

  for(k in 1:nrow(X)){
    for(j in 1:length(Dose)){
      X[k,j]=mean(DoseOpt[,k]==j)
    }

    j=length(Dose)+1
    X[k,j]=mean(DoseOpt[,k]==0)

  }




  for(k in 1:length(Dose)){
    colnames(X)[k]=paste0(k)
  }


  colnames(X)[length(Dose)+1]="P[Stop]"

  Z[[2]]=X



  NTREATED=Results[[5]]
  GROUPS=Results[[6]]
  GROUPS=GROUPS+1


  X=data.frame(matrix(ncol=length(Dose),nrow=length(MeanInts)))
  for(k in 1:length(MeanInts)){
    rownames(X)[k]=paste0("Subgroup ",k)
  }
  for(k in 1:nrow(X)){
    for(j in 1:length(Dose)){
      NTREATED1=NTREATED[GROUPS==k]
      X[k,j]=sum(NTREATED1==Dose[j])

    }
  }

  X=X/nSims



  for(k in 1:length(Dose)){
    colnames(X)[k]=paste0(k)
  }


  Z[[3]]=X
  Z[[4]]=colMeans(Dist)
  Z[[5]]=OptProb
  Z[[6]]=PM
  Z[[7]]=colMeans(NTox)
  Z[[8]]=mean(Results[[3]])
  names(Z)=c("True Toxicity Probabilities","Dose Selection Probabilities",
             "Average # Treated at each Dose","Delta","Psel", "+- MTD","Average # of Toxicities","Average Trial Duration")




  List1 = as.list(rep(0,6))
  List1[[1]]=DoseOpt
  List1[[2]]=Dist
  List1[[3]]=NTox
  List1[[4]]=Results[[3]]


  List1[[5]]=Results[[4]]
  List1[[6]]=Results[[5]]
  List1[[7]]=Results[[6]]

  names(List1)=c("Optimal Dose Selected","Delta Values", "Number of Toxicities","Trial Times","Toxicity Indicators",
                 "Doses Given","Subgroup Indicators")

  Z1 = as.list(c(0,0))

  Z1[[1]]=Z
  Z1[[2]]=List1
  names(Z1)=c("Simulation Summaries","Simulation Outputs")
  return(Z1)


}
