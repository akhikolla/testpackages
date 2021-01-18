#' Simulates replications from EFF-TOX phase 12 trial.
#'
#' Simulates trial replications from the EFF-TOX phase 12 trial trial using either adaptive randomization or fixed dose assignment.   Prints the true utility scores, dose selection probability, 
#' average number of patients treated at each dose, average number of responses, average number of toxicities, and Delta value. Returns trial outcomes.
#' @param NSims Number of simulations.
#' @param Dose Log-standardized doses (log(Raw Dose)-mean(log(Raw Dose))).
#' @param Contour Vector containing 4 entries used to make the desireability function. Contour(1) contains a desired toxicity probability given efficacy, Countour(2) contains a desired efficacy probability given toxicity, and (Contour(3),Contour(4)) is an equally desireable pair of efficacy and toxicity probabilities that are non zero or one.
#' @param PE True Efficacy Probability for each dose.
#' @param PT True toxicity probaiblity for each dose.
#' @param corET Correlation parameter between Efficacy and Toxicity status.
#' @param Nmax Maximum Sample size.
#' @param cohort Cohort Size.
#' @param NF Number of fixed assignment patients until adaptive randomization. If NF equals Nmax, the trial is conducted without adaptive randomization.
#' @param CutE Cutoff for efficacy acceptability.
#' @param CutT Cutoff for toxicity acceptability.
#' @param AcceptE Probability threshold for efficacy acceptability.
#' @param AcceptT Probability threshold for toxicity acceptability.
#' @param HypermeansEFF Vector containing prior hypermeans of length 6 for EFF-TOX parameters
#' @param HypervarsEFF Vector containing prior hypervariances of length 6 for EFF-TOX parameters
#' @return A list of size NSims with results from each simulated trial. Each entry contains a list with (1) the optimal dose selected, (2) the posterior mean utility for each dose level, (3) a matrix containing the dose given, the efficacy outcome and the toxicity outcome for each patient.
#' @references Thall, P.F. and Cook, J.D. (2004). Dose-finding based on efficacy-toxicity trade-offs. Biometrics 60, 684-693.
#' @references Chapple AG, Thall PF. A Hybrid Phase 123 Clinical Trial Design Allowing Dose Re-optimization in Phase III. Biometrics. Epub ahead of print 26 October 2018.
#' @importFrom Phase123 RunAdaptiveEffToxTrial
#' @examples
#' ##Trial PArameters here
#' Nmax=30 ##Number of patients to enroll
#' NF=30 ##Number until AR if NF=Nmax, there's no AR.
#' cohort=3
#' #' Raw Dose Values
#' Dose = c(1,2,3,3.5,5)
#' Dose=log(Dose)-mean(log(Dose))
#' ## Contour Vector
#' Contour = c(.35, .75,.7,.4)
#' #Starting Dose
#' DoseStart=1
#' ##Safety Parameters
#' CutE=.3
#' CutT=.4
#' AcceptE=.1
#' AcceptT=.1
#' ##Hypermeans
#' HypermeansEFF = c(.022,3.45,0,-4.23,3.1,0)
#' ##Hypervariances 
#' HypervarsEFF = c(2.6761, 2.6852, .2, 3.1304, 3.1165, 1)
#' HypervarsEFF=HypervarsEFF^2 
#' #True Efficacy and Toxicity probabilities
#' PE=c(.2,.4,.6,.7,.7)
#' PT=c(.2,.2,.2,.3,.5)
#' corET=0
#' ##Number of simulations
#' NSims=2
#' RESULTS=RunAdaptiveEFFTOX(NSims,Dose,PE, PT, corET, Nmax, cohort, 
#' NF, Contour, CutE, CutT, AcceptE, AcceptT, HypermeansEFF, HypervarsEFF )
#'@export
RunAdaptiveEFFTOX=function(
NSims, ##Number of simulations
Dose, ##logdoses
PE, ##True Efficacy Probability for each dose
PT, ##True toxicity probaiblity for each dose
corET, ##Correlation parameter between Eff and Tox
Nmax, ##MAximum Sample size
cohort, ##Cohort Size
NF, ##Number of fixed assignment patients until adaptive randomization
Contour, ##Contour vector
CutE, ##Cutoff For efficacy acceptability
CutT, ##Cutoff for toxicity acceptability
AcceptE, ##Probability threshold for eff acceptability
AcceptT, ##Probability threshold for tox acceptability
HypermeansEFF, ##Hypermeans for 
HypervarsEFF ##Hypervariances
){
  
  Hypermeans = HypermeansEFF
  Hypervars=HypervarsEFF
  
  
  
  
  nDose=length(Dose)
  
  
  
  TRUEUT=rep(NA,nDose)
  
  
  
  
  PiLim=c(CutE,CutT)
  ProbLim=c(AcceptE,AcceptT)
  B=2000
  
  
  
  
  NumTrt=matrix(rep(0,nDose*NSims),nrow=NSims)
  
  ##Setup Simulation parameters
  SIMSTORE = as.list(rep(0,NSims))
  
  DoseOpt=rep(0,NSims)
  NTox = rep(0,NSims)
  NEff=rep(0,NSims)
  
  DoseStore=rep(0,Nmax)
  YEStore = rep(0,Nmax)
  YTStore =rep(0,Nmax)
  
  
  
  ##No correlation
  for(m1 in 1:NSims){
    
    
    List=  RunAdaptiveEffToxTrial(1,Dose, Hypermeans,  Hypervars,
                                  Contour, PiLim, ProbLim,  cohort, Nmax,  NF, B, 1, PE, PT )
    
    
    
    Doses=List[[3]][,1]+1
    YE=List[[3]][,2]
    YT=List[[3]][,3]
    OptDose=List[[1]]
    
    
    ##This simulation doesn't count now.
    
    
    NTox[m1]=sum(YT)
    NEff[m1]=sum(YE)
    DoseOpt[m1]=OptDose
    
    for(j in 1:nDose){
      NumTrt[m1,j]=sum(Doses==j)
    }
    
    
    
  }
  
  
  
  
  ##Calculate True Accept
  TRUECEPT = (PE>=CutE)*(PT<=CutT)
  ##These are the doses that are truly acceptable
  
  ##Calculate the true desireability
  
  desire=rep(NA,nDose)
  
  
  for(k in 1:nDose){
    desire[k]=GetDesire(PE[k],PT[k],Contour)
    
  }
  
  
  cat("True Mean desireability Scores
      
      ")  
  
  print(desire)
  
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
  
  U1=(desire-mean(desire))/sd(desire)
  U1=c(0,U1*TRUECEPT)
  cat("Delta phi
      
      ")
  
  print(sum(prob1*(max(U1)-U1)))
  
  
  return(SIMSTORE)
  
  
  
  
  
}
