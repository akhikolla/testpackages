#' Gives operating characteristics of phase 123 and conventional design.
#'
#' Returns the probability of selecting the optimal dose, type I error, generalized power, probability of making the best decision, average number of patients treated and average trial duration.
#' @param Results List containing phase 123 and conventional design results.
#' @param Means True mean survival times for experimental agents at each dose.
#' @param CMu Mean survival time for the control therapy.
#' @param Delta Desired improvement in survival.
#' @param Hyp Null=0 or alternative=1 hypthesis
#' @references
#' [1] Chapple and Thall (2018). A Hybrid Phase 12/3 Clinical Trial Design Allowing Dose-Re-Optimization in Phase 3 Biometrics. Under Review.
#' @examples
#' ##True Mean Control
#' CMu=24
#' ##True Means Agent
#' Means = c(27,32,38,42,28)
#' ##Desired improvement in mean survival
#' Delta=12
#' ##Random Trial results
#'  Results=as.list(c(0,0))
#'  nSims=5
#'  X=matrix(rep(NA,nSims*4),nrow=nSims)
#'  ##DoseSelected
#'  X[,1]=c(2,3,4,4,3)
#'  X[,2]=c(0,1,1,1,1)
#'  X[,3]=c(270,500,500,420,400)
#'  X[,4]=c(70,85,88,70,88)
#'  Results[[1]]=X
#'  X[,1]=c(2,3,5,4,2)
#'  X[,2]=c(0,1,0,1,0)
#'  X[,3]=c(270,500,450,420,415)
#'  X[,4]=c(70,82,80,70,79)
#'Results[[2]]=X
#'Hyp=1
#'ReturnOCS(Results,Means,CMu,Delta,Hyp)
#'@export
ReturnOCS=function(Results,Means,CMu,Delta,Hyp){
PH123=Results[[1]]
PH12=Results[[2]]
  HYP=Hyp
##Best Dose
Best=which(Means==max(Means))
ALL = which(Means>(CMu+Delta))
TRUEMEAN=Means
  CHOSENDOSE = PH123[,1]
  StartingDose=PH12[,1]
  DECISION=PH123[,2]
  DECISIONREG=PH12[,2]
  Npats = PH123[,3]
  NpatsREG = PH12[,3]
  TRIALTIMES = PH123[,4]
  TRIALTIMESREG = PH12[,4]
  cat("TRIAL TIMES
")
  cat("Conventional Design")
  print(mean(TRIALTIMESREG/12, na.rm=TRUE))
  cat("Phase123 Design")
  print(mean(TRIALTIMES/12, na.rm=TRUE))
  cat("Number of Patients
")
  cat("Conventional Design")
  print(mean(NpatsREG, na.rm=TRUE))
  cat("Phase123 Design")
  print(mean(Npats, na.rm=TRUE))
  cat("Percentage of times best dose selected
")
  cat("Conventional Design")
  print(mean(StartingDose==Best, na.rm=TRUE))
  cat("Phase123 Design")
  print(mean(CHOSENDOSE==Best, na.rm=TRUE))



  GETDELTA=function(DECISION,CHOSEN,TRUEMEAN,CMU){

    MEAN=rep(0,length(DECISION))


    for(b in 1:length(MEAN)){
      if(DECISION[b]==1){
        MEAN[b]=TRUEMEAN[CHOSEN[b]]-CMU
      }



    }

    return(mean(MEAN))

  }



  if(Hyp==1){

if(length(ALL)>1){


  cat("Generalized Power
")
  cat("Conventional Design ")
  print(mean(DECISIONREG==1 & StartingDose %in% ALL, na.rm=TRUE))
  cat("Phase123 Design ")
  print(mean(DECISION==1 & CHOSENDOSE %in% ALL, na.rm=TRUE))



  cat("Probability of best decisions
")
  cat("Conventional Design ")
  print(mean(DECISIONREG==1 & StartingDose==Best, na.rm=TRUE))
  cat("Phase123 Design ")
  print(mean(DECISION==1 & CHOSENDOSE==Best, na.rm=TRUE))


}else{

  cat("Generalized Power
")
  cat("Conventional Design ")
  print(mean(DECISIONREG==1 & StartingDose==Best, na.rm=TRUE))
  cat("Phase123 Design ")
  print(mean(DECISION==1 & CHOSENDOSE==Best, na.rm=TRUE))

}




    cat("bar{W} values: Average true improvement in survival
")


    cat("Conventional Design ")
    print(GETDELTA(DECISIONREG,StartingDose,TRUEMEAN,CMu))



    cat("Phase123 Design ")
    print(GETDELTA(DECISION,CHOSENDOSE,TRUEMEAN,CMu))


  }else{

    cat("Type I error
")
    cat("Conventional Design ")
    print(mean((DECISIONREG==1 | DECISIONREG==-1 )  & StartingDose==Best, na.rm=TRUE) + mean(DECISIONREG==1 & StartingDose !=Best,na.rm=TRUE))
    cat("Phase123 Design ")
    print( mean( DECISION==1  & (CHOSENDOSE==Best & StartingDose !=Best) , na.rm=TRUE)+
   mean((DECISION==1 | DECISION==-1)  & (StartingDose==Best & CHOSENDOSE==Best), na.rm=TRUE)+
    mean(DECISION==1 & CHOSENDOSE !=Best,na.rm=TRUE))


    }

}

