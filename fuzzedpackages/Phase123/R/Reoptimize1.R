#' Gives the Optimal Dose for enrolling next patient cohort. Used in the SimPhase123 function.
#'
#' This function returns the optimal dose number to assign the next patient cohort or stops the trial if no dose is deemed acceptable.
#' @param Y Vector containing observed patient survival or follow up times.
#' @param I Vector indicating whether each patient experienced an exent.
#' @param YE   Vector containing observed efficacy indicators.
#' @param YT   Vector containing observed toxicity indicators.
#' @param Doses Vector containing standardized doses of patients in trial.
#' @param Dose Vector containing the standardized doses considered.
#' @param Hypermeans Vector containing prior hypermeans of length 6 for Eff-Tox parameters.
#' @param Hypervars Vector containing prior hypervariances of length 6 for Eff-Tox parameters.
#' @param B Number of iterations to perform in the MCMC.
#' @importFrom  stats sd
#' @references
#' [1] Chapple and Thall (2018).A Hybrid Phase I-II/III Clinical Trial Design Allowing Dose Re-Optimization in Phase III. Biometrics. In Press,
#' @examples
#'##Doses, YE,YT
#'Doses= c(1,1,1,2,2,2,1,1,1,3,3,3,1,1,1,2,2,2)
#'YE = c(0,0,1,1,1,0,0,0,0,1,1,1,0,0,1,1,1,0)
#'YT=c(0,0,0,1,1,0,1,0,0,1,1,1,0,0,0,1,0,0)
#'Y=rexp(length(YE))
#'I=rbinom(length(YE),1,.9)
#'##Vector of Numerical Doses
#'Dose = c(1,2,3,3.5,5)
#'Dose=(Dose-mean(Dose))/sd(Dose)
#'##Hypermeans for Eff-Tox
#'Hypermeans = c(.022,3.45,0,-4.23,3.1,0)
#'Hypervars = c(2.6761, 2.6852, .2, 3.1304, 3.1165, 1)
#'Hypervars=Hypervars^2
#'###Number of iterations
#'B=20000
#'Reoptimize1(Y,I,YE,YT, Doses, Dose, Hypermeans,  Hypervars,B)
#' @export
Reoptimize1=function(Y,I,YE,YT, Doses, Dose, Hypermeans,  Hypervars, B ){



  YE1=YE
  YT1=YT

  ##Use EFFTOX PROGRAM TO GET probmat quantities!!
  G2=EFFTOX(YE1, YT1, Dose[Doses],Dose,  Hypermeans,  Hypervars, B )


  probmat1=G2




  MaxObs = matrix(rep(0,length(Dose)*4),nrow=4)



  YE1=(YE1-mean(YE1))/sd(YE1)
  YT1=(YT1-mean(YT1))/sd(YT1)


  ##Have Data For trial, Run MCMC


  MaxObs=MaxObs*0+max(Y)

  G1=PieceMCMC(Y,I,YE1,YT1,Dose[Doses],Dose,B,probmat1,MaxObs)



  G2=G1[[1]][(B/2):B,]

  VEC1 = colMeans(G2,na.rm=TRUE)

  for(k in 1:length(VEC1)){
    if(is.nan(VEC1[k])){
      VEC1[k]=0
    }




  }








  ##Change
  OptDose = which(VEC1==max(VEC1))




  ## Check if two cohorts treated
  for(m in 1:length(Dose)){
    if(sum(Doses==OptDose)<6){

      ##Less than 2 cohorts treated
      if(OptDose>1){

        if(OptDose<length(Dose)){

          if(sum(Doses==(OptDose+1))<6){
            OptDose=OptDose-1
          }else{
            if(VEC1[OptDose-1]>VEC1[OptDose+1]){
              OptDose=OptDose-1
            }else{
              OptDose=OptDose+1
            }


          }

        }else{

          OptDose=OptDose-1
        }

      }else{
        OptDose=OptDose+1
      }

    }
  }






  return(OptDose)


}


