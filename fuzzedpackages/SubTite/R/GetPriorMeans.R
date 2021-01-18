#' Calibrates prior means for Dose Finding Trial
#'
#'  Uses the clinician elicited prior reference probabilities for each subgroup and dose to obtain prior means for the Bayesian logistic regression model used in the SubTite trial design.
#' @param Clinician #Groups X #Doses matrix containing the elicited prior toxicity probabilities at the reference time for each dose and subgroup.
#' @param Dose Vector containing standardized doses.
#' @return Returns the nonlinear regression model whos parameter estimates will be used as prior means for the SubTITE Design.
#' @references
#' [1] Chapple and Thall (2017), Subgroup-specific dose finding in phase I clinical trials based on time to toxicity allowing adaptive subgroup combination
#' @examples
#' ##Specify elicited reference toxicity probabilities
#' Clinician = matrix(c(.2,.3,.4,.5,.6,.1,.2,.3,.4,.5,.05,.1,.15,.2,.3),byrow=TRUE,nrow=3)
#' Dose=sort(rnorm(5))
#' GetPriorMeans(Clinician,Dose)
#' @export
GetPriorMeans = function(Clinician,Dose){

  X=Clinician
  nGroups=nrow(X)

  Y = rep(NA,length(X))

  for(m in 1:nrow(X)){

    for(k in 1:ncol(X)){

      Y[(m-1)*ncol(X)+k] = log(X[m,k]/(1-X[m,k]))

    }


  }

  ##Now we have our Y lets make X
  ##It's going to be structured intercept, group ints, slope, group slopes
  ##Number of Groups
  G=nrow(X)
  G2=G-1
  D=length(Dose)

  COV = matrix(rep(0,length(Y)*2*G),nrow=length(Y))
  COV[,1]=1


  DOSEVEC=Y

  for(m in 1:G){
    DOSEVEC[((m-1)*length(Dose)+1):(m*length(Dose))]=Dose
  }








  if(nGroups==2){

    Group=c(rep(0,length(Dose)),rep(1,length(Dose)))
    Group1=Group==1







    m1 <- nls(Y ~ alpha + alpha1*Group1+exp(beta+beta1*Group1)*DOSEVEC,
              start = list(alpha=0,alpha1=0,beta=0,beta1=0))


    print(m1)



  }






  if(nGroups==3){

    Group=c(rep(0,length(Dose)),rep(1,length(Dose)),rep(2,length(Dose)))
    Group1=Group==1
    Group2=Group==2







    m1 <- nls(Y ~ alpha + alpha1*Group1+alpha2*Group2+exp(beta+beta1*Group1+beta2*Group2)*DOSEVEC,
              start = list(alpha=0,alpha1=0, alpha2=0,beta=0,beta1=0,beta2=0))


    print(m1)



  }





  if(nGroups==4){

    Group=c(rep(0,length(Dose)),rep(1,length(Dose)),rep(2,length(Dose)),rep(3,length(Dose)))
    Group1=Group==1
    Group2=Group==2
    Group3=Group==3







    m1 <- nls(Y ~ alpha + alpha1*Group1+alpha2*Group2+alpha3*Group3+exp(beta+beta1*Group1+beta2*Group2+beta3*Group3)*DOSEVEC,
              start = list(alpha=0,alpha1=0, alpha2=0,alpha3=0,beta=0,beta1=0,beta2=0,beta3=0))


    print(m1)



  }











  if(nGroups==5){

    Group=c(rep(0,length(Dose)),rep(1,length(Dose)),rep(2,length(Dose)),rep(3,length(Dose)),rep(4,length(Dose)))
    Group1=Group==1
    Group2=Group==2
    Group3=Group==3
    Group4=Group==4







    m1 <- nls(Y ~ alpha + alpha1*Group1+alpha2*Group2+alpha3*Group3+alpha4*Group4+exp(beta+beta1*Group1+beta2*Group2+beta3*Group3+beta4*Group4)*DOSEVEC,
              start = list(alpha=0,alpha1=0, alpha2=0,alpha3=0,alpha4=0,beta=0,beta1=0,beta2=0,beta3=0,beta4=0))


    print(m1)



  }








  if(nGroups==6){

    Group=c(rep(0,length(Dose)),rep(1,length(Dose)),rep(2,length(Dose)),rep(3,length(Dose)),rep(4,length(Dose)),rep(5,length(Dose)))
    Group1=Group==1
    Group2=Group==2
    Group3=Group==3
    Group4=Group==4
    Group5=Group==5







    m1 <- nls(Y ~ alpha + alpha1*Group1+alpha2*Group2+alpha3*Group3+alpha4*Group4+alpha5*Group5+exp(beta+beta1*Group1+beta2*Group2+beta3*Group3+beta4*Group4+beta5*Group5)*DOSEVEC,
              start = list(alpha=0,alpha1=0, alpha2=0,alpha3=0,alpha4=0,alpha5=0,beta=0,beta1=0,beta2=0,beta3=0,beta4=0,beta5=0))


    print(m1)



  }


  if(nGroups>6){
    cat("Code only supports up to 6 subgroups, contact maintainer if you desire more")

  }









}
