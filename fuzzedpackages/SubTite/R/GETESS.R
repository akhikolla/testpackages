#' Determines Prior ESS for fixed values of sigma_alpha^2 and sigmabeta^2
#'
#'Uses the prior means for the intercept and slope parameters and the number of doses to obtain an approximate prior ESS for the given prior variances. The user should calibrate varint and varbeta with varint>varbeta such that the ESS value is 1.
#' @param Dose Vector containing standardized doses.
#' @param meanmu Prior mean for baseline intercept.
#' @param meanslope Prior mean for baseline slope.
#' @param MeanInts Vector of prior means for the group specific intercept parameters.
#' @param MeanSlopes Vector of prior means for the group specific slope parameters.
#' @param varint Prior variance for the intercept parameters.
#' @param varbeta Prior variance for the slope parameters.
#' @return Returns the nonlinear regression model whos parameter estimates will be used as prior means for the SubTITE Design.
#' @references
#' [1] Chapple and Thall (2017), Subgroup-specific dose finding in phase I clinical trials based on time to toxicity allowing adaptive subgroup combination.
#' @examples
#' ###Specify the prior hypermeans
#' meanmu=-.5
#' meanslope=-.05
#' MeanInts = c(-.5,-.1)
#' MeanSlopes = c(.1,0)
#' Dose=sort(rnorm(5))
#' varint=5
#' varbeta=1
#' GetESS(Dose,meanmu,meanslope,MeanInts,MeanSlopes,varint,varbeta)
#' @export
GetESS=function(Dose,meanmu,meanslope,MeanInts,MeanSlopes,varint,varbeta){

  prob1=.9

  nDose=length(Dose)

  varint1=varint
  varbeta1=varbeta

  estBetaParams <- function(mu, var) {
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(mu+beta)
  }


  if(length(MeanInts)==1){
  GROUP=as.list(c(0,0))
  B=100000
  PROBS1 = matrix(rep(NA,B*nDose),nrow=B)
  PROBS2 = matrix(rep(NA,B*nDose),nrow=B)


  for(b in 1:B){

    slope=meanslope + rnorm(1,0,sqrt(varbeta))
    mu=meanmu + rnorm(1,0,sqrt(varint))


    Ints=MeanInts + rnorm(1,0,sqrt(varint1))

    Slopes=MeanSlopes + rnorm(1,0,sqrt(varbeta1))



    I1=rbinom(1,1,prob1)
    PROBS1[b,]=exp(exp(slope)*Dose+mu)

    PROBS2[b,]=exp(exp(slope+Slopes*I1)*Dose+mu+Ints*I1)


  }

  PROBS1=PROBS1/(1+PROBS1)

  PROBS2=PROBS2/(1+PROBS2)

  GROUP[[1]]=PROBS1

  GROUP[[2]]=PROBS2



  a=estBetaParams(colMeans(GROUP[[1]],na.rm=TRUE),apply(GROUP[[1]],2,var,na.rm=TRUE))

  b=estBetaParams(colMeans(GROUP[[2]],na.rm=TRUE),apply(GROUP[[2]],2,var,na.rm=TRUE))




  return(mean((a+b)/2));

  }



  if(length(MeanInts)==2){
    GROUP=as.list(c(0,0,0))
    B=100000
    PROBS1 = matrix(rep(NA,B*nDose),nrow=B)
    PROBS2 = matrix(rep(NA,B*nDose),nrow=B)
    PROBS3=PROBS1

    for(b in 1:B){

      slope=meanslope + rnorm(1,0,sqrt(varbeta))
      mu=meanmu + rnorm(1,0,sqrt(varint))


      Ints=MeanInts + rnorm(length(MeanInts),0,sqrt(varint1))

      Slopes=MeanSlopes + rnorm(length(MeanInts),0,sqrt(varbeta1))



      I1=rbinom(1,1,prob1)
      PROBS1[b,]=exp(exp(slope)*Dose+mu)

      PROBS2[b,]=exp(exp(slope+Slopes[1]*I1)*Dose+mu+Ints[1]*I1)

      I1=rbinom(1,1,prob1)

      PROBS3[b,]=exp(exp(slope+Slopes[2]*I1)*Dose+mu+Ints[2]*I1)



    }

    PROBS1=PROBS1/(1+PROBS1)

    PROBS2=PROBS2/(1+PROBS2)

    PROBS3=PROBS3/(1+PROBS3)

    GROUP[[1]]=PROBS1

    GROUP[[2]]=PROBS2

    GROUP[[3]]=PROBS3


    a=estBetaParams(colMeans(GROUP[[1]],na.rm=TRUE),apply(GROUP[[1]],2,var,na.rm=TRUE))

    b=estBetaParams(colMeans(GROUP[[2]],na.rm=TRUE),apply(GROUP[[2]],2,var,na.rm=TRUE))


    c=estBetaParams(colMeans(GROUP[[3]],na.rm=TRUE),apply(GROUP[[3]],2,var,na.rm=TRUE))





    return(mean((a+b+c)/3));

  }





  if(length(MeanInts)==3){
    GROUP=as.list(c(0,0,0,0))
    B=100000
    PROBS1 = matrix(rep(NA,B*nDose),nrow=B)
    PROBS2 = matrix(rep(NA,B*nDose),nrow=B)
    PROBS3=PROBS1
    PROBS4=PROBS3

    for(b in 1:B){

      slope=meanslope + rnorm(1,0,sqrt(varbeta))
      mu=meanmu + rnorm(1,0,sqrt(varint))


      Ints=MeanInts + rnorm(length(MeanInts),0,sqrt(varint1))

      Slopes=MeanSlopes + rnorm(length(MeanInts),0,sqrt(varbeta1))



      I1=rbinom(1,1,prob1)
      PROBS1[b,]=exp(exp(slope)*Dose+mu)

      PROBS2[b,]=exp(exp(slope+Slopes[1]*I1)*Dose+mu+Ints[1]*I1)

      I1=rbinom(1,1,prob1)

      PROBS3[b,]=exp(exp(slope+Slopes[2]*I1)*Dose+mu+Ints[2]*I1)


      I1=rbinom(1,1,prob1)

      PROBS4[b,]=exp(exp(slope+Slopes[3]*I1)*Dose+mu+Ints[3]*I1)


    }

    PROBS1=PROBS1/(1+PROBS1)

    PROBS2=PROBS2/(1+PROBS2)

    PROBS3=PROBS3/(1+PROBS3)


    PROBS4=PROBS4/(1+PROBS4)


    GROUP[[1]]=PROBS1

    GROUP[[2]]=PROBS2

    GROUP[[3]]=PROBS3

    GROUP[[4]]=PROBS4


    a=estBetaParams(colMeans(GROUP[[1]],na.rm=TRUE),apply(GROUP[[1]],2,var,na.rm=TRUE))

    b=estBetaParams(colMeans(GROUP[[2]],na.rm=TRUE),apply(GROUP[[2]],2,var,na.rm=TRUE))


    c=estBetaParams(colMeans(GROUP[[3]],na.rm=TRUE),apply(GROUP[[3]],2,var,na.rm=TRUE))


    d=estBetaParams(colMeans(GROUP[[4]],na.rm=TRUE),apply(GROUP[[4]],2,var,na.rm=TRUE))



    return(mean((a+b+c+d)/4));

  }





  if(length(MeanInts)==4){
    GROUP=as.list(c(0,0,0,0,0))
    B=100000
    PROBS1 = matrix(rep(NA,B*nDose),nrow=B)
    PROBS2 = matrix(rep(NA,B*nDose),nrow=B)
    PROBS3=PROBS1
    PROBS4=PROBS3
    PROBS5=PROBS3


    for(b in 1:B){

      slope=meanslope + rnorm(1,0,sqrt(varbeta))
      mu=meanmu + rnorm(1,0,sqrt(varint))


      Ints=MeanInts + rnorm(length(MeanInts),0,sqrt(varint1))

      Slopes=MeanSlopes + rnorm(length(MeanInts),0,sqrt(varbeta1))



      I1=rbinom(1,1,prob1)
      PROBS1[b,]=exp(exp(slope)*Dose+mu)

      PROBS2[b,]=exp(exp(slope+Slopes[1]*I1)*Dose+mu+Ints[1]*I1)

      I1=rbinom(1,1,prob1)

      PROBS3[b,]=exp(exp(slope+Slopes[2]*I1)*Dose+mu+Ints[2]*I1)


      I1=rbinom(1,1,prob1)

      PROBS4[b,]=exp(exp(slope+Slopes[3]*I1)*Dose+mu+Ints[3]*I1)


      I1=rbinom(1,1,prob1)

      PROBS5[b,]=exp(exp(slope+Slopes[4]*I1)*Dose+mu+Ints[4]*I1)


    }

    PROBS1=PROBS1/(1+PROBS1)

    PROBS2=PROBS2/(1+PROBS2)

    PROBS3=PROBS3/(1+PROBS3)


    PROBS4=PROBS4/(1+PROBS4)

    PROBS5=PROBS5/(1+PROBS5)



    GROUP[[1]]=PROBS1

    GROUP[[2]]=PROBS2

    GROUP[[3]]=PROBS3

    GROUP[[4]]=PROBS4

    GROUP[[5]]=PROBS5


    a=estBetaParams(colMeans(GROUP[[1]],na.rm=TRUE),apply(GROUP[[1]],2,var,na.rm=TRUE))

    b=estBetaParams(colMeans(GROUP[[2]],na.rm=TRUE),apply(GROUP[[2]],2,var,na.rm=TRUE))


    c=estBetaParams(colMeans(GROUP[[3]],na.rm=TRUE),apply(GROUP[[3]],2,var,na.rm=TRUE))


    d=estBetaParams(colMeans(GROUP[[4]],na.rm=TRUE),apply(GROUP[[4]],2,var,na.rm=TRUE))

    e=estBetaParams(colMeans(GROUP[[5]],na.rm=TRUE),apply(GROUP[[5]],2,var,na.rm=TRUE))


    return(mean((a+b+c+d+e)/5));

  }






  if(length(MeanInts)==4){
    GROUP=as.list(c(0,0,0,0,0))
    B=100000
    PROBS1 = matrix(rep(NA,B*nDose),nrow=B)
    PROBS2 = matrix(rep(NA,B*nDose),nrow=B)
    PROBS3=PROBS1
    PROBS4=PROBS3
    PROBS5=PROBS3
    PROBS6=PROBS3


    for(b in 1:B){

      slope=meanslope + rnorm(1,0,sqrt(varbeta))
      mu=meanmu + rnorm(1,0,sqrt(varint))


      Ints=MeanInts + rnorm(length(MeanInts),0,sqrt(varint1))

      Slopes=MeanSlopes + rnorm(length(MeanInts),0,sqrt(varbeta1))



      I1=rbinom(1,1,prob1)
      PROBS1[b,]=exp(exp(slope)*Dose+mu)

      PROBS2[b,]=exp(exp(slope+Slopes[1]*I1)*Dose+mu+Ints[1]*I1)

      I1=rbinom(1,1,prob1)

      PROBS3[b,]=exp(exp(slope+Slopes[2]*I1)*Dose+mu+Ints[2]*I1)


      I1=rbinom(1,1,prob1)

      PROBS4[b,]=exp(exp(slope+Slopes[3]*I1)*Dose+mu+Ints[3]*I1)


      I1=rbinom(1,1,prob1)

      PROBS5[b,]=exp(exp(slope+Slopes[4]*I1)*Dose+mu+Ints[4]*I1)


      I1=rbinom(1,1,prob1)

      PROBS6[b,]=exp(exp(slope+Slopes[5]*I1)*Dose+mu+Ints[5]*I1)


    }

    PROBS1=PROBS1/(1+PROBS1)

    PROBS2=PROBS2/(1+PROBS2)

    PROBS3=PROBS3/(1+PROBS3)


    PROBS4=PROBS4/(1+PROBS4)

    PROBS5=PROBS5/(1+PROBS5)



    PROBS6=PROBS6/(1+PROBS6)


    GROUP[[1]]=PROBS1

    GROUP[[2]]=PROBS2

    GROUP[[3]]=PROBS3

    GROUP[[4]]=PROBS4

    GROUP[[5]]=PROBS5

    GROUP[[6]]=PROBS6



    a=estBetaParams(colMeans(GROUP[[1]],na.rm=TRUE),apply(GROUP[[1]],2,var,na.rm=TRUE))

    b=estBetaParams(colMeans(GROUP[[2]],na.rm=TRUE),apply(GROUP[[2]],2,var,na.rm=TRUE))


    c=estBetaParams(colMeans(GROUP[[3]],na.rm=TRUE),apply(GROUP[[3]],2,var,na.rm=TRUE))


    d=estBetaParams(colMeans(GROUP[[4]],na.rm=TRUE),apply(GROUP[[4]],2,var,na.rm=TRUE))

    e=estBetaParams(colMeans(GROUP[[5]],na.rm=TRUE),apply(GROUP[[5]],2,var,na.rm=TRUE))

    f=estBetaParams(colMeans(GROUP[[6]],na.rm=TRUE),apply(GROUP[[6]],2,var,na.rm=TRUE))


    return(mean((a+b+c+d+e+f)/5));

  }


if(length(MeanInts)>4){
  cat("Code only supports up to 6 subgroups, contact maintainer if you desire more")

}



}

