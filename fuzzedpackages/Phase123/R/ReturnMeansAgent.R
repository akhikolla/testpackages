#' Gives true mean survival times for doses considered of the experimental agent.
#'
#' Returns the dose specific mean survival times for given efficacy and toxicity dose probability vector, distribution family and linear relationship between dose, effiacy, toxicity and survival.
#' @param Dose Vector of standardized doses considered in the trial.
#' @param PE True efficacy dose-toxicity vector.
#' @param PT True toxicity dose-toxicity vector.
#' @param beta True linear term for the rate or mean parameter
#' @param Family Time to event distribution. Options include: Exponential, Gamma, Weibull, Lognormal.
#' @param alpha Shape parameter or standard deviation of a lognormal distribution.
#' @useDynLib Phase123
#' @references
#' [1] Chapple and Thall (2018). A Hybrid Phase 12/3 Clinical Trial Design Allowing Dose Re-Optimization in Phase 3 Biometrics. Under Review.
#' @examples
#'##True Efficacy and Toxicity Probabilities
#'PT = c(.1,.15,.25,.35,.5)
#'PE=c(.2,.4,.6,.65,.7)
#'##Dose Levels considered
#'Dose = c(1,2,3,3.5,5)
#'Dose=(Dose-mean(Dose))/sd(Dose)
#'###Family of Distributions
#'Family="Gamma"
#'###Shape parameter ## Doesn't matter for exponential distribution
#'alpha=2
#'###True Beta vector
#'beta = c(.75,-.5, .3, -.25,2.143)
#'ReturnMeansAgent(PE,PT,beta,Dose,Family,alpha)
#'@export
ReturnMeansAgent = function(PE,PT,beta,Dose,Family,alpha){
Doses=Dose
if(Family=="Exponential"){


  Means = rep(0,length(Doses))
  Means1 = rep(0,4)



  for(m in 1:length(Means)){

    Probs = c((1-PE[m])*(1-PT[m]),(1-PT[m])*PE[m], PT[m]*(1-PE[m]), PT[m]*PE[m])

    YE=0
    YT=0
    Means[m] = Means[m]+ Probs[1]





    Means[m] = Means[m]+ Probs[2]*exp(beta[2])



    Means[m] = Means[m]+ Probs[3]*exp(beta[3])




    Means[m] = Means[m]+ Probs[4]*exp(beta[2]+beta[3])


  }

  z = as.list(c(0,0))

  Means=Means*exp(beta[1]*Doses+beta[4]*Doses^2+beta[5])


  return(Means)
}

if(Family=="Weibull"){




  Means = rep(0,length(Doses))
  Means1 = rep(0,4)



  for(m in 1:length(Means)){

    Probs = c((1-PE[m])*(1-PT[m]),(1-PT[m])*PE[m], PT[m]*(1-PE[m]), PT[m]*PE[m])

    YE=0
    YT=0
    Means[m] = Means[m]+ Probs[1]





    Means[m] = Means[m]+ Probs[2]*exp(beta[2])



    Means[m] = Means[m]+ Probs[3]*exp(beta[3])




    Means[m] = Means[m]+ Probs[4]*exp(beta[2]+beta[3])


  }

  z = as.list(c(0,0))
  Means=Means*exp(beta[1]*Doses+beta[4]*Doses^2+beta[5])*gamma(1+1/alpha)

  return(Means)





}

if(Family=="Lognormal"){



  Means = rep(0,length(Doses))
  Means1 = rep(0,4)

  sig=alpha


  for(m in 1:length(Means)){

    Probs = c((1-PE[m])*(1-PT[m]),(1-PT[m])*PE[m], PT[m]*(1-PE[m]), PT[m]*PE[m])

    YE=0
    YT=0
    Means[m] = Means[m]+ Probs[1]





    Means[m] = Means[m]+ Probs[2]*exp(beta[2])



    Means[m] = Means[m]+ Probs[3]*exp(beta[3])




    Means[m] = Means[m]+ Probs[4]*exp(beta[2]+beta[3])


  }

  z = as.list(c(0,0))
  Means=Means*exp(beta[1]*Doses+beta[4]*Doses^2+beta[5]+sig/2)

  return(Means)





}

if(Family=="Gamma"){




  Means = rep(0,length(Doses))
  Means1 = rep(0,4)



  for(m in 1:length(Means)){

    Probs = c((1-PE[m])*(1-PT[m]),(1-PT[m])*PE[m], PT[m]*(1-PE[m]), PT[m]*PE[m])

    YE=0
    YT=0
    Means[m] = Means[m]+ Probs[1]





    Means[m] = Means[m]+ Probs[2]*exp(beta[2])



    Means[m] = Means[m]+ Probs[3]*exp(beta[3])




    Means[m] = Means[m]+ Probs[4]*exp(beta[2]+beta[3])


  }

  z = as.list(c(0,0))
  Means=Means*alpha*exp(beta[1]*Doses+beta[4]*Doses^2+beta[5])

  return(Means)




}


}
