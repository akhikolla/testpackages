#' Gives true mean survival times for the control therapy.
#'
#' Returns the mean survival times for the control given efficacy and toxicity dose probability vector, distribution family and linear relationship, effiacy, toxicity and survival.
#' @param Family Time to event distribution. Options include: Exponential, Gamma, Weibull, Lognormal.
#' @param alpha Shape parameter or standard deviation of a lognormal distribution.
#' @param ProbC Probability of efficacy and toxicity for the control therapy.
#' @param betaC Linear term for efficacy, toxicity and beta_0 for the control groupar term for efficacy, toxicity and beta_0 for the control group.
#' @importFrom  stats sd
#' @references
#' [1] Chapple and Thall (2018).A Hybrid Phase I-II/III Clinical Trial Design Allowing Dose Re-Optimization in Phase III. Biometrics. In Press,
#' @examples
#'###Family of Distributions
#'Family="Gamma"
#'###Shape parameter
#'alpha=2
#'##True beta vector for efficacy, toxicity and intercept of the control treatment
#'betaC=c(.3,-.25,2.389)
#'##True efficacy and toxicity probability for control group
#'ProbC = c(.4,.15)
#'ReturnMeanControl(ProbC,betaC,Family,alpha)
#' @export
ReturnMeanControl = function(ProbC,betaC,Family,alpha){
PE=ProbC[1]
PT=ProbC[2]

beta=betaC

  if(Family=="Exponential"){
    ##Note: beta here is a 3-vector

    Means = rep(0,1)


    m=1

    Probs = c((1-PE[m])*(1-PT[m]),(1-PT[m])*PE[m], PT[m]*(1-PE[m]), PT[m]*PE[m])

    YE=0
    YT=0
    Means[m] = Means[m]+ Probs[1]





    Means[m] = Means[m]+ Probs[2]*exp(beta[1])



    Means[m] = Means[m]+ Probs[3]*exp(beta[2])




    Means[m] = Means[m]+ Probs[4]*exp(beta[1]+beta[2])




    z = as.list(c(0,0))
    Means=Means*exp(beta[3])

    return(Means)
  }


  if(Family=="Weibull"){




    ##Note: beta here is a 3-vector

    Means = rep(0,1)


    m=1

    Probs = c((1-PE[m])*(1-PT[m]),(1-PT[m])*PE[m], PT[m]*(1-PE[m]), PT[m]*PE[m])

    YE=0
    YT=0
    Means[m] = Means[m]+ Probs[1]





    Means[m] = Means[m]+ Probs[2]*exp(beta[1])



    Means[m] = Means[m]+ Probs[3]*exp(beta[2])




    Means[m] = Means[m]+ Probs[4]*exp(beta[1]+beta[2])




    z = as.list(c(0,0))
    Means=Means*exp(beta[3])*gamma(1+1/alpha)

    return(Means)





  }

  if(Family=="Lognormal"){
    sig=alpha




    ##Note: beta here is a 3-vector

    Means = rep(0,1)


    m=1

    Probs = c((1-PE[m])*(1-PT[m]),(1-PT[m])*PE[m], PT[m]*(1-PE[m]), PT[m]*PE[m])

    YE=0
    YT=0
    Means[m] = Means[m]+ Probs[1]





    Means[m] = Means[m]+ Probs[2]*exp(beta[1])



    Means[m] = Means[m]+ Probs[3]*exp(beta[2])




    Means[m] = Means[m]+ Probs[4]*exp(beta[1]+beta[2])




    z = as.list(c(0,0))
    Means=Means*exp(beta[3]+sig/2)

    return(Means)




  }

  if(Family=="Gamma"){




    ##Note: beta here is a 3-vector

    Means = rep(0,1)


    m=1

    Probs = c((1-PE[m])*(1-PT[m]),(1-PT[m])*PE[m], PT[m]*(1-PE[m]), PT[m]*PE[m])

    YE=0
    YT=0
    Means[m] = Means[m]+ Probs[1]





    Means[m] = Means[m]+ Probs[2]*exp(beta[1])



    Means[m] = Means[m]+ Probs[3]*exp(beta[2])




    Means[m] = Means[m]+ Probs[4]*exp(beta[1]+beta[2])




    z = as.list(c(0,0))
    Means=Means*alpha*exp(beta[3])

    return(Means)



  }


}
