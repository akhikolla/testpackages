#' @title MC estimate for the remainder
#'
#' @description Standard Monte Carlo estimate for \eqn{P(max X^{-q} >threshold | max X^{q}\le threshold)} or \eqn{P(min X^{-q} <threshold | min X^{q}\ge threshold)} where X is a normal vector. It is used for the bias correction in \code{\link{ProbaMax}} and \code{\link{ProbaMin}}.
# Input:
#' @param compBdg total computational budget in seconds.
#' @param problem list defining the problem with mandatory fields: \itemize{
#'         \item muEq = mean vector of \eqn{X^{q}};
#'         \item sigmaEq = covariance matrix of \eqn{X^q};
#'         \item threshold = threshold;
#'         \item muEmq = mean vector of \eqn{X^{-q}};
#'         \item wwCondQ = ``weights'' for \eqn{X^{-q} | X^q} [ the vector \eqn{\Sigma^{-q,q}(\Sigma^q)^{-1}}];
#'         \item sigmaCondQChol = Cholesky factorization of the conditional covariance matrix \eqn{\Sigma^{-q | q}}.
#'         }
#' @param delta total proportion of budget assigned to initial estimate (default 0.1), the actual proportion used might be smaller.
#' @param type type of excursion: "m", for minimum below threshold or "M", for maximum above threshold.
#' @param trmvrnorm function to generate truncated multivariate normal samples, it must have the following signature trmvrnorm(n,mu,sigma,upper,lower,verb), where \itemize{
#'        \item \code{n}: number of simulations;
#'        \item \code{mu}: mean vector of the Normal variable of dimension \eqn{d};
#'        \item \code{sigma}: covariance matrix of dimension \eqn{d x d};
#'        \item \code{upper}: vector of upper limits of length \code{d};
#'        \item \code{lower}: vector of lower limits of length \code{d};
#'        \item \code{verb}: the level of verbosity 3 basic, 4 extended.
#' }
#' It must return a matrix \eqn{d x n} of realizations. If not specified, the rejection sampler \code{\link{trmvrnorm_rej_cpp}} is used.
#' @param typeReturn integer: 0 (only the estimate) or 1 (heavy return with variance of the estimate, parameters of the estimator and computational times).
#' @param verb the level of verbosity, also sets the verbosity of trmvrnorm (to verb-1).
#' @param params system dependent parameters (if NULL they are estimated).
#' @return A list containing the estimated probability of excursion, see \code{typeReturn} for details.
#' @references Azzimonti, D. and Ginsbourger, D. (2018). Estimating orthant probabilities of high dimensional Gaussian vectors with an application to set estimation. Journal of Computational and Graphical Statistics, 27(2), 255-267. Preprint at \href{https://hal.archives-ouvertes.fr/hal-01289126}{hal-01289126}
#'
#' Azzimonti, D. (2016). Contributions to Bayesian set estimation relying on random field priors. PhD thesis, University of Bern.
#'
#' Genz, A. (1992). Numerical computation of multivariate normal probabilities. Journal of Computational and Graphical Statistics, 1(2):141--149.
#' @export
MC_Gauss<-function(compBdg,problem,delta=0.1,type="M",trmvrnorm=trmvrnorm_rej_cpp,typeReturn=0,verb=0,params=NULL){
  sizeX<-length(problem$muEq)
  sizeY<-length(problem$muEmq)

  if(type=="M"){
    upperTmvn = rep(problem$threshold,sizeX)
    lowerTmvn = rep(-Inf,sizeX)
    gg = function(x){return(max(x)>problem$threshold)}
  }else{
    upperTmvn = rep(Inf,sizeX)
    lowerTmvn = rep(problem$threshold,sizeX)
    gg = function(x){return(min(x)<problem$threshold)}
  }

  if(verb>0)
    cat("Starting MC... \n")

  if(is.null(params)){
    # Step 1: estimate Cx0, beta0 to get
    # (possibly) reasonable n0,m0
    if(verb>1)
      cat("Initialize parameters... ")
    # estimate Cx
    timeInPart1<-get_chronotime()
    simsX<-trmvrnorm(n = 1,mu = problem$muEq,sigma = problem$sigmaEq, upper = upperTmvn, lower = lowerTmvn, verb=(verb-1))
    time1SimX<-(get_chronotime()-timeInPart1)*1e-9

    ttX<-rep(0,20)
    ii<-seq(from=1,length.out=20,by=max(1,floor((compBdg*delta*0.4/time1SimX-20)/190)))
    for(i in seq(20)){
      timeIn<-get_chronotime()
      temp<-trmvrnorm(n = ii[i],mu = problem$muEq,sigma = problem$sigmaEq,upper = upperTmvn,lower = lowerTmvn,verb=(verb-1))
      ttX[i]<-(get_chronotime()-timeIn)*1e-9*1.03
      simsX<-cbind(simsX,temp)
    }
    Cx0<-unname(lm(ttX~ii+0)$coefficients[1])


    # estimate alpha and beta with lm
    timeIn<-get_chronotime()
    muYcondX<- problem$muEmq + problem$wwCondQ%*%(simsX[,1]-problem$muEq)
    simsYcX<- mvrnormArma(n=1,mu = muYcondX,sigma=problem$sigmaCondQChol,chol=1)
    time1SimYcX<-(get_chronotime()-timeIn)*1e-9

    tt<-rep(0,20)
    ii<-seq(from=1,length.out=20,by=max(1,floor((compBdg*delta*0.5/time1SimYcX-20)/190)))
    for(i in seq(20)){
      timeIn<-get_chronotime()
      muYcondX<- problem$muEmq + problem$wwCondQ%*%(simsX[,1]-problem$muEq)
      temp<- mvrnormArma(n=ii[i],mu = muYcondX,sigma=problem$sigmaCondQChol,chol=1)
      tt[i]<-(get_chronotime()-timeIn)*1e-9*1.03
      simsYcX<-cbind(simsYcX,temp)
    }
    lmmYcX<-lm(tt[2:20]~ii[2:20])
    alpha<-unname(lmmYcX$coefficients[1])
    beta0<-unname(lmmYcX$coefficients[2])

    #  timeBeta<-microbenchmark(problem$simulatorYcX(10,simsX[,1],problem$preProcessYcX))
    # beta0<-(quantile(timeBeta$time,probs = 0.99,names = F))*1e-9

    # estimate time to evaluate function
    timeG=rep(NA,154)
    for(i in seq(154)){
      iniT<-get_chronotime()
      gg(1)
      timeG[i]<-(get_chronotime()-iniT)
    }
    tEvalG<-quantile(sort(timeG)[-c(1,2,154,155)],probs = 0.99,names = F)*1e-9

    C_adj<-compBdg*delta -time1SimX - sum(ttX) - sum(tt)-time1SimYcX - sum(timeG)*1e-9

    if(verb>1){
      cat("Time passed:",compBdg*delta-C_adj,"\n",
          "tEvalG:",tEvalG,"\n",
          "beta0: ",beta0,"\n",
          "Cx0: ",Cx0,"\n",
          "alpha: ",alpha,"\n")
    }

    #  delta1<-1/(1-totCorXY(problem$SigmaXX,problem$SigmaYY,problem$SigmaXY))^2

    n0<-ncol(simsX)
  }else{
    if(verb>1)
      cat("Parameters already initialized. ")
    timeInPart1<-get_chronotime()
    Cx0<-params$Cx
    alpha<-params$alpha
    beta0<-params$beta
    tEvalG<-params$evalG
  }
  if(verb>1)
    cat("Done.\n")

  # derive nStar
  nStar<-round(compBdg/(Cx0+(beta0+tEvalG)))

  timePart1<-(get_chronotime()-timeInPart1)*1e-9
  if(verb>1){
    cat("Computational parameters: \n")
    cat("Cx0: ",Cx0,", beta0: ",beta0, ", nStar: ",nStar, "Time Part 1: ",timePart1,"(compBdg assigned: ",compBdg*delta,")\n")
  }
  # Step 3
  # we already have n0 sims for X and m0 sims of Y|X for each x sim
  # so we only need nStar-n0 simulations for X and mStar simulations
  # of Y conditional on the nStar-n0 xs, and mStar-m0 for the n0 Xs
  simsYcondXfull<-matrix(0,nrow = sizeY,ncol = 1)

  # generate the missing nStar -n0 simulations of X
  # re-estimate also Cx for debug reasons
  #  timeIn<-proc.time()
  if(is.null(params)){
    if(nStar>n0){
      simsXfull<-cbind(simsX,trmvrnorm(n = (nStar-n0),mu = problem$muEq,sigma = problem$sigmaEq,upper = upperTmvn,lower = lowerTmvn,verb=(verb-1)))
    }else{
      simsXfull<-simsX
      nStar<-n0
    }
  }else{
    simsXfull<-trmvrnorm(n = nStar,mu = problem$muEq,sigma = problem$sigmaEq,upper = upperTmvn,lower = lowerTmvn,verb=(verb-1))
  }
  # compute hatG
  estim<-0
  expYcondXfull<-rep(0,nStar)

  ########
  # you need to save gEval in a vector and fit the ones you have already done in it!
  ########
  gEval<-rep(0,nStar)
  for(j in seq(nStar)){
    #  muYcondX<-problem$muY+wwYcondX%*%(simsXfull[,j]-problem$muX)
    muYcondX<- problem$muEmq + problem$wwCondQ%*%(simsXfull[,j]-problem$muEq)
    simsYcondXfull[,1]<-mvrnormArma(n=1,mu = muYcondX,sigma=problem$sigmaCondQChol,chol=1)

    # now we have the simulations, we can compute all estimates
    gEval[j]<-gg(simsYcondXfull[,1])

    #    hatG<-hatG+expYcondXfull[j]
    if(j%%100==0){
      if((get_chronotime()-timeInPart1)*1e-9 >= compBdg*0.96)
        break
    }

  }
  nStar<-j
  #  betaFull<-max(mean(unname(betaFull)),0.001)
  gEval<-gEval[1:nStar]
  estim<-mean(gEval)
  timeTot<-(get_chronotime()-timeInPart1)*1e-9

  if(verb>=1){
    cat("MC computation finished.\n")
    cat("Total time: ",timeTot,"(compBdg: ",compBdg,")\n")
  }

  if(typeReturn==0){
    return(estim)
  }else{
    varHatG<-var(gEval)/nStar
    params<-list(n=nStar,Cx=Cx0,alpha=alpha,beta=beta0,evalG=tEvalG)
    results<-list(estim=estim,varEst=varHatG,params=params)
    results$times = list(part1=timePart1,total=timeTot)
    return(results)
  }
  return("fail")
}
