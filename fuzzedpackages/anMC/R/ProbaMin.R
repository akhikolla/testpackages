##########
#' @title Probability of exceedance of minimum of Gaussian vector
#'
#'
#' @description Computes \eqn{P(min X \le threshold)}
#' with choice of algorithm between ANMC_Gauss and MC_Gauss.
#' By default, the computationally expensive sampling parts are computed with the Rcpp functions.
# INPUT
#' @param cBdg computational budget.
#' @param threshold threshold.
#' @param mu mean vector.
#' @param Sigma covariance matrix.
#' @param E discretization design for the field. If \code{NULL}, a simplex-lattice design {n,n} is used, with \code{n=length(mu)}. In this case the choice of method=4,5 are not advised.
#' @param q number of active dimensions, it can be either \itemize{
#'          \item an integer: in this case the optimal \code{q} active dimension are chosen;
#'          \item a numeric vector of length 2: this is the range where to search for the best number of active dimensions;
#'          \item \code{NULL}: q is selected as the best number of active dimensions in the feasible range.
#' }
# It can be passed either as an integer or as numeric vector of length 2. The vector is the range where to search for the best number of active dimensions. If \code{NULL} q is selected as the best number of active dimensions in the feasible range.
#' @param pn coverage probability function evaluated with \code{mu}, \code{Sigma}. If \code{NULL} it is computed automatically.
#' @param lightReturn boolean, if \code{TRUE} light return.
#' @param method method chosen to select the active dimensions. See \code{\link{selectActiveDims}} for details.
#' @param verb level of verbosity (0-5), selects verbosity also for \code{\link{ANMC_Gauss}} (verb-1) and \code{\link{MC_Gauss}} (verb-1).
#' @param Algo choice of algorithm to compute the remainder Rq ("ANMC" or "MC").
#' @param trmvrnorm function to generate truncated multivariate normal samples, it must have the following signature \code{trmvrnorm(n,mu,sigma,upper,lower,verb)}, where \itemize{
#'        \item \code{n}: number of simulations;
#'        \item \code{mu}: mean vector of the Normal variable of dimension \eqn{d};
#'        \item \code{sigma}: covariance matrix of dimension \eqn{d x d};
#'        \item \code{upper}: vector of upper limits of length \code{d};
#'        \item \code{lower}: vector of lower limits of length \code{d};
#'        \item \code{verb}: the level of verbosity 3 basic, 4 extended.
#' }
#' It must return a matrix \eqn{d x n} of realizations. If not specified, the rejection sampler \code{\link{trmvrnorm_rej_cpp}} is used.
#' @param pmvnorm_usr function to compute core probability on active dimensions. Inputs: \itemize{
#' \item \code{lower:} the vector of lower limits of length \code{d}.
#' \item \code{upper:} the vector of upper limits of length \code{d}.
#' \item \code{mean:} the mean vector of length \code{d}.
#' \item \code{sigma:} the covariance matrix of dimension \code{d}.
#' }
#' returns a the probability value with attribute "error", the absolute error. Default is the function \code{\link[mvtnorm]{pmvnorm}} from the package \code{mvtnorm}.
#' @return A list containing
#' \itemize{
#'    \item{\code{probability}: }{The probability estimate}
#'    \item{\code{variance}: }{the variance of the probability estimate}
#'    \item{\code{q}:}{the number of selected active dimensions}
#' }
#' If \code{lightReturn=F} then the list also contains:
#' \itemize{
#'    \item{\code{aux_probabilities}: }{ a list with the probability estimates: \code{probability} the actual probability, \code{pq} the biased estimator \eqn{p_q}, \code{Rq} the conditional probability \eqn{R_q}}
#'    \item{\code{Eq}: }{the points of the design \eqn{E} selected for \eqn{p_q}}
#'    \item{\code{indQ}: }{the indices of the active dimensions chosen for \eqn{p_q}}
#'    \item{\code{resRq}: }{The list returned by the MC method used for \eqn{R_q}}
#'  }
#' @examples
#' \dontrun{
#' # Compute probability P(X \in [0,\infty]) with X~N(0,Sigma)
#' d<-200     # example dimension
#' mu<-rep(0,d)    # mean of the normal vector
#' # correlation structure (Miwa et al. 2003, Craig 2008, Botev 2016)
#' Sigma<-0.5*diag(d)+ 0.5*rep(1,d)%*%t(rep(1,d))
#' pANMC<-ProbaMin(cBdg=20, q=min(50,d/2), E=seq(0,1,,d), threshold=0, mu=mu, Sigma=Sigma,
#'  pn = NULL, lightReturn = TRUE, method = 3, verb = 2, Algo = "ANMC")
#' proba<-1-pANMC$probability
#'
#' # Percentage error
#' abs(1-pANMC$probability-1/(d+1))/(1/(d+1))
#'
#'
#' # Implement ProbaMin with user defined function for active dimension probability estimate
#' if(!requireNamespace("TruncatedNormal", quietly = TRUE)) {
#' stop("TruncatedNormal needed for this example to work. Please install it.",
#'      call. = FALSE)
#' }
#' # define pmvnorm_usr with the function mvNcdf from the package TruncatedNormal
#' pmvnorm_usr<-function(lower,upper,mean,sigma){
#'     pMET<-TruncatedNormal::mvNcdf(l = lower-mean,u = upper-mean,Sig = sigma,n = 5e4)
#'     res<-pMET$prob
#'     attr(res,"error")<-pMET$relErr
#'     return(res)
#' }
#' pANMC<-ProbaMin(cBdg=20, q=min(50,d/2), E=seq(0,1,,d), threshold=0, mu=mu, Sigma=Sigma,
#'  pn = NULL, lightReturn = TRUE, method = 3, verb = 2, Algo = "ANMC",pmvnorm_usr=pmvnorm_usr)
#' proba<-1-pANMC$probability
#'
#' # Percentage error
#' abs(1-pANMC$probability-1/(d+1))/(1/(d+1))
#'
#' # Implement ProbaMax with user defined function for truncated normal sampling
#'
#' if(!requireNamespace("tmg", quietly = TRUE)) {
#' stop("Package tmg needed for this example to work. Please install it.",
#'      call. = FALSE)
#' }
#' trmvrnorm_usr<-function(n,mu,sigma,upper,lower,verb){
#'  M<-chol2inv(chol(sigma))
#'  r=as.vector(M%*%mu)
#'
#'  if(all(lower==-Inf) && all(upper==Inf)){
#'    f<- NULL
#'    g<- NULL
#'  }else{
#'    if(all(lower==-Inf)){
#'      f<--diag(length(mu))
#'      g<-upper
#'      initial<-(upper-1)/2
#'    }else if(all(upper==Inf)){
#'      f<-diag(length(mu))
#'      g<- -lower
#'      initial<-2*(lower+1)
#'    }else{
#'      f<-rbind(-diag(length(mu)),diag(length(mu)))
#'      g<-c(upper,-lower)
#'      initial<-(upper-lower)/2
#'    }
#'  }
#'  reals_tmg<-tmg::rtmg(n=n,M=M,r=r,initial = initial,f=f,g=g)
#'
#'  return(t(reals_tmg))
#' }
#'
#' pANMC<-ProbaMin(cBdg=20, q=min(50,d/2), E=seq(0,1,,d), threshold=0, mu=mu, Sigma=Sigma,
#'  pn = NULL, lightReturn = TRUE, method = 3, verb = 2, Algo = "ANMC",trmvrnorm=trmvrnorm_usr)
#' proba<-1-pANMC$probability
#'
#' # Percentage error
#' abs(1-pANMC$probability-1/(d+1))/(1/(d+1))
#' }
#' @references Azzimonti, D. and Ginsbourger, D. (2018). Estimating orthant probabilities of high dimensional Gaussian vectors with an application to set estimation. Journal of Computational and Graphical Statistics, 27(2), 255-267. Preprint at \href{https://hal.archives-ouvertes.fr/hal-01289126}{hal-01289126}
#'
#' Azzimonti, D. (2016). Contributions to Bayesian set estimation relying on random field priors. PhD thesis, University of Bern.
#'
#' Chevalier, C. (2013). Fast uncertainty reduction strategies relying on Gaussian process models. PhD thesis, University of Bern.
#'
#' Dickmann, F. and Schweizer, N. (2014). Faster comparison of stopping times by nested conditional Monte Carlo. arXiv preprint arXiv:1402.0243.
#'
#' Genz, A. (1992). Numerical computation of multivariate normal probabilities. Journal of Computational and Graphical Statistics, 1(2):141--149.
#' @export
ProbaMin = function(cBdg,threshold,mu,Sigma,E=NULL,q=NULL,pn=NULL,lightReturn=T,method=4,verb=0,Algo="ANMC",trmvrnorm = trmvrnorm_rej_cpp,pmvnorm_usr=pmvnorm){

  # initialize parameters
  n<-length(mu)

  if(is.null(pn))
    pn<-pnorm((mu-threshold)/sqrt(diag(Sigma)))

  if(is.null(E)){
    E<-seq(0,1,length.out = length(mu))
  }

  if(!is.matrix(E))
    E<-as.matrix(E)

  # Select q and the active dimensions (if q is a number we select q active dims)
#  q<-min(q,length(unique(pn))-1)
  indQ<-selectActiveDims(q=q,E=E,threshold=threshold,mu=mu,Sigma=Sigma,pn=(1-pn),method=method,pmvnorm_usr=pmvnorm_usr)

  # if q was given as a range here we reinitialize it as the number of active dims
  q<-length(indQ)

  if(verb>=1){
    cat("Computed indQ points, q = ",q,"\n")
  }

  Eq<-E[indQ,]

  # compute muEq
  muEq<-mu[indQ]

  # compute k(Eq,Eq)
  KEq<-Sigma[indQ,indQ]

  while(attr(chol(KEq,pivot=T),"rank")!=q){
    # ReSelect the appropriate q points
    indQ<-selectActiveDims(q=q,E=E,threshold=threshold,mu=mu,Sigma=Sigma,pn=pn,method=method,pmvnorm_usr=pmvnorm_usr)
    #  Eq<-E[indQ]
    Eq<-E[indQ,]
    # compute muEq
    muEq<-mu[indQ]
    # compute k(Eq,Eq)
    KEq<-Sigma[indQ,indQ]
    if(verb>=2){
      cat("Covariance matrix non p.d., re-Computed indQ points\n")
    }
  }

  cholKeq<-chol(KEq)

  # Compute p'
  pPrime<- 1 - pmvnorm_usr(lower= rep(threshold,length(indQ)),upper =rep(Inf,length(indQ)),mean = muEq,sigma = KEq)
  #  sSize<-N
  if(verb>=1){
    cat("Computed pPrime = ",pPrime,"\n")
  }
  if((1-pPrime)<attr(pPrime,"error") || q==nrow(Sigma)){
    if(verb>=2){
      cat("pPrime close to 1: pPrime=",pPrime,", error=",attr(pPrime,"error"),
          "or active dimensions equal to size of problem: q=",q,"length(mu)=",nrow(Sigma),
          "\n")
    }
    if(lightReturn){
      res<- list(probability=as.vector(pPrime), variance=(2/7*attr(pPrime,"error"))^2,q=q)
    }else{
      res<- list(probability=as.vector(pPrime), variance=(2/7*attr(pPrime,"error"))^2,q=q,
                 aux_probabilities=list(probability=as.vector(pPrime),pq=pPrime,Rq=0),
                 Eq=Eq,indQ=indQ,resRq=NULL)
    }
    if(verb>=2){
      cat("Early return. \n")
    }
    return(res)
  }

  # problem = list defining the problem: with mandatory fields
  #         - muEq = mean vector of X^{q}
  #         - sigmaEq = covariance matrix of X^q
  #         - threshold = threshold
  #         - muEmq = mean vector of X^{-q}
  #         - wwCondQ = ``weights'' for X^{-q} | X^q
  #         - sigmaCondQChol = Cholesky factorization of the conditional covariance matrix X^{-q} | X^q

  problem_Gauss<-list(muEq=muEq,sigmaEq=KEq,threshold=threshold,muEmq=mu[-indQ])
  invXX<-chol2inv(cholKeq)
  sigmaXY<-Sigma[indQ,-indQ]
  problem_Gauss$wwCondQ<-crossprod(sigmaXY,invXX)
  sigmaYcondX<-Sigma[-indQ,-indQ]-problem_Gauss$wwCondQ%*%sigmaXY
  problem_Gauss$sigmaCondQChol<-chol(sigmaYcondX)


  if(Algo=="ANMC"){
    if(verb>=1){
      cat("Starting ANMC. \n")
    }
    resMCQMC<-ANMC_Gauss(compBdg = cBdg,problem = problem_Gauss,delta=0.45,type="m",trmvrnorm = trmvrnorm,typeReturn=2,verb=max(0,verb-1))
  }else{
    if(verb>=1){
      cat("Starting MC. \n")
    }
    resMCQMC<-MC_Gauss(compBdg = cBdg,problem = problem_Gauss,delta=0.2,type="m",trmvrnorm = trmvrnorm,typeReturn=2,verb=max(0,verb-1))
  }


  proba<- pPrime + resMCQMC$estim*(1-pPrime)

  varpPrime<-(attr(pPrime,"error")/3)^2
  vars<- (1-resMCQMC$estim)^2*varpPrime+resMCQMC$varEst*(1-pPrime)^2 +resMCQMC$varEst*varpPrime

  if(lightReturn){
    res<-list(probability=as.vector(proba), variance=vars,q=q)
  }else{
    res<- list(probability=as.vector(proba), variance=vars,q=q,
               aux_probabilities=list(probability=as.vector(proba),pq=pPrime,Rq=resMCQMC$estim),
               Eq=Eq,indQ=indQ,resRq=resMCQMC)
  }

  return(res)
}

