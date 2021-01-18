#' Hyperparameter selection for iMOM prior density
#'
#' @description This function finds data specific hyperparameters for inverse
#' moment prior density so that the overlap between the iMOM prior and null
#' MLE density is \eqn{1/\sqrt p}. In this algorithm \code{r} is always chosen
#' to be equal to 1 and \code{tau} is found based on the mentioned overlap.
#' @param X The design matrix. \code{NA}'s should be removed and columns be
#' scaled. It is recommended that the \code{PreProcess} function is run first
#' and its output used for this argument. The columns are genes and rows
#' represent the observations. The column names are used as gene names.
#' @param resp For logistic regression models, it is the binary response
#' vector. For Cox proportional hazard model, this is a two column matrix
#' where the first column contains survival time vector and the second column
#' is the censoring status for each observation.
#' @param eff_size This is the expected effect size in the model for a
#' standardized design matrix, which is basically the coefficient value that is
#' expected to occur the most based on some prior knowledge.
#' @param nlptype Determines the type of nonlocal prior that is used in the
#' analyses. It can be "piMOM" for product inverse moment prior, or "pMOM" for
#' product moment prior. The default is set to piMOM prior.
#' @param iter The number of iteration needed to simulate from null model in
#' order to approximate the null MLE density.
#' @param mod_prior Type of prior used for model space. \code{uniform} is for
#' uniform binomial and \code{beta} is for beta binomial prior. In the former
#' case, both hyper parameters in the beta prior are equal to \code{1} but in
#' the latter case those two hyper parameters are chosen as explained in the
#' reference papers.
#' @param family Determines the type of data analysis. \code{logistic} is for
#' binary outcome data where logistic regression modeling is used whereas
#' \code{survival} is for survival outcome data using Cox proportional
#' hazard model.
#' @return It returns a list having following object:
#' \item{tau}{The hyperparameter for iMOM prior density function, calculated
#' using the proposed algorithm for the given dataset.}
#' @author Amir Nikooienejad
#' @references Nikooienejad, A., Wang, W., and Johnson, V. E. (2016). Bayesian
#' variable selection for binary outcomes in high dimensional genomic studies
#' using nonlocal priors. Bioinformatics, 32(9), 1338-1345.\cr\cr
#' Johnson, V. E., and Rossell, D. (2010). On the use of nonlocal prior
#' densities in Bayesian hypothesis tests. Journal of the Royal Statistical
#' Society: Series B (Statistical Methodology), 72(2), 143-170.
#' @examples
#' ### Simulating Logistic Regression Data
#' n <- 20
#' p <- 100
#' set.seed(321)
#' Sigma <- diag(p)
#' full <- matrix(c(rep(0.5, p * p)), ncol = p)
#' Sigma <- full + 0.5 * Sigma
#' cholS <- chol(Sigma)
#' Beta <- c(1, 1.8, 2.5)
#' X = matrix(rnorm(n*p, 1, 2), ncol = p)
#' X = X%*%cholS
#' beta <- numeric(p)
#' beta[c(1:length(Beta))] <- Beta
#' XB <- X%*%beta
#' probs <- as.vector(exp(XB)/(1+exp(XB)))
#' y <- rbinom(n,1,probs)
#' colnames(X) <- paste("gene_",c(1:p),sep="")
#' Xout <- PreProcess(X)
#' XX <- Xout$X
#' hparam <- HyperSelect(XX, y, iter = 1000, mod_prior = "beta",
#'                       family = "logistic")
#'
#' hparam$tau
HyperSelect <- function(X, resp, eff_size = 0.7, nlptype = "piMOM", 
                        iter = 10000, mod_prior=c("unif","beta"),
                        family = c("logistic", "survival")){
  
  # ================== Function Definitions =========================
  
  mydimom <- function(x, tau, r){
    out <- tau ^ (r / 2) / gamma(r / 2) * abs(x) ^ (-r - 1) * exp(-tau / x ^ 2)
    out[is.na(out)] <- 0
    return(out)
  }
  # ================================
  
  mydmom <- function(x, tau, r){
    out <- (2*pi)^(-0.5)*tau^(-r-0.5)*exp(-x^2/(2*tau))*x^(2*r)
    out[is.na(out)] <- 0
    return(out)
  }
  # ================================ 
  
  mypimom_tmp <- function(x, tau, r)
    if (x <= 0) 0.5 * pgamma(1 / x ^ 2, r / 2, tau) else 1-(0.5 * pgamma(1 / x ^ 2, r / 2, tau))
  mypimom <- Vectorize(mypimom_tmp)
  # ================================
  
  mypmom_tmp <- function(x, tau, r){
    a <- tau^(0.5-r)*(sqrt(tau)*pnorm(abs(x)/sqrt(tau))-abs(x)/sqrt(2*pi)*exp(-0.5*x^2/tau))
    if (x <= 0) 1-a else a
  }
  mypmom <- Vectorize(mypmom_tmp)
  # ================================
  
  Obj_Fun <- function(par, p, betalim, betasd, betam, mydimom, mypimom, mydmom, mypmom, nlptype){
    
    tau <- par
    r <- 1
    
    froot <- function(x, pars, nlptype){
      if(nlptype==0) mydimom(x, tau=pars, r=1) - dnorm(x, betam, betasd)
      if(nlptype==1) mydmom(x, tau=pars, r=1) - dnorm(x, betam, betasd)
    }
    
    pv1 <- tryCatch({pv <- uniroot(froot, pars = tau, nlptype = nlptype, interval = c(0.000001, 7))$root}
                   , error = function(err){
                     return(NA)
                   })
    if (is.na(pv1)) return (1000)
    
    pv <- c(-pv, pv)
    ov1 <- 1 - diff(pnorm(pv, betam, betasd))
    if (nlptype==0) ov2 <- diff(mypimom(pv, tau = tau, r = 1))
    if (nlptype==1) ov2 <- diff(mypmom(pv, tau = tau, r = 1))
    ov <- ov1 + ov2
    out <- abs(ov - 1 / sqrt(p))
    return(out)
  }
  # ========================== Main =================================
  if(nlptype=="piMOM") nlptype_int <- 0
  if(nlptype=="pMOM") nlptype_int <- 1
  XX <- X
  bincols <- which(apply(XX, 2, function(x) {all(na.omit(x) %in% 0:1)}))
  if (length(bincols)) XX <- XX[, -bincols]
  
  X <- XX
  dx <- dim(X)
  n <- dx[1]
  
  if(family=="logistic"){
    p <- dx[2] - 1
    
    cons <- 0;
    prp <- p / n
    ar <- 2 ^ n
    if (prp > 4 && ar < Inf){
      ac <- 0
      cons <- 0
      while (ar > ac) {
        cons <- cons + 1
        ac <- choose(p, cons)
      }
    }else{
      cons <- ceiling(log(p))
    }
    cons <- min(cons,ceiling(log(p)))
    sprob <- sum(resp) / n
    
  }
  
  if(family=="survival"){
    p <- dx[2]
    exmat <- cbind(resp,X)
    cons <- ceiling(log(p))
    
    csr <- 1 - mean(exmat[, 2])
  }
  
  if (mod_prior == "beta"){
    a <- cons; b <- p - a;
  }
  if (mod_prior == "unif"){
    a <- 1; b <- 1;
  }
  
  betalim <- 10
  set.seed(1234)
  
  if(family=="logistic"){ EsBeta <- null_mle_lreg(X, n, p, cons, a, b, sprob, iter)}
  if(family=="survival"){ EsBeta <- null_mle_cox(X, n, p, cons, a, b, csr, iter)}
  betasd <- sd(EsBeta)
  betam <- 0
  
  init <- 1
  res <- nlminb(init, objective = Obj_Fun, p = p, betalim = betalim,
                betasd = betasd, betam = betam, mydimom = mydimom,
                mypimom = mypimom, mydmom = mydmom, mypmom=mypmom,
                nlptype = nlptype_int, lower=0,upper=20)
  
  pars <- res$par
  pars <- round(pars, 2)
  eff <- round(eff_size^2, 2)
  pars <- min(pars, eff)
  return(list(tau = pars))
}
