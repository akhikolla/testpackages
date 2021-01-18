#' @useDynLib BayesMFSurv
#' @importFrom stats dgamma runif as.formula model.frame model.matrix model.response na.omit na.pass
#' @import grDevices
#' @import graphics
#' @import RcppArmadillo
#' @importFrom Rcpp sourceCpp
#' @importFrom MCMCpack riwish
#' @importFrom coda mcmc
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @export


#' @title mfsurv
#' @description \code{mfsurv} fits a parametric Bayesian MF model via Markov Chain Monte Carlo (MCMC) to estimate the misclassification in the first stage
#'  and the hazard in the second stage.
#' @param formula a formula in the form Y ~ X1 + X2... | C ~ Z1 + Z2 ... where Y is the duration until failure or censoring, and C is a binary indicator of observed failure.
#' @param Y0 the elapsed time since inception until the beginning of time period (t-1).
#' @param data list object of data.
#' @param N number of MCMC iterations.
#' @param burn burn-ins to be discarded.
#' @param thin thinning to prevent autocorrelation of chain of samples by only taking the n-th values.
#' @param w size of the slice in the slice sampling for (betas, gammas, lambda). The default is c(1,1,1). This value may be changed by the user to meet one's needs.
#' @param m limit on steps in the slice sampling. The default is 10. This value may be changed by the user to meet one's needs.
#' @param form type of parametric model distribution to be used. Options are "Exponential" or "Weibull". The default is "Weibull".
#' @param na.action a function indicating what should happen when NAs are included in the data. Options are "na.omit" or "na.fail". The default is "na.omit".
#' @return mfsurv returns an object of class \code{"mfsurv"}.
#'
#' A \code{"mfsurv"} object has the following elements:
#' \item{Y}{the vector of `Y'.}
#' \item{Y0}{the vector of `Y0'.}
#' \item{C}{the vector of `C'.}
#' \item{X}{matrix X's variables.}
#' \item{Z}{the vector of `Z'.}
#' \item{betas}{data.frame,  X.intercept and X variables.}
#' \item{gammas}{data.frame, Z.intercept and Z variables.}
#' \item{lambda}{integer.}
#' \item{post}{integer.}
#' \item{iterations}{number of MCMC iterations.}
#' \item{burn_in}{burn-ins to be discarded.}
#' \item{thinning}{integer.}
#' \item{betan}{integer, length of posterior sample for betas.}
#' \item{gamman}{integer, length of posterior sample for gammas.}
#' \item{distribution}{character, type of distribution.}
#' \item{call}{the call.}
#' \item{formula}{description for the model to be estimated.}
#' @examples
#' set.seed(95)
#' bgl <- Buhaugetal_2009_JCR
#' bgl <- subset(bgl, coupx == 0)
#' bgl <- na.omit(bgl)
#' Y   <- bgl$Y
#' X   <- as.matrix(cbind(1, bgl[,1:7]))
#' C   <- bgl$C
#' Z1  <- matrix(1, nrow = nrow(bgl))
#' Y0  <- bgl$Y0
#' model1 <- mfsurv(Y ~ X | C ~ Z1, Y0 = Y0,
#'                 N = 50,
#'                 burn = 20,
#'                 thin = 15,
#'                 w = c(0.1, .1, .1),
#'                 m = 5,
#'                 form = "Weibull",
#'                 na.action = 'na.omit')
#' @export

mfsurv <-function(formula, Y0, data = list(), N, burn, thin, w = c(1,1,1), m = 10,
                  form = c("Weibull", "Exponential"), na.action=c("na.omit","na.fail")){

  if (missing(na.action)){na.action <- "na.omit"}
  if (missing(Y0))warning("Y0: elapsed time since inception missing")
  if (missing(N))warning("N: number of iterations missing")
  if (missing(burn))warning("burn: number of burn-ins missing")
  if (missing(thin))warning("thin: thinning interval missing")

  equations<-as.character(formula)
  formula1 <- paste(strsplit(equations[2], "|", fixed = TRUE)[[1]][1],sep="")
  formula2 <- paste(strsplit(equations[2], "|", fixed = TRUE)[[1]][2],equations[1],equations[3],sep="")
  mf1 <- model.frame(formula = as.formula(formula1), data = data, na.action = na.pass)
  mf2 <- model.frame(formula = as.formula(formula2), data = data, na.action = na.pass)
  X <- model.matrix(attr(mf1, "terms"), data = mf1)
  Z <- model.matrix(attr(mf2, "terms"), data = mf2)
  C <- model.response(mf2)
  Y <- model.response(mf1)
  dataset <- as.data.frame(cbind(Y,C,X,Z))

  if(na.action == "na.omit"){
    dataset <- data.frame(na.omit(dataset))
    Y <- as.matrix(dataset[,1], ncol = 1)
    C <- as.matrix(dataset[,2], ncol = 1)
    X <- data.frame(dataset[,3:(ncol(X)+2)])
    names(X) <- c("X.intercept", colnames(X[,2:ncol(X)]))
    Z <- data.frame(dataset[,(ncol(X)+3):(ncol(X)+ncol(Z)+2)])
    names(Z) <- c("Z.intercept", colnames(Z[,2:ncol(Z)]))
    Y0 <- as.numeric(Y0)
    N <- as.numeric(N)
    burn <- as.numeric(burn)
    thin <- as.numeric(thin)
    m <- as.numeric(m)
    w <- as.vector(w)
    form <- as.character(form)
    na.action <- na.action
    est <- bayes.mfsurv.default(Y, Y0, C, X, Z, N, burn, thin, w, m, form, na.action)
    est$call <- match.call()
    est$formula <- formula
    est
  }else{
    if(na.action == "na.fail"){
      if(all(is.numeric(dataset))){
        Y <- as.matrix(dataset[,1], ncol = 1)
        C <- as.matrix(dataset[,2], ncol = 1)
        X <- data.frame(dataset[,3:(ncol(X)+2)])
        names(X) <- c("X.intercept", colnames(X[,2:ncol(X)]))
        Z <- data.frame(dataset[,(ncol(X)+3):(ncol(X)+ncol(Z)+2)])
        names(Z) <- c("Z.intercept", colnames(Z[,2:ncol(Z)]))
        Y0 <- as.numeric(Y0)
        dataset$Y0 <- Y0
        N <- as.numeric(N)
        burn <- as.numeric(burn)
        thin <- as.numeric(thin)
        m <- as.numeric(m)
        w <- as.vector(w)
        form <- form
        na.action <- na.action
        est <- bayes.mfsurv.default(Y, Y0, C, X, Z, N, burn, thin, w, m, form, na.action)
        est$call <- match.call()
        est$formula <- formula
        est
      }else{
        Y <- as.numeric(Y)
        Y0 <- as.numeric(Y0)
        C <- as.numeric(C)
        X <- as.matrix(X)
        Z <- as.matrix(Z)
        if(any(is.na(Y))) warning("Time indicator contains missing values")
        if(any(is.na(C))) warning("Censoring indicator contains missing values")
        if(any(is.na(X))) warning("Explanatory variable(s) in misclassification stage contains missing values")
        if(any(is.na(Z))) warning("Explanatory variable(s) in survival stagef contains missing values")
      }
    }
  }
}


#' @title mfsurv.stats
#' @description A function to calculate the deviance information criterion (DIC) for fitted model objects of class \code{mfsurv}
#' for which a log-likelihood can be obtained, according to the formula \emph{DIC = -2 * (L - P)},
#' where \emph{L} is the log likelihood of the data given the posterior means of the parameter and
#' \emph{P} is the  estimate of the effective number of parameters in the model.
#' @param object an object of class \code{mfsurv}, the output of \code{mfsurv()}.
#' @return list.
#' @examples
#' set.seed(95)
#' bgl <- Buhaugetal_2009_JCR
#' bgl <- subset(bgl, coupx == 0)
#' bgl <- na.omit(bgl)
#' Y   <- bgl$Y
#' X   <- as.matrix(cbind(1, bgl[,1:7]))
#' C   <- bgl$C
#' Z1  <- matrix(1, nrow = nrow(bgl))
#' Y0  <- bgl$Y0
#' model1 <- mfsurv(Y ~ X | C ~ Z1, Y0 = Y0,
#'                 N = 50,
#'                 burn = 20,
#'                 thin = 15,
#'                 w = c(0.1, .1, .1),
#'                 m = 5,
#'                 form = "Weibull",
#'                 na.action = 'na.omit')
#'
#' mfsurv.stats(model1)
#' @export

mfsurv.stats = function(object){

  #Calculate L
  X <- as.matrix(object$X)
  Z <- as.matrix(object$Z)
  Y <- as.matrix(object$Y)
  Y0 <- as.matrix(object$Y0)
  C <- as.matrix(object$C)
  data <- as.data.frame(cbind(object$Y, object$Y0, object$C, object$X, object$Z))
  theta_post = cbind(object$gammas, object$betas, object$lambda)
  theta_hat = apply(theta_post, 2, mean)
  L = llFun(theta_hat,Y, Y0,C,X,Z,data)$llik
  #Calculate P
  S = nrow(theta_post) #S = number of iterations
  #Add up the log likelihoods of each iteration
  llSum = 0
  sum1 = 0
  sum2 = 0
  sum3 = 0
  #l <- as.matrix(NA, nrow=S, ncol=1)
  for (s in 1:S) {
    theta_s = as.matrix(theta_post[s,])
    ll <- llFun(theta_s,Y,Y0,C,X,Z,data)
    llSum <- llSum + ll$llik
    sum1 <- sum1 + ll$one
    sum2 <- sum2 + ll$two
    sum3 <- sum3 + ll$three
    #l[s,] <- llFun(theta_s,Y,Y0,C,X,Z,data)
  }
  P = 2 * (L - (1 / S * llSum))
  #Calculate DIC
  DIC <- -2 * (L - P)
  all <- sum1/S
  finite <- sum2/S
  small <- sum3/S
  #Return the results
  list(DIC = DIC, Loglik = L)
}



#' @title mfsurv.summary()
#' @description Returns a summary of a mfsurv object via \code{\link[coda]{summary.mcmc}}.
#' @param object an object of class \code{mfsurv}, the output of \code{\link{mfsurv}}.
#' @param parameter one of three parameters of the mfsurv output. Indicate either "betas", "gammas" or "lambda".
#' @return list. Empirical mean, standard deviation and quantiles for each variable.
#' @examples
#' set.seed(95)
#' bgl <- Buhaugetal_2009_JCR
#' bgl <- subset(bgl, coupx == 0)
#' bgl <- na.omit(bgl)
#' Y   <- bgl$Y
#' X   <- as.matrix(cbind(1, bgl[,1:7]))
#' C   <- bgl$C
#' Z1  <- matrix(1, nrow = nrow(bgl))
#' Y0  <- bgl$Y0
#' model1 <- mfsurv(Y ~ X | C ~ Z1, Y0 = Y0,
#'                 N = 50,
#'                 burn = 20,
#'                 thin = 15,
#'                 w = c(0.1, .1, .1),
#'                 m = 5,
#'                 form = "Weibull",
#'                 na.action = 'na.omit')
#'
#' mfsurv.summary(model1, "betas")
#' @export

mfsurv.summary <- function(object, parameter = c("betas", "gammas", "lambda")){

  if (parameter == "betas"){
    sum <- summary(mcmc(object$betas))
    return(sum)
  }
  if (parameter == "gammas"){
    sum <- summary(mcmc(object$gammas))
    return(sum)
  }
  if (parameter == "lambda"){
    sum <- summary(mcmc(object$lambda))
    return(sum)
  }
  class(res) <- "summary.bayesmf"
  res
}


