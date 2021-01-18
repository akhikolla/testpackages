#' @title Bivariate probit with partial observability
#' 
#' @description \code{BiProbitPartial} estimates a bivariate probit with partial 
#' observability model. 
#' 
#' The bivariate probit with partial observability model is defined as follows.
#' Let \eqn{i} denote the \eqn{i}th observation which takes values from \eqn{1}
#' to \eqn{N}, \eqn{X_1}{X[1]} be a covariate matrix of dimension 
#' \eqn{N \times k_1}{N x k[1]}, \eqn{X_2}{X[2]} be a covariate matrix of 
#' dimension \eqn{N \times k_2}{N x k[2]}, \eqn{X_{1i}}{X[1i]} be the \eqn{i}th
#' row of \eqn{X_1}{X[1]}, \eqn{X_{2i}}{X[2i]} be the \eqn{i}th row of 
#' \eqn{X_2}{X[2]}, \eqn{\beta_1}{\beta[1]} be a coefficient vector of length 
#' \eqn{k_1}{k[1]} and \eqn{\beta_2}{\beta[2]} be a coefficient vector of length 
#' \eqn{k_2}{k[2]}. Define the latent response for stage one to be 
#' \deqn{y_{1i}^\star = X_{1i} \beta_1 + \epsilon_{1i}}{y[1i]* = X[1i] \beta[1] + \epsilon[1i]} 
#' and stage two to be 
#' \deqn{y_{2i}^\star = X_{2i} \beta_2 + \epsilon_{2i}.}{y[2i]* = X[2i] \beta[2] + \epsilon[2i].} 
#' Note the stages do not need to occur sequentially. Define the outcome of 
#' the first stage to be
#' \eqn{y_{1i} = 1}{y[1i] = 1} if \eqn{y_{1i}^\star > 0}{y[1i]* > 0} and 
#' \eqn{y_{1i} = 0}{y[1i] = 0} if \eqn{y_{1i}^\star \leq 0}{y[1i]* <= 0}. 
#' Define the outcome of the second stage to be 
#' \eqn{y_{2i} = 1}{y[2i] = 1} if \eqn{y_{2i}^\star > 0}{y[2i]* > 0} and 
#' \eqn{y_{2i} = 0}{y[2i] = 0} if \eqn{y_{2i}^\star \leq 0}{y[2i]* <= 0}.
#'  The observed outcome is the product of the outcomes from the two stages 
#'  \deqn{z_{i} = y_{1i} y_{2i}.}{z[i] = y[1i] y[2i].} The pair 
#'  \eqn{(\epsilon_{1i},\epsilon_{2i})}{(\epsilon[1i],\epsilon[2i])} is distributed independently 
#'  and identically multivariate normal with means 
#'  \eqn{E[\epsilon_{1i}] = E[\epsilon_{2i}] = 0}{E(\epsilon[1i]) = E(\epsilon[2i]) = 0}, 
#'  variances 
#'  \eqn{Var[\epsilon_{1i}] = Var[\epsilon_{2i}] = 1}{Var(\epsilon[1i]) = Var(\epsilon[2i]) = 1},
#'  and correlation (or equivalently covariance) 
#'  \eqn{Cov(\epsilon_{1i},\epsilon_{2i}) = \rho}{Corr(\epsilon[1i],\epsilon[2i]) = \rho}.
#'  A more general structural representation is presented in Poirier (1980).
#'  
#'  The model can be estimated by Bayesian Markov Chain Monte Carlo (MCMC) or 
#'  frequentist maximum likelihood methods. The correlation parameter \eqn{\rho} can be
#'  estimated or fixed. The MCMC algorithm used is a
#'  block Gibbs sampler within Metropolis-Hastings scheme developed by 
#'  Rajbhandari (2014). The default maximum likelihood method is based off the 
#'  Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm. 
#'  A modification of the algorithm is 
#'  used to include box constraints for when \eqn{\rho} is estimated. See 
#'  \link[optimr]{optimr} for details.
#'
#' @param formula an object of class \link[Formula]{Formula}: a symbolic
#'  description of the model to be fitted. The details of model specification
#'  are given under 'Details'.
#' @param data an optional data frame, list or environment (or object coercible
#'  by \link{as.data.frame} to a data frame) containing the variables in the 
#'  model. If not found in \code{data}, the variables are taken from 
#'  \code{environment(formula)}, typically the environment from which 
#'  \code{BiProbitPartial} is called.
#' @param subset an optional vector specifying a subset of observations to be
#'  used in the fitting process.
#' @param na.action a function which indicates what should happen when the data
#'  contain \code{NA} observations. The default is set by the \code{na.action} setting of
#'  \link{options}, and is \link{na.fail} if that is unset. The 'factory-fresh'
#'  default is \link{na.omit}. Another possible value is \code{NULL}, no action.
#'  Value \link{na.exclude} can be useful.
#' @param philosophy a character string indicating the philosophy to be used
#'  for estimation. For Bayesian MCMC estimation \code{philosophy =
#'   "bayesian"} should be used. For frequentist maximum likelihood estimation 
#'   \code{philosophy = "frequentist"} should be used. The default is Bayesian
#'    MCMC estimation.
#' @param control a list of control parameters. See 'Details'.
#'
#' @details 
#'
#' Models for \code{BiProbitPartial} are specified symbolically. A typical
#'  model has the form \code{response ~ terms1 | terms2} where \code{response}
#'  is the name of the (numeric binary) response vector and \code{terms1} and \code{terms2}
#'  are each a series of terms which specifies a linear predictor for latent response 
#'  equations 1 and 2. A \code{terms1} specification of the form \code{first + second} 
#'  indicates all the terms in \code{first} together with all the terms
#'  in \code{second} with duplicates removed. A specification of the form
#'  \code{first:second} indicates the set of terms obtained by taking the
#'  interactions of all terms in \code{first} with all terms in \code{second}.
#'  The specification \code{first*second} indicates the cross of
#'  \code{first} and \code{second}. This is the same as 
#'  \code{first + second + first:second}. Likewise for \code{terms2}.
#'
#' A Formula has an implied intercept term for both equations. To remove the
#'  intercept from equation 1 use either \code{response ~ terms1 - 1 | terms2}
#'  or \code{response ~ 0 + terms1 | terms2}. It is analgous to remove the
#'  intercept from the equation 2.
#'  
#' If \code{philosophy = "bayesian"} is specified then the model is
#'  estimated by MCMC methods based on Rajbhandari (2014). The prior for the
#'  parameters in equations 1 and 2 is multivariate normal with mean \code{beta0}
#'  and covariance \code{B0}. The prior for \eqn{\rho} is truncated normal on the 
#'  interval \eqn{[-1,1]} with mean parameter \code{rho0} and variance parameter
#'  \code{v0} and is assumed to be apriori independent of the parameters in 
#'  equations 1 and 2. 
#'  
#'  If \code{philosophy = "frequentist"} then the model is
#'  estimated by frequentist maximum likelihood using \code{optimr} from the package
#'  \pkg{optimr}.
#' 
#' The \code{control} argument is a list that can supply the tuning parameters
#'  of the Bayesian MCMC estimation and frequentist maximum likelihood
#'  estimation algorithms. For frequentist maximum likelihood the \code{control} 
#'  argument is passed directly to \code{control} in the function \code{optimr}
#'   from the package \pkg{optimr}. If one wants to specify the \code{method} for
#'   the function \code{optimr} then \code{method} must be passed as an element 
#'   of \code{control}. See \link[optimr]{optimr} for further details.
#'   
#'  The 
#'  \code{control} argument can supply any of the following components for 
#'  Bayesian MCMC estimation.
#'  
#'  \describe{
#'  
#'  \item{beta}{Numeric vector or list of \code{nchains} 
#'   elements each a numeric vector supplying starting values for the coefficients
#'   in equations 1 and 2. For each vector, the first \eqn{k_1}{k[1]} values are for the 
#'   coefficients in the first equation. The second \eqn{k_2}{k[2]} values are for 
#'   the coefficients in the second equation. Default is \code{beta = numeric( k1 + k2 )}, 
#'   a vector of zeros.}
#'   
#'  \item{rho}{Numeric or list of \code{nchains} elements each a numeric starting 
#'  value for \eqn{\rho}. Default is \code{rho = 0}.}
#'  
#'  \item{fixrho}{Logical value to determine if \eqn{\rho} is estimated. If 
#'   \code{fixrho = TRUE} then \eqn{\rho} is fixed at value \code{rho}. Default 
#'   is \code{fixrho = FALSE}.}
#'  
#'  \item{S}{Number of MCMC iterations. Default is \code{S = 1000}. For
#'   \code{philosophy = "bayesian"} only.}
#'   
#'  \item{burn}{Number of initial pre-thinning MCMC iterations to remove after
#'   estimation. Default is \code{burn = floor(S/2)}, the floor of the number 
#'   of MCMC iterations divided by 2. For \code{philosophy = "bayesian"} only.}
#'   
#'  \item{thin}{Positive integer to keep every \code{thin} post-burn in MCMC draw
#'   and drop all others. Default is \code{thin = 1}, keep all post burn-in draws. 
#'   For \code{philosophy = "bayesian"} only.}
#'   
#'   \item{seed}{Positive integer for \code{nchains = 1} or list of \code{nchains} 
#'   elements each a positive integer fixing the seed of the random number generator. 
#'   Typically used for replication. Default is \code{seed = NULL}, no seed. For 
#'   \code{philosophy = "bayesian"} only.}
#'   
#'   \item{nchains}{Positive integer specifying the number of MCMC chains. Default is 
#'   \code{nchains = 1}. For \code{philosophy = "bayesian"} only.}
#'   
#'  \item{beta0}{Numeric vector supplying the prior mean for the coefficients of
#'   equations 1 and 2. The first \eqn{k_1}{k[1]} components are for the coefficients
#'   of equation 1. The second \eqn{k_2}{k[2]} components are for the coefficients of
#'   equation 2. Default is \code{beta0 = numeric( k1 + k2 )}, a vector of zeros. 
#'   For \code{philosophy = "bayesian"} only.}
#'  
#'  \item{B0}{Numeric matrix supplying the prior covariace of the parameters of
#'   equations 1 and 2. The first \eqn{k_1}{k[1]} rows are for the parameters of 
#'   equation 1. The second \eqn{k_2}{k[2]} rows are for the parameters of equation 
#'   2. Likewise for columns. If unspecified the default is set such that the
#'   inverse of \eqn{B0} is a zero matrix of dimension \eqn{(k_1+k_2) \times 
#'   (k_1+k_2)}{(k[1]+k[2]) x (k[1]+k[2])}, a 'flat' prior. For 
#'   \code{philosophy = "bayesian"} only.}
#'  
#'  \item{rho0}{Numeric value supplying a prior parameter for \eqn{\rho} which is
#'   the mean of a normal distribution that is truncated to the interval
#'   \eqn{[-1,1]}. Default is \code{rho0 = 0}. For
#'   \code{philosophy = "bayesian"} only.}
#'   
#'   \item{v0}{Numeric value supplying a prior parameter for \eqn{\rho} which is
#'   the variance of a normal distribution that is truncated to the interval
#'   \eqn{[-1,1]}. Default is \code{v0 = 1}. For
#'   \code{philosophy = "bayesian"} only.}
#'   
#'  \item{nu}{Numeric degrees of freedom parameter for setting the degrees of 
#'  freedom for \eqn{\rho}'s proposal t-distribution. Default is \code{nu = 10}.}
#'  
#'  \item{tauSq}{Numeric scaling parameter for scaling \eqn{\rho}'s proposal t-distribution.
#'  Default is \code{tauSq = 1}.}
#'  
#'  \item{P}{Determines how aggressive proposal draws for \eqn{\rho} are.
#'   Set to \code{P = 0} normal or \code{P = -1} for aggresive. See Rajbhandari
#'   (2014) and for details. Default is \code{P = 0}. For \code{philosophy =
#'   "bayesian"} only.}
#'   
#'   \item{trace}{Numeric value determining the value of intermediate reporting.
#'   A negative value is no reporting, larger positive values provide higher 
#'   degrees of reporting.}
#'  
#'  }
#'  
#'  Note: If the Bayesian MCMC chains appear to not be converging and/or frequentist
#'  maximum likelihood produces errors with \code{summary}, the model may be
#'  unidentified. One possible solution is to add regressors to the first equation 
#'  that are exluded from the second equation or visa-versa. See Poirier (1980) for
#'  more details.
#'  
#' 
#' @references
#' \cite{Poirier, Dale J. (1980). "Partial Observability in bivariate probit models" Journal of Econometrics 12, 209-217. (Identification)}
#' 
#' \cite{Rajbhandari, Ashish (2014). "Identification and MCMC estimation of bivariate probit model with partial observability." Bayesian Inference in Social Sciences (eds I. Jeliazkov and X. Yang). (MCMC algorithm)}
#'
#' @examples 
#' data('Mroz87',package = 'sampleSelection')
#' Mroz87$Z = Mroz87$lfp*(Mroz87$wage >= 5)
#'
#' f1 = BiProbitPartial(Z ~ educ + age + kids5 + kids618 + nwifeinc | educ + exper + city, 
#'      data = Mroz87, philosophy = "frequentist")
#' summary(f1)
#'
#' # Use the estimates from the frequenist philosophy as starting values
#' b1 = BiProbitPartial(Z ~ educ + age + kids5 + kids618 + nwifeinc | educ + exper + city, 
#'     data = Mroz87, philosophy = "bayesian", 
#'     control = list(beta = f1$par[1:(length(f1$par)-1)], rho = tail(f1$par,1)))
#' summary(b1)
#' 
#' \dontrun{#The example used in the package sampleSelection is likely unidentified for 
#' this model
#' f2 = BiProbitPartial(Z ~ educ + age + kids5 + kids618 + nwifeinc | educ, 
#'      data = Mroz87, philosophy = "frequentist") #crashes
#' summary(f2) #crashes (f2 non-existent)
#'
#' # Bayesian methods typically still work for unidentified models
#' b2 = BiProbitPartial(Z ~ educ + age + kids5 + kids618 + nwifeinc | educ, 
#'     data = Mroz87, philosophy = "bayesian", 
#'     control = list(beta = f1$par[1:(length(f1$par)-3)], rho = tail(f1$par,1)))
#' summary(b2)   
#' }
#'
#' @export
#' 
#'
#' @return \code{BiProbitPartial} returns an \eqn{S \times (k_1+k_2+1)
#'  \times nchains}{S x (k[1]+k[2]+1) x nchains} array of MCMC draws of primary 
#'  class \code{mcmc.list} and secondary class \code{BiProbitPartialb}, if 
#'  \code{philosophy = "bayesian"}. Each element in the 
#'  first dimension represents a MCMC draw. The first \eqn{k_1}{k[1]} elements 
#'  in the second dimension are draws for the coefficientss in the first 
#'  equation. The next \eqn{k_2}{k[2]} elements of the second dimension are draws 
#'  for the coefficients in the second equation. The last element of the second 
#'  dimension are draws for the correlation parameter. The elements of the 
#'  third dimension are the chains. If \eqn{\rho} was fixed 
#'  (\code{fixrho = TRUE}) then each draw for the last element in the second 
#'  dimension is returned as the value it was fixed at (the starting value, 
#'  \code{rho}).
#' 
#' If \code{philosophy = "frequentist"} a list equivalent to the
#'  output \link[optimr]{optimr} with primary class \code{optimrml} and secondary
#'  class \code{BiProbitPartialf}. 
#' 
#' @author \code{BiProbitPartial} was written by Michael Guggisberg. The majority of
#'  the MCMC estimation was written by Amrit Romana based on Rajbhandari
#'   (2014). The development of this package was partially funded by the 
#'   Institude for Defense Analyses (IDA).
#'   
#' 
BiProbitPartial <- function(formula, data, subset, na.action, philosophy = "bayesian", control = list()){
  
  mf <- match.call(expand.dots = FALSE)
  print(mf)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  
  f <- Formula::Formula(formula)
  mf[[1]] <- as.name("model.frame")
  mf$formula <- f
  mf <- eval(mf, parent.frame())
   
  Z <- stats::model.response(mf)
  X1 <- stats::model.matrix(f, data = mf, rhs = 1)
  X2 <- stats::model.matrix(f, data = mf, rhs = 2)
  par1names = paste0(colnames(X1),"_1")
  par2names = paste0(colnames(X2),"_2")
  par12names = c(par1names,par2names)
  
  k1 = ncol(X1)
  k2 = ncol(X2)
  k = k1+k2
  n = length(Z)
  ncontrol = names(control)
  
  # Check data values
  if(any(!is.element(Z,c(0,1))))stop("The response can only take on values of 0 and 1")
  
  # Need to also remove non alphabetic characters
  philosophy = tolower(philosophy)
  
  if(is.null(control$fixrho)){control$fixrho = F}
  else if(!is.logical(control$fixrho))stop("fixrho needs to be a logical")
  
  if(philosophy == "frequentist"){
      if(is.null(control$beta)){
        control$beta1 = numeric(k1)
        control$beta2 = numeric(k2)
      }else if(any(!is.numeric(control$beta)) | length(control$beta) != k)
        stop("Starting values for beta needs be a numeric vector of length k1 + k2")
      else if(length(control$beta) == k){
        control$beta1 = control$beta[1:k1]
        control$beta2 = control$beta[(k1+1):k]
        }else stop("Something is wrong with the starting values of beta")
      if(is.null(control$rho)){control$rho = 0}
      else if(abs(control$rho)>1)stop("Starting value for rho needs to be between -1 and 1")
      parnames = if(control$fixrho)c(par1names,par2names)else c(par1names,par2names,"rho")
  }else{
    if(is.null(control$nchains)){control$nchains = 1}
    else{
      if(abs(control$nchains - round(control$nchains))>=.Machine$double.eps^0.5 | control$nchains<1)
        stop("nchains should be a whole number greater than 0")
    }
    if(is.null(control$seed))control$seed = rep(0,control$nchains)
    else if(!is.element(length(control$seed),c(1,control$nchains)))
      stop("seed must be a vector of length 1 or nchains. Or seed must be unspecified.")
    if(control$seed[[1]] != 0 & length(control$seed) == 1)
      warning("If nchains>1 seed should be unspecified or a vector of length nchains")
    if(is.null(control$seed)){control$seed = 0}
    else if(any(sapply(1:length(control$seed),function(i){abs(control$seed[[i]] - round(control$seed[[i]]))>=.Machine$double.eps^0.5 | control$seed[[i]]<0})))
      stop("each element of seed should be a whole number greater than or equal to 0")
    if(length(control$seed)==1)
      control$seed = rep(control$seed,control$nchains)
    

    if(is.null(control$beta)){
      control$beta1 = rep(list(numeric(k1)),control$nchains)
      control$beta2 = rep(list(numeric(k2)),control$nchains)
    }else{
      if(!is.list(control$beta)){
        if(any(!is.numeric(control$beta))|length(control$beta) != k1 + k2)stop("Starting values for beta needs be numeric vector of lenght k1 + k2")  
        else{
          control$beta1 = rep(list(control$beta[1:k1]),control$nchains)
          control$beta2 = rep(list(control$beta[(k1+1):k]),control$nchains)
          control$beta = rep(list(control$beta),control$nchains)
        }
      }else if(!is.element(length(control$beta),c(1,control$nchains)))
        stop("Starting values for beta need to be a list of length 1 or nchains")
      else if(length(control$beta)==1){
        control$beta = rep(control$beta,control$nchains)
        control$beta1 = rep(list(control$beta[[1]][1:k1]),control$nchains)
        control$beta2 = rep(list(control$beta[[1]][(k1+1):k]),control$nchains)
      }else{
        control$beta1 = lapply(1:control$nchains,function(i)control$beta[[i]][1:k1])
        control$beta2 = lapply(1:control$nchains,function(i)control$beta[[i]][(k1+1):k])
      }
    }
    if(!is.list(control$rho)){
      if(is.null(control$rho)){control$rho = rep(0,control$nchains)}
      else if(abs(control$rho)>1)stop("Starting value for rho needs to be between -1 and 1")
    }else if(!is.element(length(control$rho),c(1,control$nchains)))stop("Starting value for rho needs to be a list of length 1 or nchains")
    else if(length(control$rho)==1)
      control$rho = rep(control$rho,control$nchains)
  
    parnames = if(control$fixrho)c(par1names,par2names)else c(par1names,par2names,"rho")
    
  }
  
  if(is.null(control$trace)){control$trace=0}
  

  if(philosophy == "bayesian"){
    # Set MCMC parameters
    if(is.null(control$nu)){control$nu = 10}
    else if(control$nu<=0)stop("nu needs to be greater than 0")
    if(is.null(control$P)){control$P = 0}
    else if(!is.element(control$P,c(0,-1)))warning("P should be 0 or -1")
    if(is.null(control$tauSq)){control$tauSq = 1}
    else if(control$tauSq)stop("tauSq should be greater than 0")
    if(is.null(control$S)){control$S = 1000}
    else if(abs(control$S - round(control$S))>=.Machine$double.eps^0.5 | control$S<1)
      stop("S should be a whole number greater than 0")
    if(is.null(control$burn)){control$burn = floor(control$S/2)}
    else{ 
      if(abs(control$burn - round(control$burn))>=.Machine$double.eps^0.5 | control$burn<1)
        stop("burn should be a whole number greater than 0")
      if(control$burn>= control$S)
        stop("burn should be less than S")
    }
    if(is.null(control$thin)){control$thin = 1}
    else if(abs(control$thin - round(control$thin))>=.Machine$double.eps^0.5 | control$thin<1)
      stop("thin should be a whole number greater than 0")
    # Set prior hyperparameters
    if(is.null(control$beta0))control$beta0 = numeric(k)
    else if(any(!is.numeric(control$beta0)) | length(control$beta0) != k)
      stop("beta0 needs be a numeric vector of length k1 + k2")
    if(is.null(control$B0)){control$B0inv = diag(0,k)}
    else{
      if(any(eigen(control$B0)$values <= 0 ))stop("B0 needs to be positive definite")
      if(!isSymmetric(control$B0))stop("B0 needs to be symmetric")
      if(any(!is.element(dim(control$B0),k)))stop("B0 needs to be of dimension k x k")
      control$B0inv = chol2inv(control$B0)}
    if(is.null(control$rho0)){control$rho0 = 0}
    else if(abs(control$rho0)>1)stop("rho0 needs to be between -1 and 1")
    if(is.null(control$v0)){control$v0 = 0.25^2}
    else if(control$v0<=0)stop("v0 needs to be greater than 0")
    } 
  else if(philosophy == "frequentist"){
    bmatch = match(ncontrol,c("nu","P","tauSq","beta0","B0","rho0","v0"))
      if(any(!is.na(bmatch)))
        warning(paste("The following objects were provided and will not be used
                      if philosophy is set to frequentist: ",paste(nargs[bmatch], collapse=", ")))
    if(is.null(control$hessian)){control$hessian=TRUE}
    else if(!is.logical(control$hessian)){stop("hessian needs to be logical")}
    if(is.null(control$method)){control$method = "Rvmmin"}
    control$maximize = T
    initialTheta = c(control$beta1,control$beta2)
    if(!control$fixrho)
      initialTheta = c(initialTheta,control$rho)
    names(initialTheta) = parnames
    }
  else stop("Philosophy needs to be bayesian or frequentist")
  
  # Estimate Model
  if(philosophy == "bayesian"){
    if(control$trace>=0)  
      cat("Estimating with Bayesian MCMC","\n")
    if(control$trace>2){
      cat("Control parameters are: ","\n")
      print(control)
    }
    
    out = sapply(1:control$nchains,function(i){
      if(control$trace>1)
        cat(paste0("Beginning chain ",i," of ", control$nchains),"\n")
      MCMC1(X1=as.matrix(X1),X2=as.matrix(X2),Z = as.matrix(Z),
                beta1 = as.matrix(control$beta1[[i]]), beta2 = as.matrix(control$beta2[[i]]), 
                rho = control$rho[[i]],fixrho = control$fixrho, S = control$S, 
                beta0 = as.matrix(control$beta0), B0inv = as.matrix(control$B0inv), 
                rho0 = control$rho0, v0=control$v0, 
                nu = control$nu, P = control$P, tauSq = control$tauSq, seed = control$seed[[i]])},simplify = "array")
    dimnames(out) = list(NULL,parnames,NULL)
    out = out[(control$burn+1):nrow(out),,]
    if(control$nchains==1){
      out = coda::mcmc.list(coda::mcmc(out[(1:floor(nrow(out)/control$thin))*control$thin,],start=control$burn+1,thin=control$thin))
    }else
      out = coda::mcmc.list(lapply(1:control$nchains,function(i)
        coda::mcmc(out[(1:floor(nrow(out)/control$thin))*control$thin,,i],start=control$burn+1,thin=control$thin)))
    
    class(out) = c(class(out),"BiProbitPartialb")
    
  }
  else{
    if(control$trace>=0)
      cat("Estimating with frequentist maximum likelihood","\n")
    if(control$trace>2){
      cat("Control parameters are: ","\n")
      print(control)
    }
  
    out = optimr::optimr(par=initialTheta,fn=llhood1,gr=grad1,
                 lower = if(control$fixrho)-Inf
                         else c(rep(-Inf,k),-1),
                 upper = if(control$fixrho) Inf
                         else c(rep(Inf,k),1),
                 method = control$method, hessian  = control$hessian,
                 control = control,
                 X1=as.matrix(X1),X2=as.matrix(X2),Z = as.matrix(Z), fixrho = control$fixrho, 
                 rho = control$rho)

    out$k1 = k1
    out$k2 = k2
    out$n = n
    if(control$fixrho)
      out$fixrhoVal = control$rho
    class(out) = c("optimrml","BiProbitPartialf")
  }
  
  return(out)
  
}


#' @title Summary method for class 'optimrml'
#'
#' @param object Object of class \code{optimrml}
#' @param ... unused
#'
#' @return matrix summary of estimates. The columns are \describe{
#' \item{Estimate}{Maximum likelihood point estimate}
#' 
#' \item{Std. Error}{Asymptotic standard error estimate of maximum likelihood point estimators using numerical hessian}
#' 
#' \item{z value}{z value for zero value null hypothesis using asymptotic standard error estimate}
#' 
#' \item{Pr(>|z|)}{P value for a two sided null hyptothesis test using the z value}
#' }
#' 
#' @examples 
#' data('Mroz87',package = 'sampleSelection')
#' Mroz87$Z = Mroz87$lfp*(Mroz87$wage >= 5)
#' 
#' f1 = BiProbitPartial(Z ~ educ + age + kids5 + kids618 + nwifeinc | educ + exper + city, 
#'      data = Mroz87, philosophy = "frequentist")
#' summary(f1)
#' 
#' b1 = BiProbitPartial(Z ~ educ + age + kids5 + kids618 + nwifeinc | educ + exper + city, 
#'     data = Mroz87, philosophy = "bayesian", 
#'     control = list(beta = f1$par[1:(length(f1$par)-1)], rho = tail(f1$par,1)))
#' summary(b1)
#' 
#' @export
#'
summary.optimrml = function(object,...){
  
  dots = match.call(expand.dots = TRUE)
  if(is.null(dots$digits))
    dots$digits = max(3, getOption("digits") - 3)
  
  coeff = object$par
  if(!is.null(object$hessian))
    se = sqrt(-diag(solve(object$hessian)))
  else
    se = rep(NA,length(coeff))
  zstat = abs(coeff/se)
  pval = stats::pnorm(zstat, lower.tail = FALSE)*2
  out = cbind(coeff,se,zstat,pval)
  colnames(out) = c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(out) = names(object$par)
  out
  
}

#' @title predict method for class 'BiProbitPartialb' 
#' 
#' @description Note this produces a 
#' Bayesian posterior predictive distribution. This accounts for estimation
#' uncertainty. If you desire a simple prediction that does not account for estimation
#' uncertainty then the frequentist philosophy should be used. If \code{nchains} is 
#' greater than 1 then the chains are combined. 
#'
#' @param object a object of class \code{BiProbitPartialb}
#' @param newdata a matrix of column dimension k1 + k2 where the first k1 columns
#' correspond to the predictors of the first equations and the second k2 
#' columns correspond to predictors of the second equation. If intercepts were
#' used they need to be explicitly input.
#' @param k1 a numeric declaring the number of covariates (including intercept) in the first equation
#' @param k2 a numeric declaring the number of covariates (including intercept) in the second equation
#' @param mRule a vector of length 1 or 2. This is the marginal decision rule
#' for classifying the outcomes for stages 1 and 2. Stage 1 is 
#' classified as 1 if the probability of stage 1 being 1 is greater
#' than or equal to \code{mRule[1]}. Likewise for stage 2.  If 
#' length of \code{mRule} is 1 then that value is recycled. The values of 
#' \code{mRule} must be between 0 and 1. The default value is \code{mRule = c(0.5,0.5)}.
#' @param jRule an optional numerical value between 0 and 1. If specified 
#' then the observable outcome (both stages being 1) is 1 if the joint
#' probability of both stages being 1 is greater than jRule. If jRule
#' is unspecified or set to \code{NULL} then the observable outcome is the product
#' of the marginal outcomes. The default value is \code{jRule = NULL}. Note, if
#' \code{jRule} is specified then the observable outcome might not equal the product of stages 1 and 2.
#' @param ... unused
#'
#' @return method \code{predict.bBiProbitPArtial} returns a data.frame with columns
#' \describe{
#' 
#' \item{linPredict1}{Predicted mean of the first stage latent outcome. This is tyically not 
#' interesting for a Bayesian analysis.}
#' 
#' \item{linPredict2}{Predicted mean of the second stage latent outcome. This is tyically not 
#' interesting for a Bayesian analysis.}
#' 
#' \item{p1.}{Probability the outcome of the first stage is 1}
#' 
#' \item{p.1}{Probability the outcome of the second stage is 1}
#'
#' \item{p00}{Probability the outcome of both stages is 0}
#' 
#' \item{p01}{Probability the outcome of the first stage is 0 and the second stage is 1}
#' 
#' \item{p10}{Probability the outcome of stage 1 is 1 and stage 2 is 0}
#' 
#' \item{p11}{Probability the outcome of both stages are 1 }
#' 
#' \item{yHat1}{Classification of the outcome for stage 1. This value
#' is 1 if \code{p1 >= mRule[1]} and 0 else}
#' 
#' \item{yHat2}{Classification of the outcome for stage 2. This value
#' is 1 if \code{p2 >= mRule[2]} and 0 else}
#' 
#' \item{ZHat}{Classification of the observable outcome.
#' If \code{jRule} is specified then this value is 1 if \code{p12 >= jRule}
#'  and 0 else. If \code{jRule} is unspecified then this value is the element-wise product of
#'  yHat1 and yHat2.}
#' 
#' }
#' 
#' @examples 
#' ##
#' # Perform a prediction with the same covariates the model is estimated with
#' ##
#' 
#' data('Mroz87',package = 'sampleSelection')
#' Mroz87$Z = Mroz87$lfp*(Mroz87$wage >= 5)
#' 
#' # Run the frequentist version first to get starting values
#' f1 = BiProbitPartial(Z ~ educ + age + kids5 + kids618 + nwifeinc | educ + exper + city, 
#'     data = Mroz87, philosophy = "frequentist")
#' 
#' b1 = BiProbitPartial(Z ~ educ + age + kids5 + kids618 + nwifeinc | educ + exper + city, 
#'     data = Mroz87, philosophy = "bayesian", 
#'     control = list(beta = f1$par[1:(length(f1$par)-1)], rho = tail(f1$par,1)))
#' 
#' library(Formula)
#' eqn = Formula::Formula( ~ educ + age + kids5 + kids618 + nwifeinc | educ + exper + city)
#' matrix1 = model.matrix(eqn, lhs = 0, rhs=1, data= Mroz87)
#' matrix2 = model.matrix(eqn, lhs = 0, rhs=2, data= Mroz87)
#' newdat = cbind(matrix1,matrix2) 
#' preds1 = predict(b1,newdat,k1 = dim(matrix1)[2],k2 = dim(matrix2)[2])
#' head(preds1)
#' preds2 = predict(b1,newdat,k1 = dim(matrix1)[2],k2 = dim(matrix2)[2], jRule = .25)
#' 
#' # Compare predicted outcome with realized outcome
#' head(cbind(Mroz87$Z,preds1$ZHat,preds2$ZHat),20)
#' 
#' @export
#'
predict.BiProbitPartialb = function(object, newdata, k1, k2, mRule = c(.5,.5), jRule = NULL,...){
  
  newdata = as.matrix(newdata)
  if(dim(newdata)[2]==1)
    newdata = t(newdata)
  
  if(dim(newdata)[2] != k1 + k2)
    stop("newdata must be a matrix with column length k1 + k2")
  
  if(length(mRule)>2 | any(mRule>1) | any(mRule<0))
    stop("mRule must be a value between 0 and 1 of length 1 or 2 or unspecified")
  else if(length(mRule)==1)
    mRule = rep(mRule,2)
  
  pardraws = as.matrix(Reduce(rbind,object))
  meandraws = colMeans(pardraws)
  
  linPredict1 = newdata[,1:k1] %*% meandraws[1:k1]
  linPredict2 = newdata[,1:k2 + k1] %*% meandraws[1:k2 + k1]
  
  mu1draws = newdata[,1:k1] %*% t(pardraws[,1:k1])
  mu2draws = newdata[,1:k2 + k1] %*% t(pardraws[,1:k2 + k1])
  
  p1. = rowMeans(stats::pnorm(0,mu1draws,lower.tail=F))
  p.1 = rowMeans(stats::pnorm(0,mu2draws,lower.tail=F))
  
  p11star = sapply(1:dim(pardraws)[1],
    function(i)pbivnorm::pbivnorm(as.vector(mu1draws[,i]),as.vector(mu2draws[,i]),rho = as.vector(pardraws[i,k1+k2+1])))
  if(is.vector(p11star))
    p11 = mean(p11star)
  else
    p11 = rowMeans(p11star)
  
  p10 = p1. - p11
  p01 = p.1 - p11
  p00 = 1-p11-p10-p01
  yHat1 = 1*(p1.>=mRule[1])
  yHat2 = 1*(p.1>=mRule[2])
  if(is.null(jRule)){
    ZHat = yHat1 * yHat2
  }else{
    ZHat = 1*(p11 >= jRule)
  }
  
  out = data.frame(linPredict1,linPredict2,p1.,p.1,p00,p10,p01,p11,yHat1,yHat2,ZHat)
  
  return(out)
  
}

#' @title predict method for class 'BiProbitPartialf'
#' 
#' @description Note, this is a simple 
#'  frequentist prediction and does not account for estimation uncertainty. If 
#'  one wants to account for estimation uncertainty it is reccomended to use the
#'  Bayesian philosophy.
#'
#' @param object a object of class \code{BiProbitPartialf}
#' @param newdata a matrix of column dimension k1 + k2 where the first k1 columns
#' correspond to the predictors of the first equations and the second k2 
#' columns correspond to predictors of the second equation. If intercepts were
#' used they need to be explicitly input.
#' @param mRule a vector of length 1 or 2. This is the marginal decision rule
#' for classifying the outcomes for stages 1 and 2. Stage 1 is 
#' classified as 1 if the probability of stage 1 being 1 is greater
#' than or equal to \code{mRule[1]}. Likewise for stage 2.  If 
#' length of \code{mRule} is 1 then that value is recycled. The values of 
#' \code{mRule} must be between 0 and 1. The default value is \code{mRule = c(0.5,0.5)}.
#' @param jRule an optional numerical value between 0 and 1. If specified 
#' then the observable outcome (both stages being 1) is 1 if the joint
#' probability of both stages being 1 is greater than jRule. If jRule
#' is unspecified or set to \code{NULL} then the observable outcome is the product
#' of the marginal outcomes. The default value is \code{jRule = NULL}. Note, if
#' \code{jRule} is specified then the observable outcome might not equal the product of stages 1 and 2.
#' @param ... unused
#'
#' @return method \code{predict.fBiProbitPArtial} returns a data.frame with columns
#' \describe{
#' 
#' \item{linPredict1}{Predicted mean of the first stage latent outcome}
#' 
#' \item{linPredict2}{Predicted mean of the second stage latent outcome}
#' 
#' \item{p1.}{Probability the outcome of the first stage is 1}
#' 
#' \item{p.1}{Probability the outcome of the second stage is 1}
#'
#' \item{p00}{Probability the outcome of both stages is 0}
#' 
#' \item{p01}{Probability the outcome of the first stage is 0 and the second stage is 1}
#' 
#' \item{p10}{Probability the outcome of stage 1 is 1 and stage 2 is 0}
#' 
#' \item{p11}{Probability the outcome of both stages are 1 }
#' 
#' \item{yHat1}{Classification of the outcome for stage 1. This value
#' is 1 if \code{p1 >= mRule[1]} and 0 else}
#' 
#' \item{yHat2}{Classification of the outcome for stage 2. This value
#' is 1 if \code{p2 >= mRule[2]} and 0 else}
#' 
#' \item{ZHat}{Classification of the observable outcome.
#' If \code{jRule} is specified then this value is 1 if \code{p12 >= jRule}
#'  and 0 else. If \code{jRule} is unspecified then this value is the element-wise product of
#'  yHat1 and yHat2.}
#' 
#' }
#' 
#' @examples 
#' ##
#' # Perform a prediction with the same covariates the model is estimated with
#' ##
#' 
#' data('Mroz87',package = 'sampleSelection')
#' Mroz87$Z = Mroz87$lfp*(Mroz87$wage >= 5)
#' 
#' f1 = BiProbitPartial(Z ~ educ + age + kids5 + kids618 + nwifeinc | educ + exper + city, 
#'      data = Mroz87, philosophy = "frequentist")
#' 
#' library(Formula)
#' eqn = Formula::Formula( ~ educ + age + kids5 + kids618 + nwifeinc | educ + exper + city)
#' matrix1 = model.matrix(eqn, lhs = 0, rhs=1, data= Mroz87)
#' matrix2 = model.matrix(eqn, lhs = 0, rhs=2, data= Mroz87)
#' newdat = cbind(matrix1,matrix2) 
#' preds1 = predict(f1,newdat)
#' head(preds1)
#' preds2 = predict(f1,newdat, jRule = .25)
#' 
#' # Compare predicted outcome with realized outcome
#' head(cbind(Mroz87$Z,preds1$ZHat,preds2$ZHat),20)

#' 
#' @export
#'
predict.BiProbitPartialf = function(object, newdata, mRule = c(.5,.5), jRule = NULL,...){
  
  newdata = as.matrix(newdata)
  if(dim(newdata)[2]==1)
    newdata = t(newdata)
  
  if(dim(newdata)[2] != object$k1 + object$k2)
    stop("newdata must be a matrix with column length k1 + k2")
  
  if(is.null(object$fixrhoVal))
    rho = utils::tail(object$par,1)
  else
    rho = object$fixrhoVal
  
  if(length(mRule)>2 | any(mRule>1) | any(mRule<0))
    stop("mRule must be a value between 0 and 1 of length 1 or 2 or unspecified")
  else if(length(mRule)==1)
    mRule = rep(mRule,2)
    
  
  linPredict1 = newdata[,1:object$k1] %*% object$par[1:object$k1]
  linPredict2 = newdata[,1:object$k2 + object$k1] %*% object$par[1:object$k2 + object$k1]
  
  p1. = stats::pnorm(0,linPredict1,lower.tail = FALSE)
  p.1 = stats::pnorm(0,linPredict2,lower.tail = FALSE)
  p11 = pbivnorm::pbivnorm(as.vector(linPredict1),as.vector(linPredict2),rho)
  p10 = p1. - p11
  p01 = p.1 - p11
  p00 = 1-p11-p10-p01
  yHat1 = 1*(p1.>=mRule[1])
  yHat2 = 1*(p.1>=mRule[2])
  if(is.null(jRule)){
    ZHat = yHat1 * yHat2
  }else{
    ZHat = 1*(p11 >= jRule)
  }
  
  out = data.frame(linPredict1,linPredict2,p1.,p.1,p00,p10,p01,p11,yHat1,yHat2,ZHat)
  
  return(out)
  
}

#' This is data to be included in my package
#'
#' @name SimDat
#' @docType data
#' @author Michael Guggisberg \email{mguggisb@@ida.org}
#' @keywords data
#' @description Simulated data of 10,000 observations from a multivariate
#' normal distribution. The true coefficients for equation 1 are 0 and 1. The
#' true coefficients for equation 2 are 0 and -1. The true \eqn{\rho} is 0.
NULL


#'@aliases BiProbitPartial-package
"_PACKAGE"
 


