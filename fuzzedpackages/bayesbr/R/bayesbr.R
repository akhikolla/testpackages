#'@title Bayesian Beta Regression with RStan
#'@aliases bayesbr
#'@name bayesbr
#'@usage bayesbr(formula = NULL, data = NULL,na.action=c("exclude",
#'"replace"),mean_betas = NULL,variance_betas = NULL,mean_gammas = NULL,
#' variance_gammas = NULL,iter = 10000, warmup = iter/2,
#' chains = 1, pars = NULL, a = NULL, b = NULL,
#' resid.type = c("quantile","sweighted", "pearson","ordinary"), ...)
#'@param formula symbolic description of the model (of type \code{y ~ x} or \code{y ~ x | z};). See more at \code{\link{formula}}
#'@param data data frame or list with the variables passed in the formula parameter, if \code{data = NULL} the function will use the existing variables in the global environment.
#'@param na.action Characters provided or treatment used in NA values. If \code{na.action} is equal to exclude (default value), the row containing the NA will be excluded in all variables of the model. If \code{na.action} is equal to replace, the row containing the NA will be replaced by the average of the variable in all variables of the model.
#'@param mean_betas,variance_betas vectors including a priori information of mean and variance for the estimated beta respectively, beta is the name given to the coefficient of each covariate that influences theta. PS: the size of the vectors must equal p + 1, p being the number of covariates for theta.
#'@param mean_gammas,variance_gammas vectors including a priori information of mean and variance for the estimated ranges respectively, gamma is the name given to the coefficient of each covariate that influences zeta. PS: the size of the vectors must be equal to q + 1, q being the number of covariates for zeta.
#'@param iter A positive integer specifying the number of iterations for each chain (including warmup). The default is 10000.
#'@param warmup A positive integer specifying the number of iterations that will be in the warm-up period,
#'will soon be discarded when making the estimates and inferences. Warmup must be less than \code{iter} and its default value is \code{iter/2}.
#'@param chains A positive integer specifying the number of Markov chains. The default is 1.
#'@param pars A vector of character strings specifying parameters of interest. The default is NULL indicating all parameters in the model.
#'@param a,b Positive integer specifying the a priori information of the parameters of the gamma distribution for the zeta, if there are covariables explaining zeta \code{a} and \code{b} they will not be used. The default value for \code{a} is 1 and default value for \code{b} is 0.01 .
#'@param resid.type A character containing the residual type returned by the model among the possibilities. The type of residue can be \emph{quantile}, \emph{sweighted}, \emph{pearson} or \emph{ordinary}. The default is \emph{quantile}.
#'@param ... 	Other optional parameters from RStan
#'@return   \code{bayesbr} return an object of class \emph{bayesbr}, a list of the following items.
#'\describe{\item{coefficients}{a list with the mean and precision elements containing the estimated coefficients of model and table with the means, medians, standard deviations and the Highest Posterior Density (HPD) Interval,}
#'\item{call}{the original function call,}
#'\item{formula}{the original formula,}
#'\item{y}{the response proportion vector,}
#'\item{stancode}{lines of code containing the .STAN file used to estimate the model,}
#'\item{info}{a list containing model information such as the argument pars passed as argument, name of variables, number of: iterations, warmups, chains, covariables for theta, covariables for zeta and observations of the sample. In addition there is an element called samples, with the posterior distribution of the parameters of interest,}
#'\item{fitted.values}{a vector containing the estimates for the values corresponding to the theta of each observation of the variable response, the estimate is made using the mean of the a prior theta distribution,}
#'\item{model}{the full model frame,}
#'\item{residuals}{a vector of residuals}
#'\item{residuals.type}{the type of returned residual,}
#'\item{loglik}{log-likelihood of the fitted model(using the mean of the parameters in the posterior distribution),}
#'\item{AIC}{a value containing the Akaike's Information Criterion (AIC) of the fitted model,}
#'\item{BIC}{a value containing the Bayesian Information Criterion (BIC) of the fitted model,}
#'\item{DIC}{a value containing the Deviance Information Criterion (DIC) of the fitted model,}
#'\item{WAIC}{a vector containing the Widely Applicable Information Criterion (WAIC) of the fitted model and their standard error, see more in \code{\link{waic}}}
#'\item{LOOIC}{a vector containing the LOO (Efficient approximate leave-one-out cross-validation) Information Criterion of the fitted model and their standard error, see more in \code{\link{loo}}}
#'\item{pseudo.r.squared}{pseudo-value of the square R (correlation to the square of the linear predictor and the a posteriori means of theta).}}
#'@description Fit of beta regression model under the view of Bayesian statistics,
#'using the mean of the posterior distribution as estimates for the mean (theta) and the precision parameter (zeta).
#'@seealso \code{\link{summary.bayesbr}}, \code{\link{residuals.bayesbr}}, \code{\link{formula}}
#'@references
#'  \doi{10.1080/0266476042000214501} Ferrari, S.L.P., and Cribari-Neto, F. (2004).
#'Beta Regression for Modeling Rates and Proportions. \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
#'@references
#'\href{https://arxiv.org/abs/1111.4246}{arXiv:1111.4246} Hoffman, M. D., and Gelman, A. (2014). The No-U-Turn sampler: adaptively setting path lengths in Hamiltonian Monte Carlo. \emph{Journal of Machine Learning Research}, \bold{15}(1), 1593-1623.
#'@references
#'\doi{10.18637/jss.v076.i01} Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B., Betancourt, M., ... & Riddell, A. (2017). Stan: A probabilistic programming language. \emph{Journal of statistical software}, \bold{76}(1).
#'@examples
#'data("StressAnxiety",package="bayesbr")
#'
#'bbr = bayesbr(anxiety ~ stress | stress, data = StressAnxiety,
#'              iter = 100)
#'summary(bbr)
#'residuals(bbr, type="ordinary")
#'print(bbr)
#'
#'
#'\donttest{
#'data("StressAnxiety", package = "bayesbr")
#'bbr2 <- bayesbr(anxiety ~ stress | stress,
#'                data = StressAnxiety, iter = 1000,
#'                warmup= 450, mean_betas = c(0,1),
#'                variance_betas = 15)
#'
#'envelope(bbr2,sim=100,conf=0.95)
#'loglikPlot(bbr2$loglik)
#'}
#'@details   Beta Regression was suggested by Ferrari and Cribari-Neto (2004), but with the look of classical statistics, this package makes use of the \code{Rstan} to, from the prior distribution of the data, obtain the posterior distribution and the estimates from a Bayesian perspective. Beta regression is useful when the response variable is in the range between 0 and 1, being used for adjusting probabilities and proportions.
#'
#'It is possible to estimate coefficients for the explanatory covariates for the theta and zeta parameters of the Beta distribution. Linear predictors are passed as parameters for both zeta and zeta, from these linear predictors a transformation of scales is made.
#'
#'Hamiltonian Monte Carlo (HMC) is a Markov chain Monte Carlo (MCMC) algorithm, from the HMC there is an extension known as the No-U-Turn Sampler (NUTS) that makes use of recursion to obtain its calculations and is used by \code{RStan}. In the context of the \code{bayesbr} package, NUTS was used to obtain a posteriori distribution from model data and a priori distribution.
#'
#'See \code{\link{predict.bayesbr}}, \code{\link{residuals.bayesbr}},\code{\link{summary.bayesbr}},\code{\link{logLik.bayesbr}} and \code{\link{pseudo.r.squared}} for more details on all methods. Because it is in the context of Bayesian statistics, in all calculations that were defined using maximum verisimilitude, this was sub-replaced by the mean of the posterior distribution of the parameters of interest of the formula.
#'@export
bayesbr = function(formula=NULL,data=NULL,na.action=c("exclude","replace"),mean_betas = NULL,
                    variance_betas = NULL,mean_gammas = NULL,
                    variance_gammas = NULL ,iter = 10000,warmup = iter/2,
                    chains = 1,pars=NULL,a = NULL,b = NULL, resid.type = c("quantile","sweighted", "pearson","ordinary"),...){
  cl = match.call()
  r_mc_aux = T
  resid.type = match.arg(resid.type)
  if(!is.null(formula)){
    dados = formula(as.formula(formula),data)
    Y = dados[[1]]
    X = dados[[2]]
    W = dados[[3]]
    name_y  = dados[[4]]
    names_x = dados[[5]]
    names_w = dados[[6]]
    model = data.frame(cbind(Y,X,W))
    char = FALSE
    for (cnames in colnames(model)){
          if(!is.numeric(model[,cnames])){
            char = TRUE
            break;
          }
    }
    if(char == TRUE){
      stop("It is not possible to work with categorical variables, to use this variable you must qualitative variables of interest in quantitative variables.",call.=TRUE)
    }
  }
  else{
    stop("formula not informed",call.=TRUE)
  }
  if(!is.null(a)){
    if(a<0){
      stop("zeta 'a' priori parameter cannot be negative",call.=TRUE)
    }
  }
  if(!is.null(b)){
    if(b<0){
      stop("zeta 'b' priori parameter cannot be negative",call.=TRUE)
    }
  }
  if(is.numeric(warmup)){
    warmup = as.integer(warmup)
  }

  if(!is.null(X)){
    if(is.matrix(X)){
      p = ncol(X)
    }
    else{
      p=1
    }
    if(is.null(mean_betas)){
      mean_betas = rep(0,p)
    }
    if(is.null(variance_betas)){
      variance_betas = rep(10,p)
    }
  }
  else{
    p = 0
  }
  if(!is.null(W)){
    if(is.matrix(W)){
      q = ncol(W)
    }
    else{
      q=1
    }
    if(is.null(mean_gammas)){
      mean_gammas = c(rep(0,q))
    }
    if(is.null(variance_gammas)){
      variance_gammas = c(rep(10,q))
    }
  }
  else{
    q=0
  }

  na.action = match.arg(na.action)
  if(na.action == "exclude"){
    model = data.frame(cbind(Y,X,W))
    na_values = which(is.na(model), arr.ind=TRUE)
    if(nrow(na_values)>0){
      model = drop_na(model)
      Y = model[,1]
      aux1 = 1+p
      aux2 = 2+p
      aux3 = 1+p+q
      if(p>0){
      X = model[,(2:aux1)]
      }
      if(q>0){
        W = model[,(aux2:aux3)]
      }
      warning("The model variables may have changed, for more details check the complete model returned in the item model."
                ,call.=TRUE)
  }
  }
  if(na.action == "replace"){
    model = cbind(Y,X,W)
    na_values = which(is.na(model), arr.ind=TRUE)
    if(nrow(na_values)>0){
      for(i in 1:nrow(na_values)){
        row = na_values[i,1]
        col = na_values[i,2]
        mean = mean(model[,col],na.rm=TRUE)
        model[row,col] = mean
      }
      Y = model[,1]
      aux1 = 1+p
      aux2 = 2+p
      aux3 = 1+p+q
      if(p>0){
        X = model[,(2:aux1)]
      }
      if(q>0){
        W = model[,(aux2:aux3)]
      }
      warning("The model variables may have changed, for more details check the complete model returned in the item model.")
      }
  }
  n = length(Y)
  if(max(Y)>=1 || min(Y)<=0){
    warning("Some of your data is outside the range between 0 and 1 (extremes not included in the range), so beta regression cannot be applied.")
  }
  if(length(mean_betas)==1 && p>1){
    mean_betas = rep(mean_betas,p)
  }
  if(length(variance_betas)==1 && p>1){
    variance_betas = rep(variance_betas,p)
  }
  if(length(mean_gammas)==1 && q>1){
    mean_gammas = rep(mean_gammas,q)
  }
  if(length(variance_gammas)==1 && q>1){
    variance_gammas = rep(variance_gammas,q)
  }
  if(length(mean_betas)!=p){
    stop("The number of a priori specifications for betas averages must equal the p + 1 (p being the number of covariates for theta)",call.=TRUE)
  }
  if(length(variance_betas)!=p){
    stop("The number of a priori specifications for betas variances must equal the p + 1 (p being the number of covariates for theta)",call.=TRUE)
  }
  if(length(mean_gammas)!=q){
    stop("The number of a priori specifications for gammas averages must equal the q + 1 (q being the number of covariates for zeta)",call.=TRUE)
  }
  if(length(variance_gammas)!=q){
    stop("The number of a priori specifications for gammas variances must equal the q + 1 (q being the number of covariates for zeta)",call.=TRUE)
  }

  if(is.null(a) || is.null(b)){
    a = 1
    b = 0.01
  }
  if(p==0){
    a = 2
    b = 2
  }
  data = list(n=n, p = p, q = q, Y=Y,a=a,b=b)
  pars_aux = c()
  if(p==0){
    data$X = matrix(1,n,0)
    data$mean_betas = vector("double",0)
    data$variance_betas = vector("double",0)
    pars_aux = c(pars_aux,"theta_e")
  }
  else{
    data$X = X
    data$mean_betas = array(mean_betas)
    data$variance_betas = array(variance_betas)
    pars_aux = c("betas","theta",pars_aux)
  }
  if(q==0){
    data$W = matrix(1,n,0)
    data$mean_gammas = vector("double",0)
    data$variance_gammas = vector("double",0)
    pars_aux = c(pars_aux,"zeta_e")
  }
  else{
    data$W = W
    data$mean_gammas = array(mean_gammas)
    data$variance_gammas = array(variance_gammas)
    pars_aux = c("gammas","zeta",pars_aux)
  }
  if(is.null(pars)){
    pars = c("betas","zeta_e","theta","theta_e","gammas","zeta")
  }

  if(!("betas" %in% pars) && p>0){
    warning('"betas" has to be included in the pars argument, so that the coefficients are calculated',
            call.= F)
  }
  if(!("gammas" %in% pars) && q>0){
    warning('"gammas" has to be included in the pars argument, so that the coefficients are calculated',
            call.= F)
  }


  object = sampling(stanmodels$bayesbr, data=data,
                    iter=iter, warmup=warmup,pars = pars_aux, chains=chains, ...)

  if(p==0){
    object@sim$samples[[1]]$theta = object@sim$samples[[1]]$theta_e
  }
  if(q==0){
    object@sim$samples[[1]]$zeta = object@sim$samples[[1]]$zeta_e
  }
  betas = values("beta",object,iter,warmup,n,p)
  gammas = values("gamma",object,iter,warmup,n,q)
  theta = values("theta",object,iter,warmup,n,p)
  zeta = values("zeta",object,iter,warmup,n,q)
  model = model.bayesbr(Y,X,W,name_y,names_x,names_w)
  names(Y) = 1:n

  rval = list()
  class(rval) = "bayesbr"

  rval$call = cl
  rval$formula = formula
  rval$y = Y
  rval$stancode = object@stanmodel
  rval$info = list(n = n, iter = iter, warmup = warmup, chains = chains, p = p, q = q)
  rval$info$names = list(name_y=name_y,names_x = names_x, names_w = names_w)
  rval$info$samples$beta = betas
  rval$info$samples$gamma = gammas
  rval$info$samples$theta = theta
  rval$fitted.values = fitted.values(rval)
  rval$info$samples$zeta = zeta
  rval$pars = pars_aux
  rval$model = model
  if(p>0){
    rval$residuals.type = resid.type
    res = residuals.bayesbr(rval,rval$residuals.type)
    rval$residuals = res
  }
  list_mean = summary_mean(rval)
  list_precision = summary_precision(rval)

  rval$coefficients = list(mean = list_mean[['betas']],
                           precision = list_precision[['gammas']],
                           summary_mean = list_mean[['table']],
                           summary_precision = list_precision[['table']])
  if(p>0){
    list_loglik = logLik.bayesbr(rval)
    rval$loglik = list_loglik$loglik
    waic_estimates = suppressWarnings(waic(list_loglik$matrix_loglik)$estimates)
    looic_estimates = suppressWarnings(loo(list_loglik$matrix_loglik)$estimates)
    rval$AIC   = AIC_bayesbr(rval)
    rval$BIC   = BIC_bayesbr(rval)
    rval$DIC   = DIC_bayesbr(rval)
    rval$WAIC  = waic_estimates[3,]
    rval$LOOIC = looic_estimates[3,]
    rval$pseudo.r.squared. = pseudo.r.squared(rval)
  }
  if(!("theta" %in% pars)){
    rval$fitted.values = NULL
    rval$info$samples$theta = NULL
  }
  if(!("zeta" %in% pars)){
    rval$info$samples$zeta = NULL
  }
  return(rval)
}
