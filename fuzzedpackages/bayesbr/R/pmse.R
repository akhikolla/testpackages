#'@title Prediction Mean Squared Error in a Beta Regression on a Bayesian Model
#'@name pmse
#'@aliases pmse
#'@description A function that selects a part of the database to fit a beta regression model and another part of this database to test the built model, returning the PMSE (prediction mean squared error) that reports the quality of the estimation for that database. In addition, the function also contains all the information that the bayesbr function returns, making it possible to do all analyzes on the fitted model.
#'@usage pmse(formula = NULL, data = NULL, test.set = 0.3,
#'na.action = c("exclude", "replace"),mean_betas = NULL,
#'    variance_betas = NULL,mean_gammas = NULL, variance_gammas = NULL,
#'     iter = 10000, warmup = iter/2,chains = 1, pars = NULL,
#'      a = NULL, b = NULL, resid.type = c("quantile",
#'      "sweighted", "pearson", "ordinary"), ...)
#'@param formula symbolic description of the model (of type \code{y ~ x} or \code{y ~ x | z};). See more at \code{\link{formula}}
#'@param data data frame or list with the variables passed in the formula parameter, if \code{data = NULL} the function will use the existing variables in the global environment.
#'@param test.set Defines the proportion of the database that will be used for testing the adjusted model and calculating the PMSE. The rest of the database will be used for modeling. Test.set must be less than "0.5", so that more than 50\% of the database is used to adjust the model.
#'@param na.action Characters provided or treatment used in NA values. If \code{na.action} is equal to exclude (default value), the row containing the NA will be excluded in all variables of the model. If \code{na.action} is equal to replace, the row containing the NA will be replaced by the average of the variable in all variables of the model.
#'@param mean_betas,variance_betas vectors including a priori information of mean and variance for the estimated beta respectively, beta is the name given to the coefficient of each covariate that influences theta. PS: the size of the vectors must equal p + 1, p being the number of covariates for theta.
#'@param mean_gammas,variance_gammas vectors including a priori information of mean and variance for the estimated ranges respectively, gamma is the name given to the coefficient of each covariate that influences zeta. PS: the size of the vectors must be equal to q + 1, q being the number of covariates for zeta.
#'@param iter A positive integer specifying the number of iterations for each chain (including warmup). The default is 10000.
#'@param warmup A positive integer specifying the number of iterations that will be in the warm-up period,
#'will soon be discarded when making the estimates and inferences. Warmup must be less than \code{iter} and its default value is \code{iter/2}.
#'@param chains A positive integer specifying the number of Markov chains. The default is 1.
#'@param pars A vector of character strings specifying parameters of interest. The default is NULL indicating all parameters in the model.
#'@param a,b Positive integer specifying the a priori information of the parameters of the gamma distribution for the zeta, if there are covariables explaining zeta \code{a} and \code{b} they will not be used.
#'@param resid.type A character containing the residual type returned by the model among the possibilities. The type of residue can be \emph{quantile}, \emph{sweighted}, \emph{pearson} or \emph{ordinary}. The default is \emph{quantile}.
#'@param ... 	Other optional parameters from RStan
#'@return \code{pmse} return an object of class \emph{pmse_bayesbr} containing the value of the prediction mean squared error and an object of the class \emph{bayesbr} with the following items:
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
#'\item{BIC}{a value containing the Bayesian Information Criterion (BIC) of the fitted model,}
#'\item{pseudo.r.squared}{pseudo-value of the square R (correlation to the square of the linear predictor and the a posteriori means of theta).}}
#'@references
#'  \doi{10.1080/0266476042000214501} Ferrari, S.L.P., and Cribari-Neto, F. (2004).
#'Beta Regression for Modeling Rates and Proportions. \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
#'@references
#'\doi{10.1016/j.jeconom.2005.07.014} Clark, T. E., & West, K. D. (2006). Using out-of-sample mean squared prediction errors to test the martingale difference hypothesis. \emph{Journal of econometrics}, \bold{135}(1-2), 155-186.
#'@seealso \code{\link{bayesbr}},\code{\link{residuals.bayesbr}},\code{\link{predict.bayesbr}}
#'@examples
#'data("bodyfat",package="bayesbr")
#'\dontshow{
#' lines = sample(1:251,15)
#' bodyfat = bodyfat[lines,]
#' }
#'bbr = pmse(siri ~ age + weight| biceps + forearm, data = bodyfat,
#'              test.set = 0.25, iter = 100)
#'
#'pmse = bbr$PMSE
#'model = bbr$model
#'summary(model)
#'residuals(model,type="sweighted")
#'\donttest{
#'bbr2 = pmse(siri ~ age + weight + height +
#'               wrist | biceps + forearm, data = bodyfat,
#'               test.set = 0.4, iter = 1000,
#'               mean_betas = 3,variance_betas = 10)

#'pmse2 = bbr2$PMSE
#'model2 = bbr2$model
#'residuals(model2,type="sweighted")
#'}
#'@export
pmse = function(formula=NULL,data=NULL,test.set=0.3,na.action=c("exclude","replace"),mean_betas = NULL,
                   variance_betas = NULL,mean_gammas = NULL,
                   variance_gammas = NULL ,iter = 10000,warmup = iter/2,
                   chains = 1,pars=NULL,a = NULL,b = NULL, resid.type = c("quantile","sweighted",
                                                                          "pearson","ordinary"),...){
  cl = match.call()
  if(test.set<0.5){
    if(!is.null(formula)){
      dados = formula(as.formula(formula),data)
      Y = dados[[1]]
      X = dados[[2]]
      W = dados[[3]]
      name_y  = dados[[4]]
      names_x = dados[[5]]
      names_w = dados[[6]]
    }
    else{
      stop("formula not informed")
    }
    if(is.null(data)){
      data = model.bayesbr(Y,X,W,name_y,names_x,names_w)
    }

    if(!is.null(X)){
      if(is.matrix(X)){
        p = ncol(X)
      }
      else{
        p=1
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
    }
    else{
      q=0
    }

    na.action = match.arg(na.action)
    if(na.action == "exclude"){
      model = data.frame(cbind(Y,X,W))
      na_values = which(is.na(model), arr.ind=TRUE)
      if(nrow(na_values)>0){
        model =drop_na(model)
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


    tam = nrow(data)
    data_model = sample_frac(data,1-test.set)
    data_analyze = suppressMessages(anti_join(data, data_model))

    bayesbr = bayesbr(formula = formula,data = data_model,na.action = na.action,
                      mean_betas =  mean_betas,variance_betas = variance_betas,
                      mean_gammas = mean_gammas,variance_gammas = variance_gammas,
                      iter = iter,warmup = warmup,chains = chains,pars = pars,
                      a = a,b = b,resid.type = resid.type,...)

    predict = as.vector(predict(bayesbr,newdata = data_analyze,type="response"))
    Y = as.vector(as.matrix(data_analyze[,name_y]))
    PMSE = mean((Y-predict)^2)
    bayesbr$PMSE = round(PMSE,5)
    bayesbr$call = cl
    rval = list()
    rval$PMSE = round(PMSE,5)
    rval$model = bayesbr
    class(rval) = "pmse_bayesbr"
    return(rval)
  }
  else{
    stop("The percentage of the data frame to be adjusted must be greater than 50%.")
  }
}
#'@export
print.pmse_bayesbr = function(x,...){
  print(x$PMSE)
}
