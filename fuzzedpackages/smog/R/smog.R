#' Generalized linear model constraint on hierarchical structure
#' by using overlapped group penalty
#' 
#' \code{smog} fits a linear non-penalized phynotype (demographic) variables such as 
#' age, gender, treatment, etc, and penalized groups of prognostic effect (main effect)
#' and predictive effect (interaction effect), by satisfying the hierarchy structure:
#' if a predictive effect exists, its prognostic effect must be in the model. It can deal
#' with continuous, binomial or multinomial, and survival response variables, underlying 
#' the assumption of Gaussian, binomial (multinomial), and Cox proportional hazard models,
#' respectively. It can accept \code{\link[stats]{formula}}, and output coefficients table,
#' fitted.values, and convergence information produced in the algorithm iterations.   
#' 
#' @param x a model matrix, or a data frame of dimensions n by p, 
#'          in which the columns represents the predictor variables. 
#' @param y response variable, corresponds to the family description. 
#'          When family is ''gaussian'' or ''binomial'', \code{y} ought to
#'          be a numeric vector of observations of length n; when family 
#'          is ''coxph'', \code{y} represents the survival objects, containing the 
#'          survival time and the censoring status. See \code{\link[survival]{Surv}}.
#' @param g a vector of group labels for the predictor variables.
#' @param v a vector of binary values, represents whether or not the 
#'          predictor variables are penalized. Note that 1 indicates 
#'          penalization and 0 for not penalization.
#' @param label a character vector, represents the type of predictors in terms of treatment,
#'              prognostic, and predictive effects by using ''t'', ''prog'', and ''pred'',
#'              respectively.
#' @param lambda1 penalty parameter for the L2 norm of each group of prognostic and predictive effects.
#' @param lambda2 ridge penalty parameter for the squared L2 norm of each group of prognostic and predictive effects.
#' @param lambda3 penalty parameter for the L1 norm of predictive effects.              
#' @param family a description of the distribution family for the response 
#'               variable variable. For continuous response variable,
#'               family is ''gaussian''; for multinomial or binary response
#'               variable, family is ''binomial''; for survival response
#'               variable, family is ''coxph'', respectively.
#' @param subset an optional vector specifying a subset of observations to be 
#'               used in the model fitting. Default is \code{NULL}.
#' @param rho   the penalty parameter used in the alternating direction method 
#'              of multipliers (ADMM) algorithm. Default is 10.
#' @param scale whether or not scale the design matrix. Default is \code{TRUE}.
#' @param eabs  the absolute tolerance used in the ADMM algorithm. Default is 1e-3.
#' @param erel  the reletive tolerance used in the ADMM algorithm. Default is 1e-3.
#' @param LL    initial value for the Lipschitz continuous constant for 
#'              approximation to the objective function in the Majorization-
#'              Minimization (MM) (or iterative shrinkage-thresholding algorithm 
#'              (ISTA)). Default is 1.
#' @param eta   gradient stepsize for the backtrack line search for the Lipschitz
#'              continuous constant. Default is 1.25. 
#' @param maxitr the maximum iterations for convergence in the ADMM algorithm. 
#'               Default is 1000.
#' @param formula an object of class ''formula'': a symbolic description of the
#'                model to be fitted. Should not include the intercept. 
#' @param data    an optional data frame, containing the variables in the model. 
#' @param ...   other relevant arguments that can be supplied to smog.
#' 
#' @return \code{smog} returns an object of class inhering from ''smog''. The 
#'         generic accessor functions \code{coef}, \code{coefficients}, 
#'         \code{fitted.value}, and \code{predict} can be used to extract
#'         various useful features of the value returned by \code{smog}.
#'         An object of ''smog'' is a list containing at least the following 
#'         components: 
#'         \subsection{}{
#'         \describe{
#'         \item{coefficients}{Data frame containing the nonzero predictor
#'                             variables' indexes, names, and estimates. When
#'                             family is ''binomial'', the estimates have K-1 
#'                             columns, each column representing the weights for the 
#'                             corresponding group. The last group behaves the
#'                             ''pivot''.}
#'         \item{fitted.values}{The fitted mean values for the response variable,
#'                              for family is ''gaussian''. When family is 
#'                              ''binomial", the fitted.values are the probabilies
#'                              for each class; when family is ''coxph'', 
#'                              the fitted.values are risk scores.}
#'         \item{residuals}{The residual is trivial for family = "gaussian".
#'                          For family = "binomial", Pearson residuals is returned; and
#'                          for family = "coxph", it yields deviance residuals, i.e., 
#'                          standardized martingale residuals.}
#'         \item{model}{A list of estimates for the intercept, treatment effect, 
#'                      and prognostic and predictive effects for the selectd
#'                      biomarkers.}
#'         \item{weight}{The weight of predictors resulted from the penalty funciton,
#'                       is used to calculate the degrees of freedom.}
#'         \item{DF}{the degrees of freedom. When family = ''gaussian'', 
#'                   \eqn{DF = tr(x_{\lambda}'(x_{\lambda}'x_{\lambda}+W)x_{\lambda})}. 
#'                   For other families, DF is approximated by \eqn{diag(1/(1+W))}.}
#'         \item{criteria}{model selection criteria, including the correction Akaike's Information 
#'                         Criterion (AIC), AIC, Bayesian Information Criterion (BIC), and the generalized 
#'                         cross-validation score (GCV), respectively. See also \code{\link{cv.smog}}.}
#'         \item{llikelihood}{the log-likelihood value for the converged model.}
#'         \item{loglike}{the penalized log-likelihood values for each 
#'                        iteration in the algorithm.}
#'         \item{PrimalError}{the averged norms \eqn{||\beta-Z||/\sqrt{p}} for each iteration,
#'                            in the ADMM algorithm.}
#'         \item{DualError}{the averaged norms \eqn{||Z^{t+1}-Z^{t}||/\sqrt{p}} for 
#'                          each iteration, in the ADMM algorithm.}
#'         \item{converge}{the number of iterations processed in the ADMM algorithm.}
#'         \item{call}{the matched call.}
#'         \item{formula}{the formula supplied.}
#'         }
#'         }
#' 
#' @details The formula has the form \code{response ~ 0 + terms} where \code{terms} is
#'          a series of predictor variables to be fitted for \code{response}. For \code{gaussian} 
#'          family, the response is a continuous vector. For \code{binomial} family, 
#'          the response is a factor vector, in which the last level denotes the ''pivot''.
#'          For \code{coxph} family, the response is a \code{\link[survival]{Surv}} 
#'          object, containing the survival time and censoring status.
#' 
#' @section Penalized regression model: The regression function contains the non-penalized predictor variables, 
#'  and many groups of prognostic and predictive terms, where in each group the prognostic term comes first, 
#'  followed by the predictive term.
#'          
#'  * `Penalty function`: Different hierachical structures within groups can result from adjusting 
#'    the penalty parameters in the penalty function:        
#'      
#'      \deqn{\Omega(\mathbf{\beta}) = \lambda_1||\mathbf{\beta}|| + \lambda_2||\mathbf{\beta}||^2+\lambda_3|\beta_2|}
#'    
#'    Where \eqn{\mathbf{\beta}=(\beta_1,\beta_2)}. Note that \eqn{\beta_1} denotes the prognostic effect 
#'    (main effect), and \eqn{\beta_2} for the predictive effect (interactive effect), respectively. 
#'    When \eqn{\lambda_2 = 0} and \eqn{\lambda_3 = 0}, it indicates no structure within groups. When 
#'    \eqn{\lambda_2 \ne 0}, the penalty function honors the structure within groups such that: 
#'    predictive effect \eqn{\ne 0 \Longrightarrow} prognostic effect \eqn{\ne 0}.
#'         
#'  * `Tuning parameters`: \code{rho,eabs,erel,LL,eta} are the corresponding parameters used in the 
#'    itervative shrinkage-thresholding algorithm (ISTA) and the alternating direction method of 
#'    multipliers algorithm (ADMM).
#'
#'  
#' @references \insertRef{ma2019structural}{smog}
#'  
#' @examples  
#' 
#' n=100;p=20
#' set.seed(2018)
#' # generate design matrix x
#' s=10
#' x=matrix(0,n,1+2*p)
#' x[,1]=sample(c(0,1),n,replace = TRUE)
#' x[,seq(2,1+2*p,2)]=matrix(rnorm(n*p),n,p)
#' x[,seq(3,1+2*p,2)]=x[,seq(2,1+2*p,2)]*x[,1]
#' 
#' g=c(p+1,rep(1:p,rep(2,p)))  # groups 
#' v=c(0,rep(1,2*p))           # penalization status
#' label=c("t",rep(c("prog","pred"),p))  # type of predictor variables
#' 
#' # generate beta
#' beta=c(rnorm(13,0,2),rep(0,ncol(x)-13))
#' beta[c(2,4,7,9)]=0
#' 
#' # generate y
#' data1=x%*%beta
#' noise1=rnorm(n)
#' snr1=as.numeric(sqrt(var(data1)/(s*var(noise1))))
#' y1=data1+snr1*noise1
#' lfit1=smog(x,y1,g,v,label,lambda1=8,lambda2=0,lambda3=8,family = "gaussian")
#' 
#' ## generate binomial data
#' prob=exp(as.matrix(x)%*%as.matrix(beta))/(1+exp(as.matrix(x)%*%as.matrix(beta)))
#' y2=ifelse(prob<0.5,0,1)
#' lfit2=smog(x,y2,g,v,label,lambda1=0.03,lambda2=0,lambda3=0.03,family = "binomial")
#' 
#' ## generate survival data
#' # Weibull latent event times
#' lambda = 0.01; rho = 1
#' V = runif(n)
#' Tlat = (- log(V) / (lambda*exp(x %*% beta)) )^(1/rho)
#' C = rexp(n, 0.001)  ## censoring time
#' time = as.vector(pmin(Tlat, C))
#' status = as.numeric(Tlat <= C)
#' y3 = as.matrix(cbind(time = time, status = status))
#' 
#' lfit3=smog(x,y3,g,v,label,lambda1=0.2,lambda2=0,lambda3=0.2,family = "coxph")
#' 
#' @importFrom tidyr spread
#' @importFrom magrittr %>% %<>%
#' @importFrom stats model.matrix
#' @import dplyr
#' 
#' @export
smog.default <- function(x, y, g, v, label, lambda1, lambda2, lambda3, family = "gaussian", subset = NULL, rho = 10, 
                         scale = TRUE, eabs = 1e-3, erel = 1e-3, LL = 1, eta = 1.25, maxitr = 1000, ...){
  
  lambda=c(lambda1,lambda2,lambda3)
  hierarchy = ifelse(lambda2 == 0, 
                     ifelse(lambda3 == 0, 0, 1),2)
  
  if(!is.null(subset)){
    x <- as.matrix(as.matrix(x)[subset,])
    y <- as.matrix(as.matrix(y)[subset,])
  }else{
    x <- as.matrix(x)
    y <- as.matrix(y)
  }
  
  g <- as.numeric(as.factor(g))
  est <- glog(y,x,g,v,lambda,hierarchy,family,rho,
              scale,eabs,erel,LL,eta,maxitr)
  
  if(nrow(est$coefficients)){  # continue for some variables are selected  
    if(family == "gaussian"){
      wx <- cbind(rep(1,nrow(x)),x)
      est$fitted.value = as.vector(as.matrix(wx[,est$coefficients$Id+1])%*%
                                     as.matrix(est$coefficients$Estimate))
      est$residuals = as.vector(y - est$fitted.value)
      
      if(!is.null(colnames(x))){
        est$coefficients$Beta = c("Intercept",colnames(x))[est$coefficients$Id+1]
        est$coefficients = est$coefficients[,c("Id","Beta","Estimate")]
      }
    }
    
    if(family == "binomial"){
      probTAB = exp(as.matrix(x[,est$coefficients$Id])%*%as.matrix(est$coefficients$Estimate))/
        (1+exp(as.matrix(x[,est$coefficients$Id])%*%as.matrix(est$coefficients$Estimate)))
      probTAB = cbind(1-rowSums(as.matrix(probTAB)),probTAB)
      predClass = apply(probTAB,1,which.max)
      predProb = apply(probTAB,1,max)
      
      est$levels = sort(unique(y))
      est$fitted.value = data.frame(Class = est$levels[predClass],
                                    Prob = predProb)
      
      # calculate the pearson residuals 
      y_mat = model.matrix( ~ 0 + as.factor(y)) %>% data.frame %>% as.matrix
      if(length(est$levels) > 2){
        est$residuals = rowSums((y_mat[,-1] - probTAB[,-1])/sqrt(probTAB[,-1]*(1 - probTAB[,-1])))
      } else{
        est$residuals = (y_mat[,-1] - probTAB[,-1])/sqrt(probTAB[,-1]*(1 - probTAB[,-1]))
      }
      
      
      if(!is.null(colnames(x))){
        est$coefficients$Beta = colnames(x)[est$coefficients$Id]
        est$coefficients = est$coefficients[,c(1,ncol(probTAB)+1,2:ncol(probTAB))]
      }
    }
    
    if(family == "coxph"){
      
      risk = as.vector(exp(as.matrix(x[,est$coefficients$Id]) %*% as.matrix(est$coefficients$Estimate)))
      
      # calculate baseline hazards
      y %<>% data.frame
      y_names = colnames(y)
      
      # sort the survival time
      y_sum = cbind(y, risk) %>%
        dplyr::arrange(get(y_names[1])) %>%
        dplyr::filter(get(y_names[2]) != 0) %>%
        dplyr::mutate(cumrisk = rev(cumsum(rev(risk))))
      
      y_expected = sapply(y[,1], function(z) sum(1/y_sum$cumrisk[y_sum[,1] <= z]))*risk
      martingale = y[,2] - y_expected
      est$residuals = sign(martingale)*sqrt(-2*(martingale + y[,2]*log(y[,2] - martingale)))
      est$fitted = risk
      
      if(!is.null(colnames(x))){
        est$coefficients$Beta = colnames(x)[est$coefficients$Id]
        est$coefficients = est$coefficients[,c("Id","Beta","Estimate")]
      }
    }
  }
  
  # calculate the degrees of freedom 
  idx = est$coefficients$Id
  est$model$intercept = ifelse(0 %in% idx,
                               est$coefficients$Estimate[idx==0], NA)
  est$model$treatment = ifelse(1 %in% idx,
                               est$coefficients$Estimate[idx==1], NA)
  GId = estimate = elabel = NULL
  if(any(idx>1)){
    est$model$biomarker = data.frame(GId = g[idx[idx>1]],
                                     estimate = est$coefficients$Estimate[idx>1],
                                     elabel = label[idx[idx>1]])
    est$model$biomarker = tidyr::spread(est$model$biomarker,elabel,estimate)
    est$model$biomarker = as.data.frame(apply(est$model$biomarker,1:2,
                                              function(t) ifelse(is.na(t),0,t)))
    
    if("prog" %in% colnames(est$model$biomarker) & 
       (!("pred" %in% colnames(est$model$biomarker)))){
      est$model$biomarker$pred = 0
    }
    
    if("pred" %in% colnames(est$model$biomarker) & 
       (!("prog" %in% colnames(est$model$biomarker)))){
      est$model$biomarker$prog = 0
    }
    
    est$weight = data.frame(t(apply(est$model$biomarker[,c("prog","pred")],1,
                                    function(t) lambda[1]/sqrt(sum(t^2))+
                                      c(0,ifelse(sign(t[2]),lambda[3]/abs(t[2]),0)))))
    colnames(est$weight) = c("prog","pred")
    
    if(!is.null(colnames(x))){
      est$model$biomarker$marker = colnames(x)[label=="prog" & (1:ncol(x) %in% idx[idx>1])]
      est$model$biomarker = est$model$biomarker[,c("GId","marker","prog","pred")]
      est$weight = cbind(est$model$biomarker[,c("GId","marker")],est$weight)
    }else{
      est$model$biomarker = est$model$biomarker[,c("GId","prog","pred")]
      est$weight = cbind(GId = est$model$biomarker$GId,est$weight)
    }
    
    Weight = sapply(idx[idx >= 1],function(t) ifelse(t == 1,0,
                                                     est$weight[est$weight$GId == g[t],label[t]]))
    if(family == "gaussian"){
      est$DF = sum(diag(solve(as.matrix(t(x[,idx[idx>=1]]))%*%as.matrix(x[,idx[idx >= 1]])+diag(Weight))
                        %*%as.matrix(t(x[,idx[idx>=1]]))%*%as.matrix(x[,idx[idx >= 1]])))
    }else{
      est$DF = sum(1/(1+Weight))
    }
    
  }else{
    if(!is.null(colnames(x))){
      est$model$biomarker = data.frame(GId = NA, marker = NA, prog = NA, pred = NA)
    }else{
      est$model$biomarker=data.frame(GId = NA, prog = NA, pred = NA)
    }
    
    est$weight = NULL
    est$DF = 1
  }
  
  # AIC, BIC, GCV criteria by using whole data 
  n <- nrow(x)
  k = est$DF
  
  aic1 = tryCatch({n/2*log(abs(2*est$llikelihood)) + n/2*((1+k/n)/(1-k+2/n))},
                  error=function(e) NA)
  aic2 = tryCatch({log(abs(2*est$llikelihood)/n) + 2*k/n}, error=function(e) NA)
  bic = tryCatch({log(abs(2*est$llikelihood)/n) + log(n)*k/n}, error=function(e) NA)
  gcv = tryCatch({abs(2*est$llikelihood)/(n*(1-k/n)^2)}, error=function(e) NA)
  
  
  est$criteria = list(aic1,aic2,bic,gcv)
  names(est$criteria) = c("cAIC","AIC","BIC","GCV")
  
  
  
  est$call <- match.call()
  class(est) <- "smog"
  est
}

#' @rdname smog.default
#' 
#' @seealso \code{\link{cv.smog}}, \code{\link{predict.smog}}, \code{\link{plot.smog}}.
#' @author Chong Ma, \email{chongma8903@@gmail.com}.
#' 
#' @importFrom stats model.frame model.matrix model.response
#' 
#' @export
smog.formula <- function(formula, data=list(), g, v, label, lambda1, lambda2, lambda3, ...){
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf,"terms"),data = mf)
  y <- model.response(mf)
  
  est <- smog.default(x,y,g,v,label,lambda1,lambda2,lambda3,...)
  est$call <- match.call()
  est$formula <- formula
  est
}


#' predict method for objects of the class smog
#' 
#' \code{predict.smog} can produce the prediction for user-given new data, based on the
#' provided fitted model (\code{object}) in the S3method of \code{smog}. If the \code{newdata} omitted,
#' it would output the prediction for the fitted model itself. The yielded result should
#' match with the family in the provided model. See \code{\link{smog}}.
#' 
#' @param object a fitted object of class inheriting from smog.
#' @param newdata a data frame containing the predictor variables, which are
#'                used to predict. If omitted, the fitted linear predictors 
#'                are used. 
#' @param family  a description of distribution family for which the response 
#'                variable is to be predicted.  
#' @param ... additional arguments affecting the predictions produced.
#' 
#' @details If \code{newdata = NULL}, the fitted.value based on the \code{object}
#'          is used for the prediction. For family = "coxph", the returned prediction value
#'          is the risk score. 
#' 
#' @return If \code{family} = "gaussian", a vector of prediction for the response is returned.
#'         For \code{family} = "coxph", a vector of predicted risk score is returned. 
#'         When \code{family} = ''binomial'', it outputs a data frame containing the predicted 
#'         group labels and the corresponding probabilies. 
#' 
#' @seealso \code{\link{smog.default}}, \code{\link{smog.formula}}, \code{\link{cv.smog}}, \code{\link{plot.smog}}.
#' @author Chong Ma, \email{chongma8903@@gmail.com}.
#' 
#' @references \insertRef{ma2019structural}{smog}
#' 
#' @importFrom stats coef fitted 
#' 
#' @export
predict.smog <- function(object, newdata = NULL, family = "gaussian",...){
  if(is.null(newdata)){
    y <- fitted(object)
  }else{
    if(!is.null(object$formula)){
      x = model.matrix(object$formula, newdata)
    }else{
      x = newdata
    }
    
    if(nrow(coef(object))){
      if(family == "gaussian"){
        wx = cbind(rep(1,nrow(x)),x)
        y <- as.vector(as.matrix(wx[,coef(object)$Id+1]) %*% as.matrix(coef(object)$Estimate))
      }
      
      if(family == "binomial"){
        probTAB = exp(as.matrix(x[,coef(object)$Id])%*%as.matrix(coef(object)$Estimate))/
          (1+exp(as.matrix(x[,coef(object)$Id])%*%as.matrix(coef(object)$Estimate)))
        probTAB = cbind(1-rowSums(probTAB),probTAB)
        predClass = apply(probTAB,1,which.max)
        predProb = apply(probTAB,1,max)
        
        y = data.frame(Class = object$levels[predClass],
                       Prob = predProb)
      }
      
      if(family == "coxph"){
        # risk score 
        y = round(exp(as.matrix(x[,coef(object)$Id])%*%as.matrix(coef(object)$Estimate)),2)
      }
      
    }else{
      stop("the model is not appropriate for making predictions")
    }
  }
  
  y
}



#' plot method for objects of `smog` class
#' 
#' \code{plot.smog} can produce a panel of plots for the primal errors, dual errors, 
#' and the penalized log-likelihood values, based on the provided fitted model 
#' (\code{x}) in the S3method of \code{smog}. 
#' 
#' @param x a fitted object of class inheriting from smog.
#' @param type,xlab default line types and x axis labels for the panel of plots.
#' @param caption a list of y axes labels for the panel of plots. 
#' @param ... additional arguments that could be supplied to \code{\link[graphics]{plot.default}}
#'            and \code{\link[graphics]{par}}.
#' 
#' @details For the panel of three plots, the \code{xlab} is ''iterations'' and the
#'          \code{type} is ''l'', by default. The \code{ylab} are ''primal error'',
#'          ''dual error'',''log-likelihood'', respectively. This panel of plots can
#'          reflect the convergence performance for the algorithm used in \code{\link{smog}}.
#' 
#' @seealso \code{\link[graphics]{par}}, \code{\link[graphics]{plot.default}},
#'          \code{\link{predict.smog}}, \code{\link{smog.default}}, \code{\link{cv.smog}}.
#'          
#' @author Chong Ma, \email{chongma8903@@gmail.com}.
#' 
#' @references \insertRef{ma2019structural}{smog}
#' 
#' @importFrom graphics plot par 
#' 
#' @export
plot.smog <- function(x,type = "l",xlab="iteration",caption=list("primal error","dual error","log-likelihood"),...){
  op<-graphics::par(mfrow=c(2,2),mar=c(2.1,4.1,2.1,2.1),no.readonly = TRUE,...)
  graphics::plot(x$PrimalError, type=type, xlab = xlab, ylab = caption[[1]],...)
  graphics::plot(x$DualError, type=type, xlab = xlab, ylab = caption[[2]],...)
  graphics::plot(x$loglike, type=type, xlab = xlab, ylab = caption[[3]],...)
  on.exit(par(op))
}











