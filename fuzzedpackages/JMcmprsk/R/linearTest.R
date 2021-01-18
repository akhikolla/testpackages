##' Joint modelling for longitutal and censored data with competing risks
##' @title Linear hypothesis testing of joint models
##' @param object  The JMcmprsk object returned by either jmo or jmc function.
##' @param coeff  Types of coefficients selected for Wald. Note "alpha" is only avaiable to jmo type JMcmprsk object.
##' @param La Linear contrast of the fixed effects of non-proportional odds covariates * (# of levels of the outcome - 2) in the longitudinal part. 
##' Default is "identity", i.e., all the fixed effects equal to zero. Otherwise, La must be a matrix.
##' @param Lb Linear contrast of the fixed effects of proportional odds covariates in the longitudinal part. 
##' Default is "identity", i.e., all the fixed effects equal to zero. Otherwise, Lb must be a matrix.
##' @param Lg Linear contrast of the fixed effects of covariates * # of competing risks in the survival part. 
##' Default is "identity", i.e., all the fixed effects equal to zero. Otherwise, Lg must be a matrix.
##' @param  Ca The hypothesized value of linear combination of the fixed effects of non-proportional odds covariates 
##' * (# of levels of the outcome - 2) in the longitudinal part. Default is 0. Otherwise, Ca must be a number / vector.
##' @param Cb The hypothesized value of linear combination of the fixed effects of proportional odds covariates 
##' in the longitudinal part. Default is 0. Otherwise, Cb must be a number / vector.
##' @param Cg The hypothesized value of linear combination of the fixed effects of covariates * # of competing risks in the survival part.
##' Default is 0. Otherwise, Cg must be a number / vector.
##' @param digits number of digits to be printed out.
##' @param ... further arguments passed to or from other methods.
##' 
##' @details Wald test statistic is used for hypothesis testing on multiple parameters: 
##' 
##'   \eqn{H_0: L\theta = C} vs: \eqn{H_1: L\theta \neq C}
##' 
##'   The test statistic is:
##'   
##'   \eqn{(L\hat{\theta} - C)'(L\hat{V_{\theta}}L)^{-1}(L\hat{\theta} - C) \sim \chi_q^2,}
##'    
##'   where \eqn{\hat{V_{\theta}}} is the estimate of covariance of the parameter \eqn{\theta} and q is the rank of the linear contrast \eqn{L}.
##' 
##' @return Return a Wald test statistic and the p value
##'   \tabular{ll}{
##'       \code{beta}    \tab  The Wald test for fixed effects for the longitutal part,i.e.  \eqn{\beta} in jmo or jmc output. \cr
##'       \code{gamma}    \tab  The Wald test for fixed effects for the survival part,i.e.  \eqn{\gamma} in jmo or jmc output.  \cr
##'       \code{alpha}    \tab  The Wald test for non-proportional odds covariates,i.e.  \eqn{\alpha} in jmo output. \cr
##'   }
##' @export
linearTest <-
  function (object,coeff=c("beta","gamma","alpha"), La = "identity", Lb = "identity", Lg = "identity",
            Ca = 0, Cb = 0, Cg = 0, digits = 4, ...) {
    if (!inherits(object, "JMcmprsk"))
      stop("Use only with 'JMcmprsk' objects.\n")
    if (object$type == "jmo") {
      if (coeff=="beta") {
        betas <- object$betas
        nbeta=length(betas)
        betacov<- object$vcmatrix[1:nbeta,1:nbeta]
        betas <- as.matrix(betas)
        if (!is.numeric(Lb)) {
          if (Lb == "identity") {
            W=t(betas - Cb)%*%(solve(betacov))%*%(betas - Cb)
            pval = pchisq(W,nbeta,lower.tail = FALSE)
            pval = formatC(pval, digits = digits, width=10,format = "f",flag="-")
            WaldTest <- data.frame(Chisq = W, df = nbeta, 
                                   "Pr(>|Chi|)" = pval, check.names = FALSE, row.names ="L*beta=Cb")
          } else {
            stop("Unexpected argument: do you mean identity?")
          }
        } else if (is.matrix(Lb) & (ncol(Lb) == nbeta)) {
          Trbeta <- Lb %*% betas - Cb
          Trbetacov <- Lb %*% betacov %*% t(Lb)
          W=t(Trbeta)%*%(solve(Trbetacov))%*%Trbeta
          df=pracma::Rank(Lb)
          pval = pchisq(W,df,lower.tail = FALSE)
          pval = formatC(pval, digits = digits, width=10,format = "f",flag="-")
          WaldTest <- data.frame(Chisq = W, df = df, 
                                 "Pr(>|Chi|)" = pval, check.names = FALSE, row.names ="L*beta=Cb")
        } else {
          stop("Linear contrast 'Lb' wasn't set up properly. 
                        Check if 'Lb' is matrix type object, dimension of Cb, or the number of columns in 'Lb' mataches.\n")
        }

        

      }else if (coeff=="alpha") {
        nbeta=length(object$betas)
        #number of nonpropotional hazards 
        alphas=as.vector(t(object$alphamatrix))
        nalpha=length(alphas)
        alphas <- as.matrix(alphas)
        alphacov<- object$vcmatrix[(nbeta+1):(nbeta+nalpha),
                                   (nbeta+1):(nbeta+nalpha)]
        if (!is.numeric(La)) {
          if (La == "identity") {
            W=t(alphas - Ca)%*%(solve(alphacov))%*%(alphas - Ca)
            pval = pchisq(W,nalpha,lower.tail = FALSE)
            pval = formatC(pval, digits = digits, width=10,format = "f",flag="-")
            WaldTest <- data.frame(Chisq = W, df = nalpha, 
                                   "Pr(>|Chi|)" = pval, check.names = FALSE, row.names ="L*alpha=Ca")
          } else {
            stop("Unexpected argument: do you mean identity?")
          }
        } else if (is.matrix(La) & (ncol(La) == nalpha)) {
          Tralpha = La %*% alphas - Ca
          Tralphacov = La %*% alphacov %*% t(La)
          W=t(Tralpha)%*%(solve(Tralphacov))%*%(Tralpha)
          df=pracma::Rank(La)
          pval = pchisq(W,df,lower.tail = FALSE)
          pval = formatC(pval, digits = digits, width=10,format = "f",flag="-")
          WaldTest <- data.frame(Chisq = W, df = df, 
                                 "Pr(>|Chi|)" = pval, check.names = FALSE, row.names ="L*alpha=Ca")
        } else {
          stop("Linear contrast 'La' wasn't set up properly. 
                        Check if 'La' is matrix type object, dimension of Ca, or the number of columns in 'La' mataches.\n")
        }
        
        
      } else if (coeff=="gamma") {
        nbeta=length(object$betas)
        alphas=as.vector(t(object$alphamatrix))
        thetas=object$Nlevels - 1
        nalpha=length(alphas)
        #number of competing risks
        gammas=as.vector(t(object$gamma_matrix))
        ngamma=length(gammas)
        gammacov<- object$vcmatrix[(nbeta+nalpha+thetas+1):(nbeta+nalpha+thetas+ngamma),
                                   (nbeta+nalpha+thetas+1):(nbeta+nalpha+thetas+ngamma)]
        if (!is.numeric(Lg)) {
          if (Lg == "identity") {
            W=t(gammas - Cg)%*%(solve(gammacov))%*%(gammas - Cg)
            pval = pchisq(W,ngamma,lower.tail = FALSE)
            pval = formatC(pval, digits = digits, width=10,format = "f",flag="-")
            WaldTest <- data.frame(Chisq = W, df = ngamma, 
                                   "Pr(>|Chi|)" = pval, check.names = FALSE, row.names ="L*gamma=Cg")
          } else {
            stop("Unexpected argument: do you mean identity?")
          }
        } else if (is.matrix(Lg) & (ncol(Lg) == ngamma)) {
          Trgamma = Lg %*% gammas - Cg
          Trgammacov = Lg %*% gammacov %*% t(Lg)
          W=t(Trgamma)%*%(solve(Trgammacov))%*%(Trgamma)
          df=pracma::Rank(Lg)
          pval = pchisq(W,df,lower.tail = FALSE)
          pval = formatC(pval, digits = digits, width=10,format = "f",flag="-")
          WaldTest <- data.frame(Chisq = W, df = df, 
                                 "Pr(>|Chi|)" = pval, check.names = FALSE, row.names ="L*gamma=Cg")
        } else {
          stop("Linear contrast 'Lg' wasn't set up properly. 
                        Check if 'Lg' is matrix type object, dimension of Cg, or the number of columns in 'Lg' mataches.\n")
        }
        
        
        
  } else {#change from print to message
    message("Error or Methods not developed yet!")
    WaldTest=NULL
   }      
      
      WaldTest
    } else if (object$type == "jmc")  {
      if (coeff=="beta") {
        #skip the first intercept
        nbeta=length(object$betas)
        betas <- object$betas[-1]
        betacov<- object$vcmatrix[2:nbeta,2:nbeta]
        nbeta=nbeta-1
        betas <- as.matrix(betas)
        if (!is.numeric(Lb)) {
          if (Lb == "identity") {
            W=t(betas - Cb)%*%(solve(betacov))%*%(betas - Cb)
            pval = pchisq(W,nbeta,lower.tail = FALSE)
            pval = formatC(pval, digits = digits, width=10,format = "f",flag="-")
            WaldTest <- data.frame(Chisq = W, df = nbeta, 
                                   "Pr(>|Chi|)" = pval, check.names = FALSE, row.names ="L*beta=Cb")
          } else {
            stop("Unexpected argument: do you mean identity?")
          }
          
        } else if (is.matrix(Lb) & (ncol(Lb) == nbeta)) {
          Trbeta <- Lb %*% betas - Cb
          Trbetacov <- Lb %*% betacov %*% t(Lb)
          W=t(Trbeta)%*%(solve(Trbetacov))%*%Trbeta
          df=pracma::Rank(Lb)
          pval = pchisq(W,df,lower.tail = FALSE)
          pval = formatC(pval, digits = digits, width=10,format = "f",flag="-")
          WaldTest <- data.frame(Chisq = W, df = df, 
                                 "Pr(>|Chi|)" = pval, check.names = FALSE, row.names ="L*beta=Cb")
        } else {
          stop("Linear contrast 'Lb' wasn't set up properly. 
                        Check if 'Lb' is matrix type object, dimension of Cb, or the number of columns in 'Lb' mataches.\n")
        }
        
        
      } else if (coeff=="gamma") {
        nbeta=length(object$betas)
        #number of competing risks
        gammas=as.vector(t(object$gamma_matrix))
        ngamma=length(gammas)
        gammacov<- object$vcmatrix[(nbeta+1+1):(nbeta+1+ngamma),
                                   (nbeta+1+1):(nbeta+1+ngamma)]
        if (!is.numeric(Lg)) {
          if (Lg == "identity") {
            W=t(gammas - Cg)%*%(solve(gammacov))%*%(gammas - Cg)
            pval = pchisq(W,ngamma,lower.tail = FALSE)
            pval = formatC(pval, digits = digits, width=10,format = "f",flag="-")
            WaldTest <- data.frame(Chisq = W, df = ngamma, 
                                   "Pr(>|Chi|)" = pval, check.names = FALSE, row.names ="L*gamma=Cg")
          } else {
            stop("Unexpected argument: do you mean identity?")
          }
        } else if (is.matrix(Lg) & (ncol(Lg) == ngamma)) {
          Trgamma = Lg %*% gammas - Cg
          Trgammacov = Lg %*% gammacov %*% t(Lg)
          W=t(Trgamma)%*%(solve(Trgammacov))%*%(Trgamma)
          df=pracma::Rank(Lg)
          pval = pchisq(W,df,lower.tail = FALSE)
          pval = formatC(pval, digits = digits, width=10,format = "f",flag="-")
          WaldTest <- data.frame(Chisq = W, df = df, 
                                 "Pr(>|Chi|)" = pval, check.names = FALSE, row.names ="L*gamma=Cg")
        } else {
          stop("Linear contrast 'Lg' wasn't set up properly. 
                        Check if 'Lg' is matrix type object, dimension of Cg, or the number of columns in 'Lg' mataches.\n")
        }
        
      }    else {#change from print to message
        message("Anova Error or Methods not developed yet!")
        WaldTest=NULL
      }     
      WaldTest
    }
}
