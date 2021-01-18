#' Generalized Multi-SNP Mediation Intersection-Union Test
#' 
#' Testing the mediation effect of multiple SNPs on an outcome following an exponential family distribution or a survival outcome through a continuous mediator.
#' 
#' @param G n by p matrix (n rows and p columns). Each row is one individual; each column is one SNP.
#' @param mediator a vector length of n. It is the mediator variable.
#' @param outcome a vector length of n. It is the outcome variable.
#' @param covariates n by r matrix (n rows and r columns). Each row is one individual; each column is one covariate. Default is NULL.
#' @param outcome_type Type of the outcome variable. "continuous" for a continuous outcome; "binary" for a binary outcome; "count" for a count outcome; "survival" for a survival outcome. 
#' @param approxi a boolean value. This is an indicator whether the approximation of computing derivatives is applied to save computation time. Default is TRUE.
#' @param verbose a boolean value. If TRUE a lot of computing details is printed. Default is FALSE.
#' 
#' @return A list of the following components:
#' \itemize {
#'   \item \code{p_value_IUT}: The p value for testing the mediation effect (theta*beta) based on intersection-union test. 
#'   \item \code{p_value_theta}: The p value for testing theta in the outcome model. The outcome model is the following. outcome ~ intercept + covariates iota + G gamma + mediator*theta 
#'   \item \code{theta_hat}: The point estimate of theta (coefficient of mediator) in the outcome model. 
#'   \item \code{p_value_beta}: The p value for testing beta in the mediator model. The mediator model is the following. mediator ~ intercept + covariates iota + G beta + error
#' }
#' 
#' @author Wujuan Zhong
#' 
#' 
#' 
GSMUT=function(G,mediator,outcome,covariates=NULL,outcome_type,approxi=TRUE,verbose=FALSE){

  # p_value_theta is the p value for testing theta in the outcome model
  # outcome model: outcome ~ intercept + covariates*iota + G*gamma + mediator*theta 
  res=Generalized_Testing_coefficient_of_mediator(G,mediator,outcome,covariates,outcome_type,approxi,verbose)
  
  p_value_theta=res$pvalue
  theta_hat=res$theta_hat
  
  # p_value_beta is the p value for testing beta in the mediator model
  # mediator model: mediator ~ intercept + covariates*iota + G*beta + error
  if (is.null(covariates)){
    obj=SKAT_Null_Model(mediator~1,out_type="C")
    p_value_beta=SKAT( as.matrix(G), obj,impute.method="bestguess")$p.value
  }else{
    obj=SKAT_Null_Model(mediator~covariates,out_type="C")
    p_value_beta=SKAT( as.matrix(G), obj,impute.method="bestguess")$p.value
  }
  

  # p_value_IUT is the p value for testing theta*beta based on intersection-union test
  p_value_IUT=max(c(p_value_theta,p_value_beta))
  return(list("p_value_IUT"=p_value_IUT,"p_value_theta"=p_value_theta,"theta_hat"=theta_hat,"p_value_beta"=p_value_beta))
}


