#' Coefficient estimation for a specific set of covariates
#'
#' @description This function estimates coefficient vector for a given set of
#' covariates in a logistic regression and Cox proportional hazard models. It
#' uses the inverse moment nonlocal prior (iMOM) for non zero coefficients.
#' @param X The design matrix. It is assumed that the preprocessing steps have
#' been done on this matrix. It is recommended that to use the output of
#' \code{\link{PreProcess}} function of the package. Also note that the
#' \code{X} should NOT have a vector of $1$'s as the first column. If the
#' coefficients of a selected model by \code{\link{bvs}} is to be estimated, it
#' is highly recommended that the design matrix that is one of the outputs of
#' the \code{bvs} function and is reported as \code{des_mat} to be used here.
#' @param resp For logistic regression models, this variable is the binary
#' response vector. For Cox proportional hazard models this is a two column
#' matrix where the first column contains the survival time vector and the
#' second column is the censoring status for each observation.
#' @param mod_cols A vector of column indices of the design matrix,
#' representing the selected model.
#' @param nlptype Determines the type of nonlocal prior that is used in the
#' analyses. It can be "piMOM" for product inverse moment prior, or "pMOM" for
#' product moment prior. The default is set to piMOM prior.
#' @param tau Hyperparameter \code{tau} of the iMOM prior.
#' @param r Hyperparameter \code{r} of the iMOM prior.
#' @param family Determines the type of data analysis. \code{logistic} is for
#' binary outcome and logistic regression model whereas,
#' \code{survival} represents survival outcomes and the Cox proportional
#' hazard model.
#' @return It returns the vector of coefficients for the given model.
#' @author Amir Nikooienejad
#' @references Nikooienejad, A., Wang, W., and Johnson, V. E. (2016). Bayesian
#' variable selection for binary outcomes in high dimensional genomic studies
#' using non-local priors. Bioinformatics, 32(9), 1338-1345.\cr\cr
#' Nikooienejad, A., Wang, W., & Johnson, V. E. (2020). Bayesian variable
#' selection for survival data using inverse moment priors. Annals of Applied
#' Statistics, 14(2), 809-828.
#' @seealso \code{\link{ModProb}}
#' @examples
#' ### Simulating Survival Data
#' n <- 400
#' p <- 1000
#' lambda <- 0.8
#' cens_rate <- 0.27
#' set.seed(123)
#' Sigma <- diag(p)
#' full <- matrix(c(rep(0.5, p*p)), ncol=p)
#' Sigma <- full + 0.5*Sigma
#' cholS <- chol(Sigma)
#' Beta <- c(-1.8, 1.2, -1.7, 1.4, -1.4, 1.3)
#' X = matrix(rnorm(n*p), ncol=p)
#' X = X%*%cholS
#' X <- scale(X)
#' beta <- numeric(p)
#' beta[c(1:length(Beta))] <- Beta

#' XB <- X%*%beta
#' uvector <- -log(runif(n));
#' times <- uvector/(lambda*exp(XB))
#' cens_time <- quantile(times,1-cens_rate)
#' status <- as.numeric(times < cens_time)
#' TS <- cbind(times,status)
#'
#' ### Estimating coeffcients of the true model and an arbitrary hyper
#' ### parameter for the iMOM prior density
#' mod <- c(1:6)
#' coef <- CoefEst(X, TS, mod, tau = 1.8, r = 2, family = "survival")
#' coef
#'
CoefEst <- function(X, resp, mod_cols, nlptype = "piMOM", tau, r,
                        family = c("logistic", "survival")){
  
  if(nlptype=="piMOM") nlptype_int <- 0
  if(nlptype=="pMOM") nlptype_int <- 1
  if (family == "logistic"){
    X <- cbind(rep(1,length(resp)),X)
    exmat <- cbind(resp, X)
    out <- lreg_coef_est(exmat, mod_cols, tau, r, nlptype_int)
  } else {
    exmat <- cbind(resp, X)
    out <- cox_coef_est(exmat, mod_cols, tau, r, nlptype_int)
  }
  return(out)
}


