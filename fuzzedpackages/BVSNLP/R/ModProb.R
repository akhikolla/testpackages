#' Logarithm of unnormalized probability of a given model
#'
#' @description This function calculates the logarithm of unnormalized
#' probability of a given set of covariates for both survival and binary
#' response data. It uses the inverse moment nonlocal prior (iMOM) for
#' non zero coefficients and beta binomial prior for the model space.
#' @param X The design matrix. It is assumed that the design matrix has
#' standardized columns. It is recommended that to use the output of
#' \code{\link{PreProcess}} function of the package.
#' @param resp For logistic regression models, this variable is the binary
#' response vector. For Cox proportional hazard models this is a two column
#' matrix where the first column contains the survival time vector and the
#' second column is the censoring status for each observation.
#' @param mod_cols A vector of column indices of the design matrix,
#' representing the model.
#' @param nlptype Determines the type of nonlocal prior that is used in the
#' analyses. It can be "piMOM" for product inverse moment prior, or "pMOM" for
#' product moment prior. The default is set to piMOM prior.
#' @param tau Hyperparameter \code{tau} of the iMOM prior.
#' @param r Hyperparameter \code{r} of the iMOM prior.
#' @param a First parameter in the beta binomial prior.
#' @param b Second parameter in the beta binomial prior.
#' @param family Determines the type of data analysis. \code{logistic} is for
#' binary outcome and logistic regression model whereas,
#' \code{survival} represents survival outcomes and the Cox proportional
#' hazard model.
#' @return It returns the unnormalized probability for the selected model.
#' @seealso \code{\link{CoefEst}}
#' @author Amir Nikooienejad
#' @references Nikooienejad, A., Wang, W., and Johnson, V. E. (2016). Bayesian
#' variable selection for binary outcomes in high dimensional genomic studies
#' using nonlocal priors. Bioinformatics, 32(9), 1338-1345.\cr\cr
#' Nikooienejad, A., Wang, W., & Johnson, V. E. (2020). Bayesian variable
#' selection for survival data using inverse moment priors. Annals of Applied
#' Statistics, 14(2), 809-828.
#' @examples
#' ### Simulating Logistic Regression Data
#' n <- 400
#' p <- 1000
#' set.seed(123)
#' Sigma <- diag(p)
#' full <- matrix(c(rep(0.5, p*p)), ncol=p)
#' Sigma <- full + 0.5*Sigma
#' cholS <- chol(Sigma)
#' Beta <- c(1,1.8,2.5)
#' X = matrix(rnorm(n*p), ncol=p)
#' X = X%*%cholS
#' beta <- numeric(p)
#' beta[c(1:length(Beta))] <- Beta
#' XB <- X%*%beta
#' probs <- as.vector(exp(XB)/(1+exp(XB)))
#' y <- rbinom(n,1,probs)
#'
#' ### Calling the function for a subset of the true model, with an arbitrary
#' ### parameters for prior densities
#' mod <- c(1:3)
#' Mprob <- ModProb(X, y, mod, tau = 0.7, r = 1, a = 7, b = 993,
#'                  family = "logistic")
#'
#' Mprob
ModProb <- function(X, resp, mod_cols, nlptype = "piMOM", tau, r, a, b,
                    family = c("logistic", "survival")){

  if(nlptype=="piMOM") nlptype_int <- 0
  if(nlptype=="pMOM") nlptype_int <- 1
  if (family == "logistic"){
    exmat <- cbind(resp, X)
    out <- lreg_mod_prob(exmat, mod_cols, tau, r, a, b, nlptype_int)
  } else {
    exmat <- cbind(resp, X)
    out <- cox_mod_prob(exmat, mod_cols, tau, r, a, b, nlptype_int)
  }
  return(out)
}
