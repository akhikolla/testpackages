#' fit a regression for given lambda with network-based regularization
#'
#' Network-based penalization regression for given values of \eqn{\lambda_{1}} and \eqn{\lambda_{2}}.
#' Typical usage is to have the cv.regnet function compute the optimal lambdas, then provide them to the
#' regnet function. Users could also use MCP or Lasso.
#'
#' @keywords models
#' @param X matrix of predictors without intercept. Each row should be an observation vector. A column of 1 will be added to the X matrix
#' by the program as the intercept.
#' @param Y response variable. For response="binary", Y should be a numeric vector with zeros and ones. For response="survival", Y should be a
#' two-column matrix with columns named 'time' and 'status'. The latter is a binary variable, with '1' indicating a event, and '0'
#' indicating censoring.
#' @param response response type. regnet now supports three types of response: "binary", "continuous" and "survival".
#' @param penalty penalty type. regnet provides three choices for the penalty function: "network", "mcp" and "lasso".
#' @param lamb.1 the tuning parameter \eqn{\lambda_{1}} that imposes sparsity.
#' @param lamb.2 the tuning parameter \eqn{\lambda_{2}} that controls the smoothness among coefficient profiles. \eqn{\lambda_{2}} is  needed
#' for network penalty.
#' @param r the regularization parameter in MCP. For binary response, r should be larger than 4.
#' @param clv a value or a vector, indexing variables that are not subject to penalty. clv only works for continuous and survival responses
#' for now, and will be ignored for other types of responses.
#' @param initiation method for initiating the coefficient vector. For binary and continuous response, the default is elastic-net,
#' and for survival response the default is zero.
#' @param alpha.i the elastic-net mixing parameter. The program can use the elastic-net for choosing initial values of
#' the coefficient vector. alpha.i is the elastic-net mixing parameter, with 0 \eqn{\le} alpha.i \eqn{\le} 1. alpha.i=1 is the
#' lasso penalty, and alpha.i=0 is the ridge penalty. If the user chooses a method other than elastic-net for initializing
#' coefficients, alpha.i will be ignored.
#' @param robust logical flag. Whether or not to use robust methods. Robust methods are only available for survival response.
#'
#' @details The current version of regnet supports two types of responses: “binary”, "continuous" and “survival”.
#' regnet(…, response="binary", penalty="network") fits a network-based penalized logistic regression;
#' regnet(…, response="continuous", penalty="network") fits a network-based least square regression;
#' regnet(…, response="survival", penalty="network") fits a robust regularized AFT model using network penalty.
#' Please see the references for more details of the models. By default, regnet uses robust methods for survival response.
#' If users would like to use non-robust methods, simply set robust=FALSE. User could also use MCP or Lasso penalty.
#'
#' The coefficients are always estimated on a standardized X matrix. regnet standardizes each columns of X to have unit variance
#' (using 1/n rather than 1/(n-1) formula). If the coefficients on the original scale are needed, the user can refit a standard model
#' using the subset of variables that have non-zero coefficients.
#'
#' @return a vector of estimated coefficients. Please note that, if there are variables not subject to penalty (indicated by clv),
#' the order of returned vector is c(Intercept, unpenalized coefficients of clv variables, penalized coefficients of other variables).
#'
#' @references Ren, J., He, T., Li, Y., Liu, S., Du, Y., Jiang, Y., and Wu, C. (2017).
#' Network-based regularization for high dimensional SNP data in the case-control study of
#' Type 2 diabetes. \href{https://doi.org/10.1186/s12863-017-0495-5}{\emph{BMC Genetics}, 18(1):44}
#'
#' Ren, J., Du, Y., Li, S., Ma, S., Jiang,Y. and Wu, C. (2019). Robust network-based regularization
#' and variable selection for high dimensional genomics data in cancer prognosis.
#' \href{https://doi.org/10.1002/gepi.22194}{\emph{Genet. Epidemiol.}, 43:276-291}
#'
#' @seealso \code{\link{cv.regnet}}
#'
#' @examples
#' ## Survival response
#' data(SurvExample)
#' X = rgn.surv$X
#' Y = rgn.surv$Y
#' clv = c(1:5) # variables 1 to 5 are clinical variables which we choose not to penalize.
#' penalty = "network"
#' b = regnet(X, Y, "survival", penalty, rgn.surv$lamb1, rgn.surv$lamb2, clv=clv, robust=TRUE)
#' index = which(rgn.surv$beta != 0)
#' pos = which(b != 0)
#' tp = length(intersect(index, pos))
#' fp = length(pos) - tp
#' list(tp=tp, fp=fp)
#'
#' @export

regnet <- function(X, Y, response=c("binary", "continuous", "survival"), penalty=c("network", "mcp", "lasso"), lamb.1=NULL, lamb.2=NULL,
                   r=NULL, clv=NULL, initiation=NULL, alpha.i=1, robust=FALSE)
{
  # intercept = TRUE
  standardize=TRUE
  response = match.arg(response)
  penalty = match.arg(penalty)
  if(penalty != "network") lamb.2 = 0
  if(missing(lamb.1)) stop("Both lambda1 and lambda2 need to be provided")
  if(missing(lamb.2)) stop("Lambda2 needs to be provided for network method")

  this.call = match.call()
  if(response == "survival"){
    if(ncol(Y) != 2) stop("Y should be a two-column matrix")
    if(!setequal(colnames(Y), c("time", "status"))) stop("Y should be a two-column matrix with columns named 'time' and 'status'")
    Y0 = Y[,"time"]
    status = Y[,"status"]
    if(sum(Y0<=0)>0) stop("Survival times need to be positive")
  }else{
    if(robust) message("Robust methods are not available for ", response, " response.")
  }
  if(alpha.i>1 | alpha.i<0) stop("alpha.i should be between 0 and 1")
  if(is.null(initiation)){
    if(response == "survival") initiation = "zero" else initiation = "elnet"
  }
  if(is.null(r)) r = 5
  alpha = alpha.i # temporary

  fit=switch (response,
              "binary" = LogitCD(X, Y, penalty, lamb.1, lamb.2, r, alpha, init=initiation, alpha.i,standardize),
              "continuous" = ContCD(X, Y, penalty, lamb.1, lamb.2, clv, r, alpha, init=initiation, alpha.i,standardize),
              "survival" = SurvCD(X, Y0, status, penalty, lamb.1, lamb.2, clv, r, init=initiation, alpha.i, robust, standardize)
  )
  # fit$call = this.call
  # class(fit) = "regnet"
  fit
}
