#'@title Formula Variables
#'@aliases formula
#'@name formula
#'@description Transforming a formula object into a list with the variables and their names for the beta regression model of the bayesbr package.
#'@usage formula(formula, data = NULL)
#'@param formula symbolic description of the model (of type \code{y ~ x} or \code{y ~ x | z};).
#'@param data Data frame with regression observations
#'@details The form of the formula used for the \code{Bayesbr} package follows the pattern proposed in \code{\link{Formula}}. The expression y ~ represents that y is the response variable of the beta regression, everything to the right of the ~ operator represents covariates or intercepts for the parameter \eqn{\theta} or \eqn{\zeta} of the variable response .
#'
#'The + operator adds one more explanatory covariate for the parameter,the operator : indicates interaction between variables adjacent to the operator,
#' operator * adds the variables adjacent to the operator as covariable and the interaction between them
#'the operator | represents that the next covariates are explanatory
#'for \eqn{\zeta} and those that were before the operator are explanatory
#'for \eqn{\theta}. So, in the formula \code{y ~ x1 + x2 | x3 + x4} x1 and x2
#'are the covariates for the parameter \eqn{\theta} and x3 and x4 are the
#' covariates of \eqn{\zeta}. \eqn{\theta} and \eqn{\zeta} are parameters
#'  of the variable y answer. The numbers 1 and 0 represent, respectively,
#'   the presence or not of the intercept in the construction of the model.
#'   By default, the intercept is included, so the number 1 is
#'   necessary only when the user wants to include only the intercept for
#'   the estimation of the parameter in question. Here are some examples:
#'
#'\code{y ~ 0 | x1}: No estimate for \eqn{\theta}
#'
#'\code{y ~ 1 | 0 + x2}: The estimation for \eqn{\theta} will be made only with the intercept, and the estimation for \eqn{\zeta} will not use the intercept only the covariable x2
#'
#'\code{y~ x3*x4 | x5:x6}: The estimation for \eqn{\theta} will be with the covariables \code{x3} and \code{x4} and the interaction between them, and the estimation for \eqn{\zeta} will be the interaction between variables \code{x5} and \code{x6}.
#'
#'The variables passed to the formula can be environment variables or columns of a dataframe, in which case the dataframe must be informed.
#'
#'@return A list containing the following items:
#'\describe{\item{Y}{A vector containing the model response variable,}
#'\item{X}{A matrix containing the covariates for theta of the model,}
#'\item{W}{A matrix containing the covariates of the model for zeta}
#'\item{name_y}{The name passed in the call to the \code{\link{bayesbr}} function for the variable response,}
#'\item{name_x}{The name passed in the call to the \code{\link{bayesbr}} function for the covariates for theta,}
#'\item{name_w}{The name passed in the call to the \code{\link{bayesbr}} function for covariates for zeta.}}
#'@seealso \code{\link{bayesbr}}
formula = function(formula,data = NULL){
  if(is.null(data)){
    data = environment()
  }
  formula <- Formula(formula)
  if(length(formula)[2]>2){
    stop('invalid template formula')
  }
  model_frame <- model.frame(formula, data = data)
  Y = model.response(model_frame)
  if(is.null(Y)){
    stop('invalid template formula')
  }
  name_y = all.vars(formula)[1]
  X = suppressWarnings(as.data.frame(model.matrix(formula, data = model_frame, rhs = 1)))
  names_x = suppressWarnings(colnames(X))
  W = suppressWarnings(as.data.frame(model.matrix(formula, data = model_frame, rhs = 2)))
  names_w = suppressWarnings(colnames(W))
  if(is.null(names_x)){
    X = NULL
  }
  else{X = as.matrix(X)}
  if(is.null(names_w)){
    W = NULL
  }
  else{W = as.matrix(W)}
  list = list(y=Y,x = X,w=W, name_y = name_y, names_x = names_x, names_w = names_w)
  return(list)
}


