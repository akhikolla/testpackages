#' Solve Henderson's Mixed Model Equation.
#'
#' Consider a linear mixed model with normal random effects,
#' \deqn{Y_{ij} = X_{ij}^T\beta + v_i + \epsilon_{ij}}
#' where \eqn{i=1,\ldots,n,\quad j=1,\ldots,m}, or it can be equivalently expressed using matrix notation,
#' \deqn{Y = X\beta + Zv + \epsilon}
#' where \eqn{Y\in \mathrm{R}^{nm}} is a known vector of observations, \eqn{X \in \mathrm{R}^{nm\times p}} and
#' \eqn{Z \in \mathrm{R}^{nm\times n} } design matrices for \eqn{\beta} and \eqn{v} respectively,
#' \eqn{\beta \in \mathrm{R}^p} and \eqn{v\in \mathrm{R}^n} unknown vectors of fixed effects and random effects where
#' \eqn{v_i \sim N(0,\lambda_i)}, and
#' \eqn{\epsilon \in \mathrm{R}^{nm}} an unknown vector random errors independent of random effects. Note that \eqn{Z}
#' does not need to be provided by a user since it is automatically created accordingly to the problem specification.
#'
#' @param X an \eqn{(nm\times p)} design matrix for \eqn{\beta}.
#' @param Y a length-\eqn{nm} vector of observations.
#' @param Mu a length-\eqn{nm} vector of initial values for \eqn{\mu_i = E(Y_i)}.
#' @param Lambda a length-\eqn{n} vector of initial values for \eqn{\lambda}, variance of \eqn{v_i \sim N(0,\lambda_i)}
#'
#' @return a named list containing \describe{
#' \item{beta}{a length-\eqn{p} vector of BLUE \eqn{\hat{beta}}.}
#' \item{v}{a length-\eqn{n} vector of BLUP \eqn{\hat{v}}.}
#' \item{leverage}{a length-\eqn{(mn+n)} vector of leverages.}
#' }
#'
#'
#' @examples
#' ## small setting for data generation
#' n = 100; m = 2; p = 2
#' nm = n*m;   nmp = n*m*p
#'
#' ## generate artifical data
#' X = matrix(rnorm(nmp, 2,1), nm,p) # design matrix
#' Y = rnorm(nm, 2,1)                # observation
#'
#' Mu = rep(1, times=nm)
#' Lambda = rep(1, times=n)
#'
#' ## solve
#' ans = SolveHMME(X, Y, Mu, Lambda)
#'
#'
#' @references
#' \insertRef{henderson_estimation_1959}{HMMEsolver}
#'
#' \insertRef{robinson_that_1991}{HMMEsolver}
#'
#' \insertRef{mclean_unified_1991}{HMMEsolver}
#'
#' \insertRef{kim_fast_2017}{HMMEsolver}
#'
#' @export
SolveHMME = function(X, Y, Mu, Lambda){
  #--------------------------------------------------------------
  ## Part 1 : Parameter Checks
  #   1. get size information of the problem
  Rn = length(Lambda)
  dimMat = dim(X)
  Rp = dimMat[2]
  Rm = dimMat[1] / Rn
  #   2. X
  if (!check_matrix(X,n=Rn*Rm,p=p)){
    stop("SolveHMME : an input design matrix X is invalid.")
  }
  #   3. Y
  Y = as.vector(Y)
  if (!check_vector(Y,n=Rn*Rm)){
    stop("SolveHMME : an input observation vector Y is invalid.")
  }
  #   4. Mu
  Mu = as.vector(Mu)
  if (!check_vector(Mu,n=Rn*Rm)){
    stop("SolveHMME : an input vector of initial values is invalid.")
  }
  #   5. Lambda
  Lambda = as.vector(Lambda)
  if ((!check_vector(Lambda,n=Rn))||(any(Lambda<=0))){
    stop("SolveHMME : an input vector of variance parameters is invalid.")
  }

  #--------------------------------------------------------------
  ## Part 2 : Main Computation by Jiwoong
  # RnModel = nModel1
  # RnModel2 = nModel12
  #
  # nModel = integer(1)
  # nModel[1]=as.integer(RnModel)
  #
  # nModel2 = integer(1)
  # nModel2[1]=as.integer(RnModel2)


  n = integer(1)
  n[1] = as.integer(Rn)

  m = integer(1)
  m[1] = as.integer(Rm)

  p = integer(1)
  p[1] = as.integer(Rp)

  nm = as.integer(Rn*Rm)
  nmp = as.integer(Rn*Rm*Rp)
  x = double(nmp)

  for(i in 1:Rn){
    for(j in 1:Rm){

      for(k in 1:Rp){
        x[(i-1)*m*p+(j-1)*p+k] = X[(i-1)*m+j, k]
      }

    }
  }

  y = double(nm)
  mu = double(nm)

  for(i in 1:nm){
    y[i] = Y[i]
    mu[i]= Mu[i]
  }

  lambda = double(Rn)

  for(i in 1:Rn){
    lambda[i] = Lambda[i]
  }


  out = double(Rn*Rm+2*Rn+Rp)

  result = .C('SolveHMM', output=out, x, y, mu, lambda, n, m, p)
  output = result$output
  #--------------------------------------------------------------
  ## Part 3 : Wrapping the results
  # first p elements - BLUE beta
  # next  n elements - BLUP v
  # rest             - leverages' diagonal

  OutList = list()
  OutList$beta = output[1:Rp]
  OutList$v    = output[(Rp+1):(Rp+Rn)]
  OutList$leverage = output[(Rp+Rn+1):length(output)]
  return(OutList)
}
