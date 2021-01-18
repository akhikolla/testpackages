#' fit a Semi-parametric Bayesian variable selection
#'
#' fit a Bayesian semi-parametric model for both linear and non-linear G×E interactions. Users can also specify all the interactions as linear and fit a Bayesian LASSO type of model.
#'
#' @keywords models
#' @param X the matrix of predictors (genetic factors) without intercept. Each row should be an observation vector. A column of 1 will be added to the X matrix
#' as the intercept.
#' @param Y the response variable. The current version of BVCfit only supports continuous response.
#' @param Z a vector of environmental factor for non-linear G×E interactions.
#' @param E a vector of environmental factor for linear G×E interactions.
#' @param clin a matrix of clinical variables. Clinical variables are not subject to penalty.
#' @param iterations the number of MCMC iterations.
#' @param burn.in the number of iterations for burn-in.
#' @param sparse logical flag. If TRUE, spike-and-slab priors will be used to shrink coefficients of irrelevant covariates to zero exactly. 'sparse' has effect only when VC=TRUE.
#' @param structural logical flag. If TRUE, the coefficient functions with varying effects and constant effects will be penalized separately. 'structural' has effect only when VC=TRUE.
#' @param VC logical flag. If TRUE, varying coefficient functions will be used for modeling the interactions between Z and X.
#' If FALSE, interactions between Z and X will be modeled as linear interactions.
#' @param kn the number of interior knots for B-spline.
#' @param degree the degree of B spline basis.
#' @param hyper a named list of hyperparameters.
#' @param debugging logical flag. If TRUE, progress will be output to the console and extra information will be returned.
#'
#' @details By default, varying coefficient functions are used for modeling the nonlinear interactions between Z and X. Assuming both E and clin are NULL, the model can be expressed as
#' \deqn{Y = \beta_{0}(Z)+\sum\beta_{j}(Z)X_{j} + \epsilon }
#' The basis expansion and changing of basis with B splines will be done automatically:
#' \deqn{\beta_{j}(\cdot)\approx \gamma_{j1} + \sum_{k=2}^{q}{B}_{jk}(\cdot)\gamma_{jk}}
#' where \eqn{B_{jk}(\cdot)} represents B spline basis. \eqn{\gamma_{j1}} and \eqn{(\gamma_{j2}, \ldots, \gamma_{jq})^\top} correspond to the constant and varying parts of the coefficient functional, respectively.
#' q=kn+degree+1 is the number of basis functions. By default, kn=degree=2. User can change the values of kn and degree to any other positive integers.
#' If E is provided, the linear interactions between E and X will be added modeled as pairwise-products:
#' \deqn{Y = \beta_{0}(Z)+\sum\beta_{j}(Z)X_{j} + \zeta_{0}E + \sum \zeta_{j}EX_{j} + \epsilon}
#' If clin is provided, clinical variables
#' will be added to the model.
#'
#' If VC=FALSE, all interactions are treated as linear and a Bayesian LASSO model will be used. With non-null values of E and clin, the full linear model is:
#' \deqn{Y \sim Z + ZX + clin + E + EX}
#' Please check the references for more details about the model.
#'
#' Users can modify the hyper-parameters by providing a named list of hyper-parameters via the argument 'hyper'.
#' The list can have the following named components
#' \itemize{
#'   \item{a.c, a.v, a.e: }{ shape parameters of the Gamma priors on \eqn{\lambda_{c}}, \eqn{\lambda_{v}} and \eqn{\lambda_{e}}, respectively.}
#'   \item{b.c, b.v, b.e: }{ rate parameters of the Gamma priors on \eqn{\lambda_{c}}, \eqn{\lambda_{v}} and \eqn{\lambda_{e}}, respectively.}
#'   \item{r.c, r.v, r.e: }{ shape parameters of the Beta priors (\eqn{\pi^{r-1}(1-\pi)^{w-1}}) on \eqn{\pi_{c}}, \eqn{\pi_{v}} and \eqn{\pi_{e}}, respectively.}
#'   \item{w.c, w.v, w.e: }{ shape parameters of the Beta priors on \eqn{\pi_{c}}, \eqn{\pi_{v}} and \eqn{\pi_{e}}, respectively.}
#'   \item{s: }{ shape parameters of the Inverse-gamma prior on \eqn{\sigma^{2}}.}
#'   \item{h: }{ scale parameters of the Inverse-gamma prior on \eqn{\sigma^{2}}.}
#' }
#' Please check the references for more details about the prior distributions.
#'
#' @return an object of class "BVCfit" is returned, which is a list with components:
#' \item{posterior}{posterior samples from the MCMC}
#' \item{coefficients}{a list of posterior estimates of coefficients}
#' \item{burn.in}{the number of iterations for burn-in}
#' \item{iterations}{the number of MCMC iterations.}
#'
#' @references
#' Ren, J., Zhou, F., Li, X., Chen, Q., Zhang, H., Ma, S., Jiang, Y., Wu, C. (2019) Semi-parametric Bayesian variable selection for gene-environment interactions.
#' \url{https://arxiv.org/abs/1906.01057}
#'
#' @examples
#' data(gExp)
#'
#' ## default method
#' spbayes=BVCfit(X, Y, Z, E, clin)
#' spbayes
#'
#' \donttest{
#' ## non-structural
#' structural=FALSE
#' spbayes=BVCfit(X, Y, Z, E, clin, structural=structural)
#' spbayes
#'
#' ## non-sparse
#' sparse=FALSE
#' spbayes=BVCfit(X, Y, Z, E, clin, sparse=sparse)
#' spbayes
#' }
#'
#' @export

BVCfit <- function(X, Y, Z, E=NULL, clin=NULL, iterations=10000, burn.in=NULL, sparse=TRUE, structural=TRUE, VC=TRUE, kn=2, degree=2, hyper=NULL, debugging=FALSE)
{
  x = as.matrix(X); y = cbind(Y); kn = as.integer(kn); degree = as.integer(degree)
  q = kn + degree + 1
  n = nrow(x); s = ncol(x)
  noClin = noE = TRUE
  EX = ZX = NULL
  Z = cbind(Z)
  CLC = Z

  if(iterations<1) stop("iterations must be a positive integer.")
  if(is.null(burn.in)){
	BI = floor(iterations)/2
	if(iterations<=BI) stop("iterations must be larger than burn.in.")
  }else if(burn.in>=1){
	BI = as.integer(burn.in)
  }else{
	stop("burn.in must be a positive integer.")
  }

  if(any(c(length(Z), length(y)) != n)){
    stop("X, Y and Z have different numbers of observations.")
  }

  if(any(kn, degree) <= 0){
	stop("incorrect smoothness specification, kn and degree must be positive integers.")
  }

  if(!is.null(clin)){
    clin = as.matrix(clin)
    if(is.null(colnames(clin))){colnames(clin)=paste("clin.", 1:ncol(clin), sep="")}
    CLC = cbind(as.matrix(clin), Z=Z)
    noClin = FALSE
    Clin.names = colnames(clin)
  }

  if(!is.null(E)){
    CLC = cbind(E=E, CLC)
    EX = as.matrix(as.numeric(E) * x)
    noE = FALSE
  }

  if(!VC){
    ZX = as.matrix(as.numeric(Z) * x);
  }

  CLC.names = colnames(CLC)
  if(is.null(colnames(x))){
    Var.names = paste("G", 1:s, sep="")
  }else{
    Var.names = colnames(x)
  }

  x = cbind(1, x) # add intercept

  if(VC){
	design = Design.matrix(Z, x, kn, degree)
	xx = design$X
  }else{
	xx = x
  }
  # design = Design.matrix(Z, x, kn, degree)
  # xx = design$X
  xxwe = cbind(xx, CLC, EX, ZX)
  lasso.cv = glmnet::cv.glmnet(xxwe,y,alpha=1,nfolds=5)
  lambda.cv = lasso.cv$lambda.min;
  lasso.fit = glmnet::glmnet(xxwe, y, family="gaussian",alpha=1,nlambda=50)
  coeff.array = as.vector(stats::predict(lasso.fit, s=lambda.cv, type="coefficients"))[-2];

  nclc = ncol(CLC)
  if(debugging) message("No. of E: ", (!noE)*1, "\nNo. of clinical covariates: ", ifelse(noClin, 0, ncol(CLC)-1-!noE), "\n")

  if(VC){
	  hat.m = coeff.array[1:q]      ## coeff for intercept
	  hat.r0 = coeff.array[(1:s)+q] ## coeff for constant
	  hat.r.star = utils::head(coeff.array, ncol(xx))[-(1:(s+q))] ## coeff for varying part
  }

  coeff.clc = utils::tail(coeff.array, -ncol(xx)) ## E CLC Z EX ZX
  hat.clc = coeff.clc[1:nclc]              ## E CLC Z
  hat.zeta = utils::tail(coeff.clc, -nclc)  ## EX ZX

  if(!VC){
	  out = BLasso(xx, y, CLC, EX, ZX, s, iterations, coeff.array[1], coeff.array[2:ncol(xx)], hat.clc, hat.zeta, hyper, debugging)
	  CC = apply(out$posterior$GS.r0[-c(1:BI),,drop=FALSE], 2, stats::median)
    LL = apply(out$posterior$GS.rs[-c(1:BI),,drop=FALSE], 2, stats::median)
	  names(CC) = names(LL) = Var.names
	  coeff = list(intercept=stats::median(out$posterior$GS.m[-c(1:BI)]), Z=stats::median(out$posterior$GS.Z[-c(1:BI)]), Main=CC, Interaction=LL)
  }else{
    if(structural){
      out = BVC_SI(xx, y, CLC, EX, s, q, iterations, hat.m, hat.r0, hat.r.star, hat.clc, hat.zeta, sparse, hyper, debugging)
      INT = apply(out$posterior$GS.m[-c(1:BI),,drop=FALSE], 2, stats::median)
      CC = apply(out$posterior$GS.r0[-c(1:BI),,drop=FALSE], 2, stats::median)
      VV = apply(out$posterior$GS.rs[-c(1:BI),,drop=FALSE], 2, stats::median)
      coeff = cbind(INT, rbind(CC, matrix(VV, nrow = q-1)))
    }else{
      hat.r = c(rbind(hat.r0, matrix(hat.r.star, nrow = (q-1))))
      out = BVC_NS(design$Xns, y, CLC, EX, s, q, iterations, hat.m, hat.r, hat.clc, hat.zeta, sparse, hyper, debugging)
      INT = apply(out$posterior$GS.m[-c(1:BI),,drop=FALSE], 2, stats::median)
      VV = apply(out$posterior$GS.rs[-c(1:BI),,drop=FALSE], 2, stats::median)
      coeff = cbind(INT, matrix(VV, nrow = q))
    }
    colnames(coeff) = c("intercept",Var.names)
    rownames(coeff) = paste("basis", 0:(q-1), sep="")
  }

  this.call = match.call()
  basis = list(q=q, kn=kn, degree=degree)

  if(noE && noClin){
    coeff.clin = NULL
    coeff.E = NULL
    coeff.zeta = NULL
  }else if(!noE && !noClin){
    coeff.clc = apply(out$posterior$GS.clc[-c(1:BI),,drop=FALSE], 2, stats::median)
    coeff.clin = coeff.clc[-1]
    coeff.E = coeff.clc[1]
    coeff.zeta = apply(out$posterior$GS.zeta[-c(1:BI),,drop=FALSE], 2, stats::median)
    names(coeff.clin) = Clin.names
    names(coeff.zeta) = paste("G", 1:s, sep="")
  }else if(noE){
    coeff.clin = apply(out$posterior$GS.clc[-c(1:BI),,drop=FALSE], 2, stats::median)
    coeff.E = NULL
    coeff.zeta = NULL
    names(coeff.clin) = Clin.names
  }else{
    coeff.clin = NULL
    coeff.E = apply(out$posterior$GS.clc[-c(1:BI),,drop=FALSE], 2, stats::median)
    coeff.zeta = apply(out$posterior$GS.zeta[-c(1:BI),,drop=FALSE], 2, stats::median)
    names(coeff.zeta) = paste("G", 1:s, sep="")
  }

  coefficient = list(E=coeff.E, clin=coeff.clin, EX=coeff.zeta, ZX=coeff)
  if(noE){
    coefficient$E = NULL
    coefficient$EX = NULL
  }
  if(noClin) coefficient$clin = NULL

  fit = list(call = this.call, posterior = out$posterior, coefficient=coefficient, burn.in = BI, iterations=iterations)

  if(debugging && sparse) fit$debugList = out$debugList;
  if(VC) fit$basis = basis;

  class(fit)=c("BVCfit", class(out))
  fit
}
