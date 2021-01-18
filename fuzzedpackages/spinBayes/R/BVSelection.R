#' Variable selection for a BVCfit object
#'
#' Variable selection for a BVCfit object
#'
#' @param obj BVCfit object.
#' @param ... other BVSelection arguments
#'
#' @details For class 'BVCSparse', the median probability model (MPM) (Barbieri and Berger 2004) is used to identify predictors that are significantly associated
#' with the response variable. For class 'BVCNonSparse', variable selection is based on 95\% credible interval.
#' Please check the references for more details about the variable selection.
#'
#' @references
#' Ren, J., Zhou, F., Li, X., Chen, Q., Zhang, H., Ma, S., Jiang, Y., Wu, C. (2019) Semi-parametric Bayesian variable selection for gene-environment interactions.
#' \url{https://arxiv.org/abs/1906.01057}
#'
#' Barbieri, M.M. and Berger, J.O. (2004). Optimal predictive model selection
#' \href{https://projecteuclid.org/euclid.aos/1085408489}{\emph{Ann. Statist}, 32(3):870–897}
#'
#' @rdname BVSelection
#' @return an object of class "BVSelection" is returned, which is a list with components:
#' \item{method}{posterior samples from the MCMC}
#' \item{indices}{a list of indices and names of selected variables}
#' \item{summary}{a summary of selected variables}
#'
#' @seealso \code{\link{BVCfit}}
#'
#' @examples
#' data(gExp)
#' ## sparse
#' spbayes=BVCfit(X, Y, Z, E, clin)
#' spbayes
#'
#' selected = BVSelection(spbayes)
#' selected$indices
#'
#' ## non-sparse
#' spbayes=BVCfit(X, Y, Z, E, clin, sparse=FALSE)
#' spbayes
#'
#' selected = BVSelection(spbayes)
#' selected
#'
#' @export
BVSelection <- function(obj,...){
  UseMethod('BVSelection', obj)
}


#' @param burn.in MCMC burn-in.
#' @param prob probability for credible interval, between 0 and 1. e.g. prob=0.95 leads to 95\% credible interval
#' @rdname BVSelection
#' @method BVSelection BVCNonSparse
#' @export
BVSelection.BVCNonSparse=function(obj, burn.in=obj$burn.in, prob=0.95,...){
  BI = ifelse(is.null(burn.in), 0, burn.in)
  GS.r0 = obj$posterior$GS.r0
  GS.rs = obj$posterior$GS.rs
  GS.zeta = obj$posterior$GS.zeta
  q = obj$basis$q
  if(BI>0){
    GS.r0 = GS.r0[-c(1:BI),]
    GS.zeta = GS.zeta[-(1:BI),]
    GS.rs = GS.rs[-c(1:BI),]
  }

  if(is.null(GS.r0)){
    SelectR0 = 0
    SelectRstar = Selection.CI(GS.rs, q, prob)
  }else{
    SelectR0 = Selection.CI(GS.r0, 1, prob)
    SelectRstar = Selection.CI(GS.rs, max(1, q-1), prob)
  }
  SelectZeta = if(is.null(GS.zeta)){ 0 }else{ Selection.CI(GS.zeta, 1, prob) }

  MPM.V = which(SelectRstar > 0)
  MPM.C = setdiff(which(SelectR0 > 0), MPM.V)
  if(is.null(q)){
    MPM.C = which(SelectR0 > 0)
  }else{
    MPM.C = setdiff(which(SelectR0 > 0), MPM.V)
  }
  MPM.Z = which(SelectZeta > 0)
  numb = matrix(c(length(MPM.C), length(MPM.V), length(MPM.Z)), ncol=1,
                dimnames=list(c("Constant effect", "Varying effect", "Linear interaction"), "#"))

  Var.names = colnames(obj$coefficient$VC)[-1]
  if(length(MPM.C)>0){
    Main = MPM.C
    names(Main) = Var.names[MPM.C]
  }else{
    Main = NULL
  }

  if(length(MPM.V)>0){
    Varying = MPM.V
    names(Varying) = Var.names[MPM.V]
  }else{
    Varying = NULL
  }

  if(length(MPM.Z)>0){
    Linear = MPM.Z
    names(Linear) = names(obj$coefficient$EG)[MPM.Z]
  }else{
    Linear = NULL
  }

  if(is.null(q)){
    sel = list(Main=Main, Linear.ZX=Varying, Linear.EX=Linear)
    rownames(numb) = c("Main effect", "Linear interaction (ZX)", "Linear interaction (EX)")
  }else{
    sel = list(Constant=Main, Varying=Varying, Linear=Linear)
  }

  method = paste(prob*100,"% credible interval", sep = "")
  out = list(method=method, indices=sel, summary=numb)
  class(out) = "BVSelection"
  out
}



#' @rdname BVSelection
#' @method BVSelection BVCSparse
#' @export
BVSelection.BVCSparse=function(obj, burn.in=obj$burn.in,...){
  BI = ifelse(is.null(burn.in), 0, burn.in)
  max_BI = obj$iterations - BI
  GS.r0 = obj$posterior$GS.r0
  GS.zeta = obj$posterior$GS.zeta
  GS.phi = obj$posterior$GS.phi
  if(BI>0){
    GS.r0 = GS.r0[-c(1:BI),]
    GS.zeta = GS.zeta[-(1:BI),]
    GS.phi = GS.phi[-c(1:BI),]
  }
  SelectR0 = if(is.null(GS.r0)){ 0 }else{ apply(GS.r0, 2, function(t) sum(t!=0))}
  SelectRstar = apply(GS.phi, 2, sum)
  SelectZeta = if(is.null(GS.zeta)){ 0 }else{ apply(GS.zeta, 2, function(t) sum(t!=0))}

  MPM.V = which(SelectRstar > max_BI/2)
  MPM.C = setdiff(which(SelectR0 > max_BI/2), MPM.V)
  MPM.Z = which(SelectZeta > max_BI/2)
  numb = matrix(c(length(MPM.C), length(MPM.V), length(MPM.Z)), ncol=1,
                dimnames=list(c("Constant effect", "Varying effect", "Linear interaction"), "#"))

  Var.names = colnames(obj$coefficient$VC)[-1]
  if(length(MPM.C)>0){
    Main = MPM.C
    names(Main) = Var.names[MPM.C]
  }else{
    Main = NULL
  }

  if(length(MPM.V)>0){
    Varying = MPM.V
    names(Varying) = Var.names[MPM.V]
  }else{
    Varying = NULL
  }

  if(length(MPM.Z)>0){
    Linear = MPM.Z
    names(Linear) = names(obj$coefficient$EG)[MPM.Z]
  }else{
    Linear = NULL
  }

  sel = list(Constant=Main, Varying=Varying, Linear=Linear)
  method = paste("Median Probability Model (MPM)", sep = "")

  out = list(method=method, indices=sel, summary=numb)
  class(out) = "BVSelection"
  out
}
