#' Variable selection for a roben object
#'
#' Variable selection for a roben object
#'
#' @param obj roben object.
#' @param ... other GxESelection arguments.
#'
#' @details For class `Sparse', the median probability model (MPM) (Barbieri and Berger, 2004) is used to identify predictors that are significantly associated
#' with the response variable. For class `NonSparse', variable selection is based on 95\% credible interval.
#' Please check the references for more details about the variable selection.
#'
#' @references
#' Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y. and Wu, C. (2020). Robust Bayesian variable selection for gene-environment interactions.
#'
#' Barbieri, M.M. and Berger, J.O. (2004). Optimal predictive model selection. {\emph{Ann. Statist}, 32(3):870–897}
#'
#' @rdname GxESelection
#' @return an object of class `GxESelection' is returned, which is a list with components:
#' \item{method}{method used for identifying important effects.}
#' \item{effects}{a list of names of selected effects.}
#' \item{summary}{a summary of selected effects.}
#' \item{indicator}{a matrix of indicators of selected effects.}
#'
#' @seealso \code{\link{roben}}
#'
#' @examples
#' data(GxE_small)
#' iter=5000
#' ## sparse
#' fit=roben(X, Y, E, clin, iterations=iter)
#' selected=GxESelection(fit)
#' selected
#'
#' \donttest{
#' ## non-sparse
#' fit=roben(X, Y, E, clin, iterations=iter, sparse=FALSE)
#' selected=GxESelection(fit)
#' selected
#' }
#'
#' @export
GxESelection <- function(obj,...){
  if(!inherits(obj, "roben")) stop("This is not a roben object")
  UseMethod('GxESelection', obj)
}


#' @rdname GxESelection
#' @method GxESelection Sparse
#' @export
GxESelection.Sparse=function(obj, burn.in=obj$burn.in,...){
  BI = ifelse(is.null(burn.in), 0, burn.in)
  max_BI = obj$iterations - BI

  GS.beta = obj$posterior$GS.beta
  s = ncol(obj$coefficient$GE)
  if(BI>0){
    GS.beta = GS.beta[-c(1:BI),]
  }
  SelectBeta = apply(GS.beta, 2, function(t) sum(t!=0))


  Ind.GE = matrix((SelectBeta > max_BI/2)*1, ncol=s)

  numb = matrix(c(sum(Ind.GE[1,]), sum(Ind.GE[-1,])), ncol=1,
                dimnames=list(c("Main G effects", "GxE interactions"), "#"))

  G.names = colnames(obj$coefficient$GE)
  E.names = rownames(obj$coefficient$GE)

  colnames(Ind.GE) = G.names
  rownames(Ind.GE) = E.names

  if(sum(Ind.GE[1,])>0){
    Main = G.names[which(Ind.GE[1,]>0)]
  }else{
    Main = NULL
  }

  if(sum(Ind.GE[-1,])>0){
    inds = which(Ind.GE>0, arr.ind = TRUE)
    inds = inds[inds[,1]!=1,]
    paste(G.names[inds[,2]], "x", E.names[inds[,1]], sep="")
    GxE = paste(G.names[inds[,2]], "x", E.names[inds[,1]], sep="")
  }else{
    GxE = NULL
  }


  sel = list(Main.G=Main, GxE=GxE)
  method = paste("Median Probability Model (MPM)", sep = "")

  out = list(method=method, effects=sel, summary=numb, indicator=Ind.GE)
  class(out) = "GxESelection"
  out
}

#' @param burn.in MCMC burn-in.
#' @param prob probability for credible interval, between 0 and 1. e.g. prob=0.95 leads to 95\% credible interval.
#' @rdname GxESelection
#' @method GxESelection NonSparse
#' @export
GxESelection.NonSparse=function(obj, burn.in=obj$burn.in, prob=0.95,...){
  BI = ifelse(is.null(burn.in), 0, burn.in)

  GS.beta = obj$posterior$GS.beta
  s = ncol(obj$coefficient$GE)

  if(BI>0){
    GS.beta = GS.beta[-c(1:BI),]
  }

  s = ncol(obj$coefficient$GE)

  lt = (1-prob)/2; ut= 1-lt
  limits = apply(GS.beta, 2, stats::quantile, probs = c(lt, ut))
  temp = matrix(abs(sign(limits[1,]) + sign(limits[2,]))==2, ncol = s)
  Ind.GE = temp*1

  numb = matrix(c(sum(Ind.GE[1,]), sum(Ind.GE[-1,])), ncol=1,
                dimnames=list(c("Main G effects", "GxE interactions"), "#"))

  G.names = colnames(obj$coefficient$GE)
  E.names = rownames(obj$coefficient$GE)

  colnames(Ind.GE) = G.names
  rownames(Ind.GE) = E.names

  if(sum(Ind.GE[1,])>0){
    Main = G.names[which(Ind.GE[1,]>0)]
  }else{
    Main = NULL
  }

  if(sum(Ind.GE[-1,])>0){
    inds = which(Ind.GE>0, arr.ind = TRUE)
    inds = inds[inds[,1]!=1,]
    paste(G.names[inds[,2]], "x", E.names[inds[,1]], sep="")
    GxE = paste(G.names[inds[,2]], "x", E.names[inds[,1]], sep="")
  }else{
    GxE = NULL
  }


  sel = list(Main.G=Main, GxE=GxE)
  method = paste(prob*100,"% credible interval", sep = "")

  out = list(method=method, effects=sel, summary=numb, indicator=Ind.GE)
  class(out) = "GxESelection"
  out

}

