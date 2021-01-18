\name{bboptim}
\alias{bboptim}
\alias{rbb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Black-box function optimization
}
\description{
\code{bboptim} implements optimization of a black-box function, possibly estimated with error, using prediction of the function by smoothing of its values in a given set of points, followed by a call to \code{optim} for optimization of the predicted function. \code{rbb} samples the parameter space of the function using a crude implementation of Expected Improvement (e.g. Bingham et al., 2014) methods: points with the highest predicted probability of improvement of the response value among a set of candidates sampled uniformly are retained.
}
\usage{
bboptim(data, ParameterNames = NULL, respName = NULL, control = list(),
        force = FALSE, optimizers = blackbox.getOption("optimizers"), precision=1e-03)
rbb(object,n=NULL,from=NULL,focus=0.75)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  A data frame including both function parameters and function values (or \dQuote{response} values).
}
  \item{ParameterNames}{
    A character vector, identifying the columns of the data that correspond to the function parameters. If NULL, all columns except the last are assumed to hold parameter values.   }
  \item{respName}{
    A character string, identifying the column of the data that corresponds to the function values. If NULL, the last column is assumed to hold function values.   }
  \item{control}{A list passed to the \code{control} argument of the \code{optim} function; e.g., \code{list(fnscale=-1)} for  maximization.}
  \item{force}{Boolean, passed to \code{calcGCV}. TRUE forces the analysis of data without pairs of response values for given parameter values. This is \emph{not} recommended as there \emph{should} be such pairs. If the response is estimated with error, this is required for good smoothing. If it is deterministic, \code{bboptim} will learn it from the information provided by the pairs.}
  \item{optimizers}{(A vector of) character strings, from which the optimization methods are selected. Default are that of \code{\link{calcGCV}} for the smoothing step, and \code{nloptr} with its own \code{"NLOPT_LN_BOBYQA"} method for the smoothed function maximization. See the source of the function for other possible methods (the latter being subject to change with little notice).}
  \item{object}{An object of class \code{bboptim}}
  \item{n}{Number of distinct points to be returned. n+1 points will be returned (see Details). If \code{NULL}, the following default value is used:
  \code{min(2^(np+1),floor(10*(1+3*log(np))))}, where \code{np} is the number of function parameters.
  }
  \item{from}{A larger (>2n) number of points from which \code{n} are selected by an expected inmprovement criterion. If NULL, a default value computed as \code{n*floor(10*(1+3*log(np)))}, where \code{np} is the number of function parameters, is used. }
  \item{focus}{A number between 0 and 1. Determines the proportion of points that are sampled closer to the currently inferred maximum (see Details).}
  \item{precision}{target value of prediction variance in inferred optimum.}
}
\value{
\code{bboptim} returns an object of class \code{bboptim}, a list which includes
\item{optr}{the result of the \code{optim} call}
\item{RMSE}{the root meant square prediction error of response at the optimum \code{optr$par}}
\item{optr_fitted}{the best of the fitted points, with its fitted response value and prediction RMSE}
\item{fit}{the predictor of the response (an \code{HLfit} object as returned by \code{\link[spaMM]{corrHLfit}}, with a \code{predict} method, etc.)}
\item{conv_crits}{Indicators of convergence (see Details)}
and some other elements.

\code{rbb} returns a data frame.
}
\details{ \code{rbb} selects a proportion \code{1-focus} of the returned points according to expected improvement, from points sampled uniformly in a space defined by a tesselation of the fitted \code{object}'s parameter points. They are completed to n-1 points, by points similarly selected but within a space defined by a selection of fitted points with the best predicted response values. Finally, two replicates of the predicted optimum (the \code{optim} \code{$par} result contained in the \code{object}) are included. A total of n+1 points (n distinct) is thus returned.

Global optimization cannot be proven, but it is tested by the following criteria: (1) the predicted optimum is close enough to the optimum among assessed parameter points (i.e. either the optimum parameters are well approached or the function is flat in some way), and (2) the prediction variance at the inferred optimum is low enough (so that the predictions used in the first criterion can be trusted). Accordingly, \code{conv_crits} has elements (1) \code{objective} that indicates whether \code{optr$value} betters \code{optr_fitted$value} by more than \code{control$reltol}, if given, or else by more than \code{sqrt(.Machine$double.eps)}; and (2) \code{precision} that indicates whether variance of prediction error at the inferred optimum is lower than the target \code{precision}. This variance is computed as described for \code{\link[spaMM]{predict.HLfit}}, with \code{variances=list(linPred=TRUE,dispVar=TRUE)}. 
}
\references{
D. Bingham, P. Ranjan, and W.J. Welch (2014) Design of Computer Experiments for Optimization, Estimation of Function Contours, and Related Objectives, pp. 109-124 in Statistics in Action: A Canadian Outlook (J.F. Lawless, ed.). Chapman and Hall/CRC.
}
\examples{
# Classical toy example with optional noise
fr <- function(v,sd) {   ## Rosenbrock Banana function 
  10 * (v["y"] - v["x"]^2)^2 + (1 - v["x"])^2 + rnorm(1,sd=sd)
}
set.seed(123)

# Initial parameter values, including duplicates. See ?init_grid.
parsp <- init_grid(lower=c(x=0,y=0),upper=c(x=2,y=2),nUnique=25)

#### Without noise
# add function values
simuls <- cbind(parsp,bb=apply(parsp,1,"fr",sd=0))

# optimization
bbresu <- bboptim(simuls)
print(bbresu)

# refine with additional points
if (blackbox.getOption("example_maxtime")>4) {
 while ( any( ! bbresu$conv_crits) ) {
  print(unlist(bbresu$optr[c("par","value")]))
  candidates <- rbb(bbresu)
  newsimuls <- cbind(candidates,bb=apply(candidates,1,"fr",sd=0))
  bbresu <- bboptim(rbind(bbresu$fit$data,newsimuls))
 }
 print(bbresu)
}

#### With noise

if (blackbox.getOption("example_maxtime")>78) {
 set.seed(123)
 simuls <- cbind(parsp,bb=apply(parsp,1,"fr",sd=0.1))

 bbresu <- bboptim(simuls, precision=0.02)
 
 while ( any( ! bbresu$conv_crits) ) {
  print(unlist(bbresu$optr[c("par","value")]))
  candidates <- rbb(bbresu)
  newsimuls <- cbind(candidates,bb=apply(candidates,1,"fr",sd=0.1))
  bbresu <- bboptim(rbind(bbresu$fit$data,newsimuls), precision=0.02)
 }
 print(bbresu)
}

# basic plot
\dontrun{
require(spaMM)
opt <- bbresu$optr$par
mapMM(bbresu$fit, decorations=points(opt[1],opt[2],cex=2,pch="+"))
}
}
\keyword{optimize}