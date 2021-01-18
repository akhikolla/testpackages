\name{mixest2}
\alias{mixest2}

\title{Computes Mixture Estimation with Normal Regression Components.}

\description{
This function estimates recursively mixtures with normal regression components with a dynamic model of switching.
}

\usage{
mixest2(y,x,mods=NULL,ftype=NULL,V=NULL,W=NULL,atype=NULL,Tvar=NULL)
}

\arguments{
\item{y}{one column {\code{\link[base]{matrix}}} of forecasted time-series, observations inserted rowwise}
\item{x}{{\code{\link[base]{matrix}}} of independent time-series (predictors), observations inserted rowwise}
\item{mods}{see {\code{\link[dynmix]{mixest1}}}}
\item{ftype}{optional, \code{\link[base]{numeric}} indicating type of forecasting, 1 represents averaging forecasts from all components by the estimated weights (i.e., forecasting based on coefficients derived from the estimated mixture), 2 represents selecting the forecast given by the model with the highest weight, 3 represents selecting the forecast from the so-called median probability model (Barbieri and Berger, 2004), by default \code{ftype=1} is taken}
\item{V}{optional, \code{\link[base]{numeric}} initial variance, the same for all components, by default \code{V=1} is taken}
\item{W}{optional, \code{\link[base]{numeric}} initial value to be put in the diagonal  matrix representing the covariance matrices for regression coefficients, the same for all components, by default \code{W=1} is taken}
\item{atype}{optional, \code{\link[base]{numeric}} indicating approximation of pdfs, 0 represents quasi-Bayesian approach, 1 represents minimization of the Kerridge inaccuracy, by default \code{atype=0} is taken}
\item{Tvar}{optional, \code{\link[base]{numeric}} indicating the number of first observations, when variance and covariance updating will not be performed, i.e., the initial values will be kept fixed, by default \code{Tvar=30} is taken}
}

\value{object of class {\code{mixest}}, i.e., {\code{\link[base]{list}}} of 
\item{$y.hat}{{\code{\link[base]{vector}}} of predictions}
\item{$rvi}{{\code{\link[base]{matrix}}} of relative variable importances}
\item{$coef}{{\code{\link[base]{matrix}}} of regression coefficients corresponding to \code{ftype} method chosen}
\item{$weights}{{\code{\link[base]{matrix}}} of estimated weights of component models}
\item{$V}{{\code{\link[base]{vector}}} of updated variances from the selected models, consistent with \code{ftype} chosen}
\item{$R}{{\code{\link[base]{matrix}}} of updated diagonal of covariances corresponding to independent variables in regressions, consistent with \code{ftype} chosen}
\item{$components}{{\code{\link[base]{matrix}}} of \code{mods}}
\item{$parameters}{{\code{\link[base]{character}}} of parameters used in the model}
}

\examples{
\donttest{
data(oil)
m1 <- mixest2(y=oil[,1,drop=FALSE],x=oil[,-1,drop=FALSE],ftype=1,V=100,W=100)
}
}

\source{
Nagy, I., Suzdaleva, E., Karny, M., Mlynarova, T., 2011, Bayesian estimation of dynamic finite mixtures. \emph{International Journal of Adaptive Control and Signal Processing} \bold{25}, 765--787.
}

\references{
Barbieri, M. M., Berger, J. O., 2004, Optimal predictive model selection. \emph{The Annals of Statistics} \bold{32}, 870--897.

Burnham, K. P., Anderson, D. R., 2002, \emph{Model Selection and Multimodel Inference}, Springer. 

Dedecius, K., 2010, \emph{Partial Forgetting in Bayesian Estimation}, Czech Technical University in Prague. 

Karny, M. (ed.), 2006, \emph{Optimized Bayesian Dynamic Advising}, Springer. 

Nagy, I., 2015, \emph{Mixture Models and Their Applications}, Czech Technical University in Prague. 

Nagy, I., Suzdaleva, E., 2017, \emph{Algorithms and Programs of Dynamic Mixture Estimation}, Springer. 
}

\seealso{\code{\link[dynmix]{mixest1}}}
