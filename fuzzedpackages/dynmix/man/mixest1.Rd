\name{mixest1}
\alias{mixest1}

\title{Computes Mixture Estimation with State-Space Components.}

\description{
This function estimates recursively mixtures with state-space components with a dynamic model of switching. The components are normal linear models. Suppose there are available \eqn{k} potentially important predictors of \code{y}, i.e., \eqn{x_1, \dots, x_k}. Then up to \eqn{2^{k}} linear models including constant term can be created by inclding or not including each of these predictors in the individual model, i.e., component of the mixture.
}

\usage{
mixest1(y,x,mods=NULL,ftype=NULL,lambda=NULL,kappa=NULL,V=NULL,W=NULL,atype=NULL)
}

\arguments{
\item{y}{one column {\code{\link[base]{matrix}}} of forecasted time-series, observations inserted rowwise}
\item{x}{{\code{\link[base]{matrix}}} of independent time-series (predictors), observations inserted rowwise}
\item{mods}{optional, \code{\link[base]{matrix}} indicating which models should be used as components, the first column indicates inclusion of a constant in a component model, by default all possible models with a constant are used, inclusion of a variable is indicated by 1, omitting by 0, component models are indexed by rows, variables (time-series) are indexed by columns}
\item{ftype}{optional, \code{\link[base]{numeric}} indicating type of forecasting, 0 represents forecasting based on coefficients derived from the estimated mixture, 1 represents averaging forecasts from all components by the estimated weights, 2 represents selecting the forecast given by the model with the highest weight, 3 represents selecting the forecast from the so-called median probability model (Barbieri and Berger, 2004), by default \code{ftype=0} is taken}
\item{lambda}{optional, \code{\link[base]{numeric}} between 0 and 1, a forgetting factor in covariance estimation method described by Raftery et al. (2010), by default the method of Nagy and Suzdaleva (2013) is used}
\item{kappa}{optional, \code{\link[base]{numeric}} between 0 and 1, a parameter for the exponentially weighted moving average estimation of components variances, described for example by Koop and Korobilis (2012), if \code{lambda} is specified but \code{kappa} is not, then the method of recursive moments described by Raftery et al. (2010) is used, by default the method of Nagy and Suzdaleva (2013) is used}
\item{V}{optional, \code{\link[base]{numeric}} initial variance for all components (output equation), by default \code{V=1} is taken}
\item{W}{optional, \code{\link[base]{numeric}} initial value to be put in the diagonal matrix representing the covariance matrices (state equation), by default \code{W=1} is taken}
\item{atype}{optional, \code{\link[base]{numeric}} indicating approximation of pdf, 0 represents quasi-Bayesian approach, 1 represents minimization of the Kerridge inaccuracy (where suitiable optimization is done with the Gauss-Newton method, still this increases the computation time greatly), by default \code{atype=0} is taken}
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
m1 <- mixest1(y=oil[,1,drop=FALSE],x=oil[,-1,drop=FALSE],ftype=1,V=100,W=100)

# Models with only one variable
mods <- diag(1,nrow=ncol(oil[,-1,drop=FALSE]),ncol=ncol(oil[,-1,drop=FALSE]))
mods <- cbind(1,mods)
m2 <- mixest1(y=oil[,1,drop=FALSE],x=oil[,-1,drop=FALSE],mods=mods,ftype=1,V=100,W=100)
}
}

\source{
Nagy, I., Suzdaleva, E., 2013, Mixture estimation with state-space components and Markov model of switching. \emph{Applied Mathematical Modelling} \bold{37}, 9970--9984.
}

\references{
Barbieri, M. M., Berger, J. O., 2004, Optimal predictive model selection. \emph{The Annals of Statistics} \bold{32}, 870--897.

Burnham, K. P., Anderson, D. R., 2002, \emph{Model Selection and Multimodel Inference}, Springer. 

Karny, M. (ed.), 2006, \emph{Optimized Bayesian Dynamic Advising}, Springer. 

Koop, G., Korobilis, D., 2012, Forecasting inflation using Dynamic Model Averaging. \emph{International Economic Review} \bold{53}, 867--886.

Nagy, I., Suzdaleva, E., 2017, \emph{Algorithms and Programs of Dynamic Mixture Estimation}, Springer. 

Raftery, A. E., Karny, M., Ettler, P., 2010, Online prediction under model uncertainty via Dynamic Model Averaging: Application to a cold rolling mill. \emph{Technometrics} \bold{52}, 52--66.
}

\seealso{\code{\link[dynmix]{mixest2}}}
