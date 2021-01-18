\name{altf4}
\alias{altf4}

\title{Computes a Time-Varying Parameters Rolling Regression Averaged over Different Window Sizes.}

\description{
It is necessary to compare a given forecast method with some alternative ones. This function computes selected forecast quality measures for a time-varying parameters rolling regression averaged over different window sizes (which might be treated as alternative forecasting method to Dynamic Model Averaging, Dynamic Model Selection, etc.). The averaging is performed as in Raftery et al. (2010). The only difference is that the state space of the models are constructed not by chosing different combinations of independent variables, but for a fixed set of independent variables various rolling windows sizes are chosen and models constructed in such a way constitute the state space.

ME (Mean Error), RMSE (Root Mean Squared Error), MAE (Mean Absolute Error), MPE (Mean Percentage Errror) and MAPE (Mean Absolute Percentage Error) are computed as \code{\link[forecast]{accuracy}}. HR (Hit Ratio) is computed as \code{\link{hit.ratio}}.
}

\usage{
altf4(y,x,windows,V=NULL,alpha=NULL,lambda=NULL,initial.period=NULL,
d=NULL,fmod=NULL,parallel=NULL,c=NULL,small.c=NULL)
}

\arguments{
\item{y}{\code{\link[base]{numeric}} or a column \code{\link[base]{matrix}} of a dependent variable}
\item{x}{\code{\link[base]{matrix}} of independent variables, different columns correspond to different independent variables}
\item{windows}{\code{\link[base]{numeric}} \code{\link[base]{vector}}, sizes of a rolling regression windows (numbers of observations)}
\item{V}{optional, \code{\link[base]{numeric}}, initial variance in the state space equation for the recursive moment estimator updating method, as in Raftery et al. (2010), if not specified \code{V=1} is taken, see \code{\link{tvp}}}
\item{lambda}{optional, \code{\link[base]{numeric}}, a forgetting factor between 0 and 1 used in variance approximations, if not specified \code{lambda=0.99} is taken, see \code{\link{tvp}}}
\item{alpha}{optional, \code{\link[base]{numeric}}, a forgetting factor \eqn{\alpha} between 0 and 1 used in probabilities estimations, if not specified \code{alpha=0.99} is taken, see \code{\link{fDMA}}}
\item{initial.period}{optional, \code{\link[base]{numeric}}, a number of observation since which forecast quality measures are computed, if not specified the whole sample is used, i.e., \code{initial.period=1}}
\item{d}{optional, \code{\link[base]{logical}}, a parameter used for HR (Hit Ratio) calculation, should be \code{d=FALSE} for level time-series and \code{d=TRUE} if time-series represent changes, if not specified \code{d=FALSE}}
\item{fmod}{optional, class \code{dma} object, a model to be compared with alternative forecast}
\item{parallel}{optional, \code{\link[base]{logical}}, indicate whether parallel computations should be used, by default \code{parallel=FALSE}}
\item{c}{optional, see \code{\link{tvp}}}
\item{small.c}{optional, see \code{\link{fDMA}}}
}

\value{class \code{altf4} object, \code{\link[base]{list}} of
\item{$summary}{\code{\link[base]{matrix}} of forecast quality measures ordered by columns}
\item{$y.hat}{\code{\link[base]{list}} of predicted values from a time-varying parameters rolling regression averaged over selected window sizes}
\item{$y}{\code{y}, forecasted time-series}
\item{$coeff.}{\code{\link[base]{list}} of coefficients from a time-varying parameters rolling regression averaged over selected window sizes}
\item{$weights}{\code{\link[base]{list}} of weights of models used in averaging}
\item{$exp.win.}{\code{\link[base]{list}} of expected window size}
}

\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]

a1 <- altf4(y=ld.wti,x=ld.drivers,d=TRUE,windows=c(36,100,150))

win <- c(36,100,150)
a2 <- altf4(y=ld.wti,x=ld.drivers,d=TRUE,windows=win,alpha=0.9,lambda=0.95)

# models without a constant term
a3 <- altf4(y=ld.wti,x=ld.drivers,d=TRUE,windows=win,alpha=0.9,lambda=0.95,c=FALSE)

# models only with a constant term
empty <- matrix(,nrow=nrow(ld.drivers),ncol=0)
a4 <- altf4(y=ld.wti,x=empty,d=TRUE,windows=win,alpha=0.9,lambda=0.95)
}
}

\references{
Pesaran, M. H., Pick, A., 2011. Forecast combination across estimation windows. \emph{Journal of Business & Economic Statistics} \bold{29}, 307--318.

Raftery, A. E., Gneiting, T., Balabdaoui, F., Polakowski, M., 2005. Using Bayesian Model Averaging to calibrate forecast ensembles. \emph{Monthly Weather Review} \bold{133}, 1155--1174.

Raftery, A. E., Karny, M., Ettler, P., 2010. Online prediction under model uncertainty via Dynamic Model Averaging: Application to a cold rolling mill. \emph{Technometrics} \bold{52}, 52--66.
}

\seealso{\code{\link{plot.altf4}}, \code{\link{print.altf4}}, \code{\link{summary.altf4}}, \code{\link{roll.reg}}, \code{\link{tvp}}, \code{\link{altf}}, \code{\link{altf2}}, \code{\link{altf3}}.}
