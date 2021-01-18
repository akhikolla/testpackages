\name{altf}
\alias{altf}

\title{Computes a Few Alternative Forecasts.}

\description{
It is necessary to compare a given forecast method with some alternative ones. This function computes selected forecast quality measures for a few selected forecast methods (which might be treated as alternative ones to Dynamic Model Averaging, Dynamic Model Selection, etc.). 

Naive forecast (naive) is computed in a way that all forecasts are set to be the value of the last observation.

For rolling OLS forecast (roll. OLS) for the first periods (until the size of a \code{window} is obtained) are estimated through recursive OLS (rec. OLS). 

Autoregressive models (AR(1) and AR(2)) are computed by ordinary least squares method. 

Time-varying parameters models (TVP, TVP-AR(1) and TVP-AR(2)) are computed as \code{\link{tvp}} with \code{V=1} and \code{lambda=0.99}.

Auto ARIMA (auto ARIMA) is computed as \code{\link[forecast]{auto.arima}}.

ME (Mean Error), RMSE (Root Mean Squared Error), MAE (Mean Absolute Error), MPE (Mean Percentage Errror) and MAPE (Mean Absolute Percentage Error) are computed as \code{\link[forecast]{accuracy}}. HR (Hit Ratio) is computed as \code{\link{hit.ratio}}.
}

\usage{
altf(y,x,window=NULL,initial.period=NULL,d=NULL,f=NULL,fmod=NULL,c=NULL)
}

\arguments{
\item{y}{\code{\link[base]{numeric}} or a column \code{\link[base]{matrix}} of a dependent variable}
\item{x}{\code{\link[base]{matrix}} of independent variables, different columns correspond to different independent variables}
\item{window}{optional, \code{\link[base]{numeric}}, a size of a rolling regression window (a number of observations), if not specified 10\% of all observations are taken}
\item{initial.period}{optional, \code{\link[base]{numeric}}, a number of observation since which forecast quality measures are computed, if not specified the whole sample is used, i.e., \code{initial.period=1}, this argument also divides the sample into in-sample and out-of-sample for non-recursive methods (OLS, AR(1), AR(2), auto ARIMA)}
\item{d}{optional, \code{\link[base]{logical}}, a parameter used for HR (Hit Ratio) calculation, should be \code{d=FALSE} for level time-series and \code{d=TRUE} if time-series represent changes, if not specified \code{d=FALSE}}
\item{f}{optional, \code{\link[base]{logical}} vector, indicating which of alternative forecasts -- naive, OLS, rec. OLS, roll. OLS, TVP, AR(1), AR(2), auto ARIMA, TVP-AR(1) and TVP-AR(2) -- should be computed, if not specified \code{f=c(\link[base]{rep}(TRUE,10))}, i.e., all alternative forecasts are computed}
\item{fmod}{optional, class \code{dma} object, a model to be compared with alternative forecast}
\item{c}{optional, \code{\link[base]{logical}}, a parameter indicating whether constant is included in models, if not specified \code{c=TRUE} is used, i.e., constant is included}
}

\value{class \code{altf} object, \code{\link[base]{list}} of
\item{$summary}{\code{\link[base]{matrix}} of forecast quality measures ordered by columns, forecast methods are ordered by rows}
\item{$y.hat}{\code{\link[base]{list}} of predicted values from all forecasting methods which were applied}
\item{$y}{\code{y}, forecasted time-series}
\item{$coeff.}{\code{\link[base]{list}} of coefficients from all forecasting methods which were applied (for naive forecast they are not computed)}
\item{$p.val.}{\code{\link[base]{list}} of p-values for t-test of statistical significance for coefficients from all forecasting methods which were applied (for naive and TVP models they are not computed, and for auto ARIMA z-test is used)}
}

\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]

a1 <- altf(y=ld.wti,x=ld.drivers,d=TRUE,initial.period=60)

# models where constant term is not included in modelled equations (if applicable)
a2 <- altf(y=ld.wti,x=ld.drivers,d=TRUE,c=FALSE,initial.period=60)

# compute just selected models
fcomp <- c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,FALSE)
a3 <- altf(y=ld.wti,x=ld.drivers,d=TRUE,f=fcomp,initial.period=60)

m1 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10)
a4 <- altf(y=ld.wti,x=ld.drivers,d=TRUE,f=fcomp,fmod=m1,initial.period=60)
}
}

\seealso{\code{\link{plot.altf}}, \code{\link{print.altf}}, \code{\link{summary.altf}}, \code{\link{rec.reg}}, \code{\link{roll.reg}}, \code{\link{altf2}}, \code{\link{altf3}}, \code{\link{altf4}}.}
