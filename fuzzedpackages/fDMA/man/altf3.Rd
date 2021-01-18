\name{altf3}
\alias{altf3}

\title{Computes a Rolling Regression Averaged over Different Window Sizes.}

\description{
It is necessary to compare a given forecast method with some alternative ones. This function computes selected forecast quality measures for a rolling regression averaged over different window sizes (which might be treated as alternative forecasting method to Dynamic Model Averaging, Dynamic Model Selection, etc.). 

ME (Mean Error), RMSE (Root Mean Squared Error), MAE (Mean Absolute Error), MPE (Mean Percentage Errror) and MAPE (Mean Absolute Percentage Error) are computed as \code{\link[forecast]{accuracy}}. HR (Hit Ratio) is computed as \code{\link{hit.ratio}}.
}

\usage{
altf3(y,x=NULL,windows,av=NULL,initial.period=NULL,d=NULL,fmod=NULL,parallel=NULL,c=NULL)
}

\arguments{
\item{y}{\code{\link[base]{numeric}} or a column \code{\link[base]{matrix}} of a dependent variable}
\item{x}{\code{\link[base]{matrix}} of independent variables, different columns correspond to different independent variables, if not specified only constant term will be included}
\item{windows}{\code{\link[base]{numeric}} \code{\link[base]{vector}}, sizes of a rolling regression windows (numbers of observations)}
\item{av}{optional, a method for model averaging, \code{av="ord"} corresponds to equal weights for each model, \code{av="aic"} corresponds to information theoretic model averaging based on Akaike Information Criterion, \code{av="aicc"} corresponds to information theoretic model averaging based on Akaike Information Criterion with a correction for finite sample sizes, \code{av="bic"} corresponds to information theoretic model averaging based on Bayesian Information Criterion, \code{av="mse"} corresponds to setting weights proportional to the inverse of the models Mean Squared Error, if \code{av} is \code{\link[base]{numeric}} then weights are computed proportional to the \code{av}-th power of window size, if not specified \code{av="ord"} is used}
\item{initial.period}{optional, \code{\link[base]{numeric}}, a number of observation since which forecast quality measures are computed, if not specified the whole sample is used, i.e., \code{initial.period=1}}
\item{d}{optional, \code{\link[base]{logical}}, a parameter used for HR (Hit Ratio) calculation, should be \code{d=FALSE} for level time-series and \code{d=TRUE} if time-series represent changes, if not specified \code{d=FALSE}}
\item{fmod}{optional, class \code{dma} object, a model to be compared with alternative forecast}
\item{parallel}{optional, \code{\link[base]{logical}}, indicate whether parallel computations should be used, by default \code{parallel=FALSE}}
\item{c}{optional, see \code{\link{roll.reg}}}
}

\details{For each \code{av} method, in the initial period equal weights for each model are taken, and then successively updated based on the chosen criterion.}

\value{class \code{altf3} object, \code{\link[base]{list}} of
\item{$summary}{\code{\link[base]{matrix}} of forecast quality measures ordered by columns}
\item{$y.hat}{\code{\link[base]{list}} of predicted values from a rolling regression averaged over selected window sizes}
\item{$y}{\code{y}, forecasted time-series}
\item{$coeff.}{\code{\link[base]{list}} of coefficients from a rolling regression averaged over selected window sizes}
\item{$weights}{\code{\link[base]{list}} of weights of models used in averaging}
\item{$p.val.}{\code{\link[base]{list}} of p-values (averaged over selected window sizes) for t-test of statistical significance for coefficients from a rolling regression}
\item{$exp.win.}{\code{\link[base]{list}} of expected window size}
}

\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]

a1 <- altf3(y=ld.wti,x=ld.drivers,d=TRUE,windows=c(36,100,150))

a2 <- altf3(y=ld.wti,x=ld.drivers,d=TRUE,av="aic",windows=c(36,100,150))

a3 <- altf3(y=ld.wti,x=ld.drivers,d=TRUE,av=-2,windows=c(36,100,150))

# models without a constant term 
a4 <- altf3(y=ld.wti,x=ld.drivers,d=TRUE,av=-2,windows=c(36,100,150),c=FALSE)

# models only with a constant term
a5 <- altf3(y=ld.wti,d=TRUE,av=-2,windows=c(36,100,150))
}
}

\references{
Pesaran, M. H., Pick, A., 2011. Forecast combination across estimation windows. \emph{Journal of Business & Economic Statistics} \bold{29}, 307--318. 
}

\seealso{\code{\link{plot.altf3}}, \code{\link{print.altf3}}, \code{\link{summary.altf3}}, \code{\link{roll.reg}}, \code{\link{altf}}, \code{\link{altf2}}, \code{\link{altf4}}.}
