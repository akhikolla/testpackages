
\name{predgc}
\alias{predgc}
\title{Prediction at Unobserved Locations in Gaussian Copula Models for Geostatistical Count Data}
\usage{
predgc(obs.y, obs.x = NULL, obs.locs, pred.x = NULL, pred.locs,
       longlat = FALSE, distscale = 1, marginal, corr, obs.effort = 1,
       pred.effort = 1, method = "GHK", estpar = NULL, corrpar0 = NULL,
       pred.interval = NULL, parallel = FALSE,
       ghkoptions = list(nrep = c(100,1000), reorder = FALSE, seed = 12345),
       paralleloptions = list(n.cores = 2, cluster.type = "SOCK"))
}
\arguments{
\item{obs.y}{a non-negative integer vector of observed response with its length
equals to the number of observed locations.}

\item{obs.x}{a numeric matrix or data frame of covariates at observed locations,
with its number of rows equals to the number of observed locations.
If no covariates then \code{obs.x = NULL}.}

\item{obs.locs}{a numeric matrix or data frame of observed locations.obs.effort
The first column is \emph{x} or longitude, the second column is \emph{y} or latitude.
The number of observed locations is equal to the number of rows.}

\item{pred.x}{a numeric matrix or data frame of covariates at prediction locations,
with its number of rows equals to the number of prediction locations.
If no covariates then \code{pred.x = NULL}.}

\item{pred.locs}{a numeric matrix or data frame of prediction locations.
First column is \emph{x} or longitude, second column is \emph{y} or latitude.
The number of prediction locations equals to the number of rows.}

\item{longlat}{if FALSE, use Euclidean distance, if TRUE use great circle distance. The default is FALSE.}

\item{distscale}{a numeric scaling factor for computing distance. If original distance is in kilometers, then
\code{distscale = 1000} will convert it to meters.}

\item{marginal}{an object of class \code{\link{marginal.gc}} specifying the marginal distribution.}

\item{corr}{an object of class \code{\link{corr.gc}} specifying the correlation function.}

\item{obs.effort}{sampling effort at observed locations.
For binomial marginal it is the size parameter (number of trials). See details.}

\item{pred.effort}{sampling effort at prediction locations.
For binomial marginal it is the size parameter (number of trials). See details.}

\item{method}{two methods are implemented. If
\code{method = 'GHK'} then the maximum simulated likelihood estimates are computed, if
\code{method = 'GQT'} then the maximum surrogate likelihood estimates are computed.}

\item{estpar}{if \code{estpar = NULL} then the likelihood estimates
will be computed first, then plug-in into the predictive density.
When all estimates are available,
it is suggested to specify \code{estpar} with the parameter estimates so
re-fitting is not needed. If so, the sequence of the
input values for \code{estpar} is: regression parameters, overdispersion (if any), and correlation parameters (range and nugget, if applicable).
}

\item{corrpar0}{the starting value of correlation parameters in optimization procedure.
        If \code{corrpar0 = NULL} then
initial range is set to be half of the median distance in distance matrix, and initial nugget
(if \code{nugget = TRUE}) is 0.2.}

\item{pred.interval}{a number between \code{0} and \code{1} representing confidence level of the prediction
interval. The program will output two types of the prediction intervals, see detail.
If \code{pred.interval = NULL} then no prediction interval will be computed.}

\item{parallel}{if \code{TRUE} then parallel computing is used to predict multiple prediction locations
simultaneously. If \code{FALSE} then a serial version will be called.}

\item{ghkoptions}{a list of three elements that only need to be specified if \code{method = 'GHK'}.

       \code{nrep} is the Monte Carlo size of the importance sampling algorithm for likelihood approximation.
        It can be a vector with increasing positive integers so that the model is
        fitted with a sequence of different Monte Carlo sizes, and the starting
         values for optimization are taken from the previous fitting.
          The default value is 100 for the first optimization and 1000 for the second and definitive optimization.

       \code{reorder} indicates whether the integral will be reordered every iteration in computation
       according to the algorithm in Gibson, etal (1994), default is FALSE.

       \code{seed} is seed of the pseudorandom generator used in Monte Carlo simulation.}

\item{paralleloptions}{a list of two elements that only need to be specified if \code{parallel = TRUE}.

     \code{n.cores} is the number of cores to be used in parallel prediction.

     \code{cluster.type} is type of cluster to be used for parallel computing; can be "SOCK", "MPI", "PVM", or "NWS".}
}
\value{
A list of class "predgc" with the following elements:
\item{obs.locs}{ observed locations.}

\item{obs.y}{ observed values at observed locations.}

\item{pred.locs}{ prediction locations.}

\item{predValue}{the expectation of the conditional predictive distribution.}

\item{predCount}{ predicted counts; the closest integer that \code{predValue} rounded to.}

\item{predVar}{ estimated variance of the prediction at prediction locations.}

\item{ConfidenceLevel}{ confidence level (between 0 to 1) if prediction interval is computed.}

\item{predInterval.EqualTail}{ equal-tail prediction interval.}

\item{predInterval.Shortest}{ shortest length prediction interval.}
}
\description{
Computes the plug-in prediction at unobserved sites. Two methods are implemented. If
\code{method = 'GHK'} then the maximum simulated likelihood estimates are computed and
the sequential importance sampling method is used in the integral evaluation. If
\code{method = 'GQT'} then the maximum surrogate likelihood estimates are computed and the
generalized quantile transform method is used in integral approximation.
}
\details{
This program implemented two methods in predicting the response at unobserved sites. See \code{\link{mlegc}}.

         The argument \code{obs.effort} and \code{pred.effort}
         are the sampling effort (known). It can be used to consider heterogeneity of
         the measurement time or area at different locations. The default is 1 for all locations.
         See Han and De Oliveira (2016) for more details.

         The program computes two types of prediction intervals at a given confidence level. The shortest
         prediction interval is obtained from evaluating the highest to lowest prediction densities;
         the equal tail prediction interval has equal tail probabilities.
}
\examples{
\dontrun{
## For fast check predict at four locations only
data(Weed95)
weedobs <- Weed95[Weed95$dummy==1, ]
weedpred <- Weed95[Weed95$dummy==0, ]
predweed1 <- predgc(obs.y = weedobs$weedcount, obs.x = weedobs[,4:5], obs.locs = weedobs[,1:2],
                     pred.x = weedpred[1:4,4:5], pred.locs = weedpred[1:4,1:2],
                     marginal = negbin.gc(link = 'log'), pred.interval = 0.9,
                     corr = matern.gc(kappa = 0.5, nugget = TRUE), method = 'GHK')
#summary(predweed1)
#plot(predweed1)

## Time consuming examples
## Weed prediction at 200 locations using parallel programming
predweed2 <- predgc(obs.y = weedobs$weedcount, obs.x = weedobs[,4:5], obs.locs = weedobs[,1:2],
                     pred.x = weedpred[,4:5], pred.locs = weedpred[,1:2],
                     marginal = negbin.gc(link = 'log'),
                     corr = matern.gc(kappa = 0.5, nugget = TRUE), method = 'GHK',
                     pred.interval = 0.95, parallel = TRUE,
                     paralleloptions = list(n.cores = 4))
#summary(predweed2)
#plot(predweed2)


## A more time consuming example for generating a prediction map at a fine grid
data(OilWell)
gridstep <- seq(0.5, 30.5, length = 40)
locOilpred <- data.frame(Easting = expand.grid(gridstep, gridstep)[,1],
                        Northing = expand.grid(gridstep, gridstep)[,2])
PredOil <- predgc(obs.y = OilWell[,3], obs.locs = OilWell[,1:2],  pred.locs = locOilpred,
                       marginal = binomial.gc(link = 'logit'),
                       corr = matern.gc(nugget = FALSE), obs.effort = 1,
                       pred.effort = 1, method = 'GHK',
                       parallel = TRUE, paralleloptions = list(n.cores = 4))
PredMat <- summary(PredOil)

## To generate better prediction maps
library(colorspace)
filled.contour(seq(0.5,30.5,length=40), seq(0.5,30.5,length=40),
               matrix(PredMat$predMean,40,), zlim = c(0, 1), col=rev(heat_hcl(12)),
               nlevels=12, xlab = "Eastings", ylab = "Northings",
               plot.axes = {axis(1); axis(2); points(OilWell[,1:2], col = 1,
               cex = 0.25 + 0.25*OilWell[,3])})

filled.contour(seq(0.5,30.5,length=40), seq(0.5,30.5,length=40),
               matrix(PredMat$predVar,40,),
               zlim = c(0, 0.3), col = rev(heat_hcl(12)), nlevels = 10,
               xlab = "Eastings", ylab = "Northings",
               plot.axes = {axis(1); axis(2); points(OilWell[,1:2], col = 1,
               cex = 0.25 + 0.25*OilWell[,3])})

}

}
\author{
Zifei Han \email{hanzifei1@gmail.com}
}
\references{
Han, Z. and De Oliveira, V. (2016) On the correlation structure of Gaussian
            copula models for geostatistical count data.
            \emph{Australian and New Zealand Journal of Statistics, 58:47-69}.

Kazianka, H. and Pilz, J. (2010)
Copula-based geostatistical modeling of continuous and discrete data including covariates.
\emph{Stoch Environ Res Risk Assess 24:661-673}.

Kazianka, H. (2013)
Approximate copula-based estimation and prediction of discrete spatial data.
\emph{Stoch Environ Res Risk Assess 27:2015-2026}.

 Masarotto, G. and Varin, C. (2012)
 Gaussian copula marginal regression.
 \emph{Electronic Journal of Statistics 6:1517-1549}.
 \url{https://projecteuclid.org/euclid.ejs/1346421603}.

 Masarotto, G. and Varin C. (2017). Gaussian Copula Regression in R. \emph{Journal of Statistical Software}, \bold{77}(8), 1--26. \doi{10.18637/jss.v077.i08}.

Han, Z. and De Oliveira, V. (2018) gcKrig: An R Package for the Analysis of Geostatistical Count Data Using Gaussian Copulas.
\emph{Journal of Statistical Software}, \bold{87}(13), 1--32. \doi{10.18637/jss.v087.i13}.
}
\seealso{
\code{\link[gcmr]{gcmr}};
         \code{\link{mlegc}}
}

\keyword{Prediction}
