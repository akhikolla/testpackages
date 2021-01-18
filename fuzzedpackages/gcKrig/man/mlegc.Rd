
\name{mlegc}
\alias{mlegc}
\title{Maximum Likelihood Estimation in Gaussian Copula Models for Geostatistical Count Data}
\usage{
mlegc(y, x = NULL, locs, marginal, corr, effort = 1, longlat = FALSE,
  distscale = 1, method = "GHK", corrpar0 = NULL, ghkoptions = list(nrep
  = c(100, 1000), reorder = FALSE, seed = 12345))
}
\arguments{
\item{y}{a non-negative integer vector of response with its length equals to the number of sampling locations.}

\item{x}{a numeric matrix or data frame of covariates,
with its number of rows equals to the number of sampling locations.
If no covariates then \code{x = NULL}.}

\item{locs}{a numeric matrix or data frame of \emph{n-D} points with row denoting points.
The first column is \emph{x} or longitude, the second column is \emph{y} or latitude.
The number of locations is equal to the number of rows.}

\item{marginal}{an object of class \code{\link{marginal.gc}} specifying the marginal distribution.}

\item{corr}{an object of class \code{\link{corr.gc}} specifying the correlation function.}

\item{effort}{the sampling effort. For binomial marginal it is the size parameter (number of trials).
See details.}

\item{longlat}{if FALSE, use Euclidean distance, if TRUE use great circle distance. The default is FALSE.}

\item{distscale}{a numeric scaling factor for computing distance. If original distance is in kilometers, then
\code{distscale = 1000} will convert it to meters.}

\item{method}{two methods are implemented. If
\code{method = 'GHK'} then the maximum simulated likelihood estimates are computed, if
\code{method = 'GQT'} then the maximum surrogate likelihood estimates are computed.}

\item{corrpar0}{the starting value of correlation parameter in the optimization procedure.
        If \code{corrpar0 = NULL} then
initial range is set to be half of the median distance in distance matrix and initial nugget
(if \code{nugget = TRUE}) is 0.2.}

\item{ghkoptions}{a list of three elements that only need to be specified if \code{method = 'GHK'}.

       \code{nrep} is the Monte Carlo size of the importance sampling algorithm for likelihood approximation.
        It can be a vector with increasing positive integers so that the model is
        fitted with a sequence of different Monte Carlo sizes, and the starting
         values for optimization are taken from the previous fitting.
          The default value is 100 for the first optimization and 1000 for the second and definitive optimization.

       \code{reorder} indicates whether the integral will be reordered every iteration in computation
       according to the algorithm in Gibson, etal (1994), default is FALSE.

       \code{seed}  is the seed of the pseudorandom generator used in Monte Carlo simulation.}
}
\value{
A list of class "mlegc" with the following elements:

\item{MLE}{the maximum likelihood estimate.}
\item{x}{the design matrix.}
\item{nug}{1 if \code{nugget = TRUE}, 0 if \code{nugget = FALSE}.}
\item{nreg}{number of regression parameters.}
\item{log.lik}{the value of the maximum log-likelihood.}
\item{AIC}{the Akaike information criterion.}
\item{AICc}{the AICc information criterion; essentially AIC with a greater penalty for extra parameters.}
\item{BIC}{ the Bayesian information criterion.}
\item{kmarg}{number of marginal parameters.}
\item{par.df}{number of parameters.}
\item{N}{number of observations.}
\item{D}{the distance matrix.}
\item{optlb}{lower bound in optimization.}
\item{optub}{upper bound in optimization.}
\item{hessian}{the hessian matrix evaluated at the final estimates.}
\item{args}{arguments passed in function evaluation.}
}
\description{
Computes the  maximum likelihood estimates. Two methods are implemented. If
\code{method = 'GHK'} then the maximum simulated likelihood estimates are computed, if
\code{method = 'GQT'} then the maximum surrogate likelihood estimates are computed.
}
\details{
This program implemented one simulated likelihood method via sequential importance
         sampling (see Masarotto and Varin 2012), which is same as the method implemented in package
         \code{gcmr} (Masarotto and Varin 2016) except an antithetic variable is used. It also implemented
         one surrogate likelihood method via distributional transform (see Kazianka and Pilz 2010), which is
         generally faster.

         The argument \code{effort} is the sampling effort (known). It can be used to consider the heterogeneity
         of the measurement time or area at different locations.
         The default is 1 for all locations. See Han and De Oliveira (2016) for more details.
}
\examples{
\dontrun{
## Fit a Simulated Dataset with 100 locations
grid <- seq(0.05, 0.95, by = 0.1)
xloc <- expand.grid(x = grid, y = grid)[,1]
yloc <- expand.grid(x = grid, y = grid)[,2]

set.seed(123)
simData1 <- simgc(locs = cbind(xloc,yloc), sim.n = 1,
                    marginal = negbin.gc(mu = exp(1+xloc), od = 1),
                    corr = matern.gc(range = 0.4, kappa = 0.5, nugget = 0))

simFit1 <- mlegc(y = simData1$data, x = xloc, locs = cbind(xloc,yloc),
                 marginal = negbin.gc(link = 'log'),
                 corr = matern.gc(kappa = 0.5, nugget = FALSE), method = 'GHK')

simFit2 <- mlegc(y = simData1$data, x = xloc, locs = cbind(xloc,yloc),
                 marginal = negbin.gc(link = 'log'),
                 corr = matern.gc(kappa = 0.5, nugget = FALSE), method = 'GQT')
#summary(simFit1);summary(simFit2)
#plot(simFit1);plot(simFit2)



## Time consuming examples
## Fit a real dataset with 70 sampling locations.
data(Weed95)
weedobs <- Weed95[Weed95$dummy==1, ]
weedpred <- Weed95[Weed95$dummy==0, ]
Weedfit1 <- mlegc(y = weedobs$weedcount, x = weedobs[,4:5], locs = weedobs[,1:2],
                     marginal = poisson.gc(link='log'),
                     corr = matern.gc(kappa = 0.5, nugget = TRUE),
                     method = 'GHK')
summary(Weedfit1)
plot(Weedfit1)


## Fit a real dataset with 256 locations
data(LansingTrees)
Treefit1 <- mlegc(y = LansingTrees[,3], x = LansingTrees[,4], locs = LansingTrees[,1:2],
                  marginal = negbin.gc(link = 'log'),
                  corr = matern.gc(kappa = 0.5, nugget = FALSE), method = 'GHK')
summary(Treefit1)
plot(Treefit1)

# Try to use GQT method
Treefit2<- mlegc(y = LansingTrees[,3], x = LansingTrees[,4],
                locs = LansingTrees[,1:2], marginal = poisson.gc(link='log'),
                corr = matern.gc(kappa = 0.5, nugget = TRUE), method = 'GQT')
summary(Treefit2)
plot(Treefit2)

## Fit a real dataset with randomized locations
data(AtlanticFish)
Fitfish <- mlegc(y = AtlanticFish[,3], x = AtlanticFish[,4:6], locs = AtlanticFish[,1:2],
                   longlat = TRUE, marginal = negbin.gc(link='log'),
                   corr = matern.gc(kappa = 0.5, nugget = TRUE), method = 'GHK')
summary(Fitfish)

## Fit a real dataset with binomial counts; see Masarotto and Varin (2016).
library(gcmr)
data(malaria)
malariax <- data.frame(netuse = malaria$netuse,
                       green = malaria$green/100,
                       phc = malaria$phc)
Fitmalaria <- mlegc(y = malaria$cases, x = malariax, locs = malaria[,1:2],
                    marginal = binomial.gc(link='logit'), corrpar0 = 1.5,
                    corr = matern.gc(kappa = 0.5, nugget = FALSE),
                    distscale = 0.001, effort = malaria$size, method = 'GHK')
summary(Fitmalaria)


## Fit a real spatial binary dataset with 333 locations using probit link
data(OilWell)
Oilest1 <- mlegc(y = OilWell[,3], x = NULL, locs = OilWell[,1:2],
                 marginal = binomial.gc(link = 'probit'),
                 corr = matern.gc(nugget = TRUE), method = 'GHK')
summary(Oilest1)
plot(Oilest1, col = 2)
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

 Masarotto, G. and Varin, C. (2012)
 Gaussian copula marginal regression.
 \emph{Electronic Journal of Statistics 6:1517-1549}.
 \url{https://projecteuclid.org/euclid.ejs/1346421603}.

Masarotto, G. and Varin C. (2017). Gaussian Copula Regression in R. \emph{Journal of Statistical Software}, \bold{77}(8), 1--26. \doi{10.18637/jss.v077.i08}.

Han, Z. and De Oliveira, V. (2018) gcKrig: An R Package for the Analysis of Geostatistical Count Data Using Gaussian Copulas.
\emph{Journal of Statistical Software}, \bold{87}(13), 1--32. \doi{10.18637/jss.v087.i13}.

}
\seealso{
\code{\link[gcmr]{gcmr}}
}

\keyword{Estimation}
