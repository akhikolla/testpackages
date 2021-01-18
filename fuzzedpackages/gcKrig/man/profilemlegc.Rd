\name{profile.mlegc}
\alias{profile.mlegc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Profile Likelihood Based Confidence Interval of Parameters
for Gaussian Copula Models in Geostatistical Count Data
}

\usage{
\method{profile}{mlegc}(fitted, par.index, alpha = 0.05, start.point = NULL,
        method = 'GQT', nrep = 1000, seed = 12345, ...)
          }

\description{
This function computes the (approximate) profile likelihood based confidence interval. The algorithm
starts by choosing two starting points at different sides of the MLE and using an iterative process
to find the approximate lower and upper bound.
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitted}{an object of class \code{mlegc}, typically inherited from function \code{\link{mlegc}}.
}

 \item{par.index}{the index of the parameter which should be profiled.  }

  \item{alpha}{the significance level, default is \code{0.05} which corresponds to \code{95}
  percent confidence interval.}

  \item{start.point}{numeric vector of length 2 indicating the starting points for finding the
  left and right bound. If \code{start.point = NULL} then the default starting points will be used.}

\item{method}{Two methods are implemented. If
\code{method = 'GHK'} then the simulated likelihood will be used, if
\code{method = 'GQT'} then the surrogate likelihood will be used.}

\item{nrep}{the Monte Carlo size of the importance sampling algorithm for likelihood approximation;
only need to be specified if \code{method = 'GHK'}.}

\item{seed}{seed of the pseudorandom generator used in Monte Carlo simulation;
only need to be specified if \code{method = 'GHK'}.}

\item{...}{other arguments passed.}

}




\value{
Lower and upper bounds of the approximate confidence interval.
}
\references{
 Masarotto, G. and Varin, C. (2012)
 Gaussian copula marginal regression.
 \emph{Electronic Journal of Statistics 6:1517-1549}.
 \url{https://projecteuclid.org/euclid.ejs/1346421603}.

 Masarotto, G. and Varin C. (2017). Gaussian Copula Regression in R. \emph{Journal of Statistical Software}, \bold{77}(8), 1--26. \doi{10.18637/jss.v077.i08}.

Han, Z. and De Oliveira, V. (2018) gcKrig: An R Package for the Analysis of Geostatistical Count Data Using Gaussian Copulas.
\emph{Journal of Statistical Software}, \bold{87}(13), 1--32. \doi{10.18637/jss.v087.i13}.
}
\author{
Zifei Han \email{hanzifei1@gmail.com}
}


\seealso{
\code{\link{mlegc}}
}
\examples{
\dontrun{
data(LansingTrees)
Treefit4 <- mlegc(y = LansingTrees[,3], x = LansingTrees[,4],
                    locs = LansingTrees[,1:2], marginal = zip.gc(link = 'log'),
                    corr = matern.gc(kappa = 0.5, nugget = TRUE), method = 'GHK')
summary(Treefit4)

profile(Treefit4, 1, 0.05,  method = 'GHK', nrep = 1000, seed = 12345)
profile(Treefit4, 2, 0.05,  method = 'GHK', nrep = 1000, seed = 12345)
profile(Treefit4, 3, 0.05,  method = 'GHK', nrep = 1000, seed = 12345)
profile(Treefit4, 4, 0.05,  method = 'GHK', nrep = 1000, seed = 12345)
profile(Treefit4, 5, 0.05, method = 'GHK', nrep = 1000, seed = 12345)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{Estimation}
