\name{corr.gc}
\alias{corr.gc}
\title{
Spatial Correlation Functions for Simulation, Likelihood Inference and Spatial Prediction in Gaussian Copula Models with
Geostatistical Count Data
}
\description{
Class of isotropic correlation functions available in the \code{gcKrig} library.
}

\details{
By default, range parameter is not provided, so
this function is used for likelihood inference and spatial prediction with function \code{mlegc}
and \code{predgc}.
Users need to specify if the correlation model includes a nugget effect
\code{nugget = TRUE} or not \code{nugget = FALSE}. For Matern and powered exponential correlation functions, the
shape parameter \code{kappa} is also required from users.

When both range and nugget parameters are given,
the function is used to specify the correlation structure in simulation
with function \code{simgc} in package \code{gcKrig}.

}


\value{
At the moment, the following three correlation functins are implemented:
\tabular{ll}{

\code{\link{matern.gc}} \tab the Matern correlation function. \cr
\code{\link{powerexp.gc}} \tab the powered exponential correlation function. \cr
\code{\link{spherical.gc}} \tab the spherical correlation function. \cr
}
}


\author{
Zifei Han \email{hanzifei1@gmail.com}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\references{
De Oliveira, V. (2013) Hierarchical Poisson models for spatial count data.
             \emph{Journal of Multivariate Analysis,122:393-408}.

Han, Z. and De Oliveira, V. (2018) gcKrig: An R Package for the Analysis of Geostatistical Count Data Using Gaussian Copulas.
\emph{Journal of Statistical Software}, \bold{87}(13), 1--32. \doi{10.18637/jss.v087.i13}.

}

\seealso{
         \code{\link{matern.gc}},
         \code{\link{powerexp.gc}},
         \code{\link{spherical.gc}}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Correlation}% use one of  RShowDoc("KEYWORDS")
\keyword{Simulation}% __ONLY ONE__ keyword per line
