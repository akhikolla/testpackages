\name{marginal.gc}
\alias{marginal.gc}
\title{
Marginals for Data Simulation,
Correlation Assessment,
Likelihood Inference and Spatial Prediction in Gaussian Copula Models
for Geostatistical Data
}
\description{
Class of marginals available in \code{gcKrig} library for geostatistical data simulation,
correlation structure assessment (both continuous and discrete marginals)
and model inferences (discrete marginals only).
In former cases parameters of
the marginals are given by users, otherwise
parameters are estimated from the data (except when doing prediction
with function \code{predgc},
users can choose to either input known estimates or
estimate the parameters with input data).
}

\details{
The package \code{gcKrig} does not include inference and prediction functionalities
for continuous marginals. For inference with continuous marginals, see Masarotto and Varin (2012).

By default, when the marginals are discrete, they are used for
estimation with function \code{mlegc} and prediction with
function \code{predgc}.
They can be used in function \code{simgc} and \code{corrTG} as well for
the purpose of data simulation and correlation computation in a transformed
Gaussian random field (Han and De Oliveira, 2016), if parameter values are specified.

For continuous marginals, they are used for simulation with function
\code{simgc} and correlation computation with \code{corrTG} only, so
parameters should always be specified.

}
\value{
At the moment, the following marginals are implemented:
\tabular{ll}{

\code{\link{beta.gc}} \tab beta marginals. \cr
\code{\link{binomial.gc}} \tab binomial marginals. \cr
\code{\link{gm.gc}} \tab gamma marginals. \cr
\code{\link{gaussian.gc}} \tab Gaussian marginals. \cr
\code{\link{negbin.gc}} \tab negative binomial marginals. \cr
\code{\link{poisson.gc}} \tab Poisson marginals. \cr
\code{\link{weibull.gc}} \tab Weibull marginals. \cr
\code{\link{zip.gc}} \tab zero-inflated Poisson marginals. \cr
}
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}

\references{
Cribari-Neto, F. and Zeileis, A. (2010) Beta regression in R. \emph{Journal of Statistical Software, 34:1-24}. \url{http://www.jstatsoft.org/v34/i02/}.

Ferrari, S.L.P. and Cribari-Neto, F. (2004) Beta regression for modeling rates and proportions. \emph{Journal of Applied Statistics, 31:799-815}.

Han, Z. and De Oliveira, V. (2016) On the correlation structure of Gaussian
            copula models for geostatistical count data.
            \emph{Australian and New Zealand Journal of Statistics, 58:47-69}.

Masarotto, G. and Varin, C. (2012) Gaussian copula marginal regression.
\emph{Electronic Journal of Statistics, 6:1517-1549}.
\url{https://projecteuclid.org/euclid.ejs/1346421603}.

Masarotto, G. and Varin C. (2017) Gaussian Copula Regression in R.
\emph{Journal of Statistical Software}, \bold{77}(8), 1--26. \doi{10.18637/jss.v077.i08}.

Han, Z. and De Oliveira, V. (2018) gcKrig: An R Package for the Analysis of Geostatistical Count Data Using Gaussian Copulas.
\emph{Journal of Statistical Software}, \bold{87}(13), 1--32. \doi{10.18637/jss.v087.i13}.

}


\seealso{\code{\link{beta.gc}},
         \code{\link{binomial.gc}},
         \code{\link{gm.gc}},
         \code{\link{gaussian.gc}},
         \code{\link{negbin.gc}},
         \code{\link{poisson.gc}},
         \code{\link{weibull.gc}},
         \code{\link{zip.gc}}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Marginal}
\keyword{Simulation}% use one of  RShowDoc("KEYWORDS")
\keyword{Estimation}
\keyword{Prediction}
