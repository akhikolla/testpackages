\name{penalty}
\alias{penalty}
\alias{penalties}

\title{Penalties}

\description{
Overview about the supported penalties. More details are also given in section 3.2 of the vignette.
}

\details{
The penalties (ways to balance different scales) can be divided into two groups: scale penalisation and balancing by weights. More precisely, the scale penalisations \code{"sqrt"}, \code{"log"} and \code{"none"} and balancing by weights called \code{"weights"} are available.

Let \code{T} be the unpenalised test statistic of the specified parametric family on an interval of length \code{l} and \code{nq} the number of observations used for the penalisation, typically the number of observations \code{n} but can also be chosen larger.

\describe{
\item{\code{"sqrt"}}{penalised statistic is \code{sqrt(2 * T) - sqrt(2 * log(exp(1) * nq / l)}. This penalisation is proposed in (Frick et al., 2014) and guarantees for most parametric families that the penalised multiscale statistic is asymptotically finite. This is not true for parametric family \code{"hsmuce"}. Hence, this penalisation is recommended and the default one for the parametric families \code{"gauss"} and \code{"mDependentPS"}, but not for \code{"hsmuce"}.}
\item{\code{"log"}}{penalised statistic is \code{T - log(exp(1) * nq / l)}. This penalisation is outdated and only still supported for comparisons.}
\item{\code{"none"}}{no penalisation, penalised statistic is equal to the unpenalised. Multiscale regression without a penalisation is not recommend.}
\item{\code{"weights"}}{critical values will be computed by weights, see section 3.2.2 in the vignette and (Pein et al., 2017, section 2) for more details. This penalty is recommend and the default one for the parametric family \code{"hsmuce"}, but can also be used for other families. Will be replaced by \code{"none"} in \code{\link{computeStat}} and \code{\link{monteCarloSimulation}}.}
}

The penalisation is selected via the \code{penalty} argument providing the corresponding string. If \code{NULL} the default penalty of the specified parametric family will be used, see \code{\link{parametricFamily}} for which one this will be.
}

\seealso{\link{parametricFamily}, \code{\link{critVal}}}

\references{
Frick, K., Munk, A., Sieling, H. (2014) Multiscale change-point inference. With discussion and rejoinder by the authors. \emph{Journal of the Royal Statistical Society, Series B} \bold{76}(3), 495--580.

Pein, F., Sieling, H., Munk, A. (2017) Heterogeneous change point inference. \emph{Journal of the Royal Statistical Society, Series B}, \bold{79}(4), 1207--1227.
} 

\examples{
set.seed(1)
y <- c(rnorm(50), rnorm(50, 2))

# penalty "sqrt"
fit <- stepFit(y, alpha = 0.5, penalty = "sqrt", jumpint = TRUE, confband = TRUE)

# default for family "gauss"
identical(stepFit(y, alpha = 0.5, jumpint = TRUE, confband = TRUE), fit)

# penalty "weights"
!identical(stepFit(y, alpha = 0.5, penalty = "weights",
                   jumpint = TRUE, confband = TRUE), fit)

# penalty "weights" is default for parametric family "hsmuce"
# by default equal weights are chosen
identical(stepFit(y, alpha = 0.5, family = "hsmuce",
                  jumpint = TRUE, confband = TRUE),
          stepFit(y, alpha = 0.5, family = "hsmuce", penalty = "weights",
                  weights = rep(1 / 6, 6), jumpint = TRUE, confband = TRUE))

# different weights
!identical(stepFit(y, alpha = 0.5, family = "hsmuce", weights = 6:1 / sum(6:1),
                   jumpint = TRUE, confband = TRUE),
           stepFit(y, alpha = 0.5, family = "hsmuce", penalty = "weights",
                   weights = rep(1 / 6, 6), jumpint = TRUE, confband = TRUE))

# penalty "sqrt is default for parametric family "mDependentPS"
identical(stepFit(y, alpha = 0.5, family = "mDependentPS", covariances = c(1, 0.5),
                  jumpint = TRUE, confband = TRUE),
          stepFit(y, alpha = 0.5, family = "mDependentPS", covariances = c(1, 0.5),
                  penalty = "sqrt", jumpint = TRUE, confband = TRUE))

}
