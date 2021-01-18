\name{intervalSystem}
\alias{intervalSystem}
\alias{intervalsystem}

\title{Interval systems}

\description{
Overview about the supported interval systems. More details are given in section 6 of the vignette.
}

\details{
The following interval systems (set of intervals on which tests will be performed) are available. Intervals are given as indices of observations / sample points.

\describe{
\item{\code{"all"}}{all intervals. More precisely, the set \eqn{\{[i, j], 1 \leq i \leq j \leq n\}}{{[i, j], 1 \le i \le j \le n}}. This system allows all lengths \code{1:n}.}
\item{\code{"dyaLen"}}{all intervals of dyadic length. More precisely, the set \eqn{\{[i, j], 1 \leq i \leq j \leq n\ s.t.\ j - i + 1 = 2^k,\ k\in N_0\}}{{[i, j], 1 \le i \le j \le n s.t. j - i + 1 = 2^k, k a natural number}}. This system allows all lengths of dyadic length \code{2^(0:as.integer(floor(log2(n)) + 1e-6))}.}
\item{\code{"dyaPar"}}{the dyadic partition, i.e. all disjoint intervals of dyadic length. More precisely, the set \eqn{\{[(i - 1) * 2^k + 1, i * 2^k], i = 1,\ldots, \lfloor n / 2^k\rfloor,\ k = 0, \ldots, \lfloor\log_2(n)\rfloor\}}{{[(i - 1) * 2^k + 1, i * 2^k], i = 1,\ldots, floor(n / 2^k), k = 0, \ldots, floor(log2(n))}}. This system allows all lengths of dyadic length \code{2^(0:as.integer(floor(log2(n)) + 1e-6))}.}
}

The interval system is selected via the \code{intervalSystem} argument, providing the corresponding string. By default (\code{NULL}) the default interval system of the specified parametric family will be used, which one this will be is described in \code{\link{parametricFamily}}. With the additional argument \code{lengths} it is possible to specify a set of lengths such that only tests on intervals with a length contained in this set will be performed. The set of lengths has to be a subset of all lengths that are allowed by the interval system and the parametric family. By default (\code{NULL}) all lengths allowed by the interval system and the parametric family are used.
}

\seealso{\link{parametricFamily}}

\examples{
y <- c(rnorm(50), rnorm(50, 2))

# interval system of all intervals and all lengths
fit <- stepFit(y, alpha = 0.5, intervalSystem = "all", lengths = 1:100,
               jumpint = TRUE, confband = TRUE)

# default for family "gauss" if number of observations is 1000 or less
identical(stepFit(y, alpha = 0.5, jumpint = TRUE, confband = TRUE), fit)

# intervalSystem "dyaLen" and a subset of lengths
\donttest{!identical(stepFit(y, alpha = 0.5, intervalSystem = "dyaLen", lengths = c(2, 4, 16),
                   jumpint = TRUE, confband = TRUE), fit)

# default for lengths are all possible lengths of the interval system
# and the parametric family
identical(stepFit(y, alpha = 0.5, intervalSystem = "dyaPar",
                  jumpint = TRUE, confband = TRUE),
          stepFit(y, alpha = 0.5, intervalSystem = "dyaPar", lengths = 2^(0:6),
                  jumpint = TRUE, confband = TRUE))

# interval system "dyaPar" is default for parametric family "hsmuce"
# length 1 is not possible for this parametric family
identical(stepFit(y, alpha = 0.5, family = "hsmuce",
                  jumpint = TRUE, confband = TRUE),
          stepFit(y, alpha = 0.5, family = "hsmuce", intervalSystem = "dyaPar",
                  lengths = 2^(1:6), jumpint = TRUE, confband = TRUE))

# interval system "dyaLen" is default for parametric family "mDependentPS"
identical(stepFit(y, alpha = 0.5, family = "mDependentPS", covariances = c(1, 0.5),
                  jumpint = TRUE, confband = TRUE),
          stepFit(y, alpha = 0.5, family = "mDependentPS", covariances = c(1, 0.5),
                  intervalSystem = "dyaLen", lengths = 2^(0:6),
                  jumpint = TRUE, confband = TRUE))
}}
