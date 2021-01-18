\name{F.type.test.statistic}
\alias{F.type.test.statistic}
\encoding{latin1}
\title{Computation of F-type test statistics}
\description{This function computes test statistics of the form
\eqn{T_w} and \eqn{T_{E, \mathsf{W}}}
as defined in \enc{Pötscher}{Poetscher} and Preinerstorfer (2016). The weights
for \eqn{T_w} and for \eqn{T_{E, \mathsf{W}}} are obtained from a kernel
function 
(Bartlett, Parzen, or Quadratic Spectral kernel, which provide nonnegative 
definite covariance 
estimators) and a bandwidth parameter. See also the description of the argument
\code{ker} below for further details concerning the weights. The class of 
test statistics of the form \eqn{T_w} or \eqn{T_{E, \mathsf{W}}}
includes F-type tests based on covariance estimators with data-independent
bandwidth parameters and without prewhitening as considered in, e.g., 
Newey and West (1987), Andrews (1991), Kiefer and Vogelsang (2002, 2005), cf.
also Preinerstorfer and \enc{Pötscher}{Poetscher} (2016).}

\usage{
F.type.test.statistic(y, R, r, X, bandwidth, ker, Eicker = FALSE, cores = 1)
}

\arguments{
\item{y}{Either an observation vector, or a matrix the columns of which are
observation vectors. The
number of rows of an observation vector must coincide with the number of rows
of the design matrix \code{X}.}

\item{R}{The restriction matrix. \code{F.type.test.statistic} computes a test 
statistic for the hypothesis \eqn{R \beta = r}. \code{R} needs to be of
full row rank, and needs to have the same number of columns as \code{X}.}

\item{r}{The restriction vector. \code{F.type.test.statistic} computes a test 
statistic for the hypothesis \eqn{R \beta = r}. \code{r} needs to be a
vector with the same number of coordinates as the number of rows of \code{R}.}
                
\item{X}{The design matrix. \code{X} needs to be of full column rank. 
The number of columns of \code{X} must be smaller than the number of rows of 
\code{X}.}

\item{bandwidth}{Bandwidth parameter used in the construction of the test 
statistic. A positive real number.}

\item{ker}{Kernel function used in the construction of the test statistic. 
\code{ker} can take one of the
values "Bartlett", "Parzen", or "Quadratic Spectral". The 
\code{\link{kweights}}
function is used to generate 
the weights. For 
the test statistic \eqn{T_{w}} this implies the weights used via 
\eqn{w(j, n) = ker(j/bandwidth)}. For  the test statistic 
\eqn{T_{E, \mathsf{W}}} this
implies the weights matrix \eqn{\mathsf{W}} via \eqn{\mathsf{W}_{ij} = 
ker((i-j)/bandwidth)}, cf. 
\enc{Pötscher}{Poetscher} and Preinerstorfer (2016) for definitions of
\eqn{T_{w}} and
\eqn{T_{E, \mathsf{W}}}.}

\item{Eicker}{Determines the test statistic computed. If \code{Eicker = TRUE},
then \eqn{T_{E, \mathsf{W}}} (with 
\eqn{\mathsf{W}_{ij}=ker((i-j)/bandwidth)}) is computed on the input observation
vector(s) \code{y}. If 
\code{Eicker = FALSE}, then \eqn{T_{w}} (with \eqn{w(j, n) = 
ker(j/bandwidth)}) is computed on the input observation
vector(s) \code{y} (cf. \enc{Pötscher}{Poetscher} and
Preinerstorfer 
(2016) for a precise definition of these test statistics). Default is 
\code{Eicker = FALSE}.}

\item{cores}{The number of CPU cores used in the (parallelized)
computation of the test statistics. Default is 1. This can be used to speed up
the computation in case \code{y} is a matrix with many columns.  Parallelized 
computation is enabled only if the compiler used to build \pkg{acrt} 
supports OpenMP.}
}

\details{For details concerning the test statistics please
see the relevant sections in \enc{Pötscher}{Poetscher} and Preinerstorfer (2016)
.
}

\value{The function returns a list consisting of:
\item{test.val}{Either a vector the entries of which correspond to the values 
of the
test statistic evaluated at each column of the input matrix \code{y}, or, 
if \code{y} is a vector, the test statistic evaluated at \code{y}.}
}

\examples{
n <- 100
y <- rnorm(n)
X <- cbind(rep(1, length = n), rnorm(n))
R <- matrix(c(0, 1), nrow = 1, ncol = 2)
r <- 0
bandwidth <- n/10
ker <- "Bartlett"
F.type.test.statistic(y, R, r, X, bandwidth, ker)
}

\references{
Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation consistent 
covariance matrix estimation. \emph{Econometrica}, \strong{59} 817-858.

Kiefer, N. M. and Vogelsang, T. J. (2002). Heteroskedasticity - autocorrelation 
robust standard errors using the Bartlett kernel without truncation. 
\emph{Econometrica}, \strong{70} 2093-2095.

Kiefer, N. M. and Vogelsang, T. J. (2005). A new asymptotic theory for 
heteroskedasticity - autocorrelation robust tests. \emph{Econometric Theory},
\strong{21} 1130-1164.

Newey, W. K. and West, K. D. (1987). A simple, positive semi-definite, 
heteroskedasticity and autocorrelation consistent covariance matrix. 
\emph{Econometrica}, \strong{55} 703-708.

\enc{Pötscher}{Poetscher}, B.M. and Preinerstorfer, D. (2016). Controlling the 
size of autocorrelation robust tests. \url{https://arxiv.org/abs/1612.06127/}

Preinerstorfer, D. and \enc{Pötscher}{Poetscher}, B. M. (2016). On size and 
power of heteroskedasticity and autocorrelation robust tests. 
\emph{Econometric Theory}, \strong{32} 261-358.
}


\seealso{
\code{\link{kweights}}.
}  
