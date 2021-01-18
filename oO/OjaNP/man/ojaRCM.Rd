\name{ojaRCM}
\alias{ojaRCM}
\title{ Oja Rank Convariance Matrix }
\description{
  The function computes the Oja rank covariance matrix of a data set \code{X}.
}
\usage{
ojaRCM(X, p = NULL, silent = FALSE, na.action = na.fail)
}

\arguments{
\item{X}{numeric data.frame or matrix containing the data points as rows.}

\item{p}{\code{NULL} or a number between 0 and 1 which specifies the fraction of hyperplanes to be used for subsampling. If \code{p = 1}, no subsampling is done.
If \code{p = NULL}, the value of \code{p} is determined based on the size of the data set. See function \code{\link{ojaRank}} for details.}
  
\item{silent}{logical, if subsampling is done or the expected computation time is too long, a warning message will be printed unless \code{silent} is \code{TRUE}. The default is \code{FALSE}.}

\item{na.action}{a function which indicates what should happen when the data
          contain 'NA's. Default is to fail.}
}



\details{
The function computes the Oja rank covariance matrix of the data set \code{X}, that is (since Oja ranks are centered) the covariance matrix of the Oja ranks of the data points 
in \code{X}, taken w.r.t. the data set \code{X}. 

For a definition of the \dfn{Oja rank covariance matrix} and its properties see references below.
The matrix \code{X} needs to have at least as many rows as columns in order to give sensible results. 
The return value is a quadratic, symmetric matrix having as many columns as \code{X}. It works also for matrices \code{X} with only one column and also vectors, 
but note that the variance of univariate ranks does not yield much information about the data.

The function offers a subsampling option in order to speed up computation for large data sets. The subsampling fraction is controlled by the parameter \code{p}. 
If \code{p} is not specified (which defaults to \code{p = NULL}), it 
is automatically determined based on the dimension of the problem. The function tries to realize a reasonable compromise between accuracy and computing time, 
that is, for sufficiently small data matrices \code{X} the sampling fraction \code{p} is set to 1. 
Subsampling is applied to hyperplanes, not data points. A sample is drawn once, all Oja ranks are then computed based on this sample. 
For further details on subsampling see function \code{\link{ojaRank}}. Subsampling is useful. Even for very small \code{p} useable results can be expected, see e.g. Example 2.

}


\value{
  a symmetric matrix with \code{ncol(X)} columns and rows.
}
\references{ 
\cite{Fischer D, Mosler K, \enc{Möttönen}{Mottonen} J, Nordhausen K, Pokotylo O and Vogel D (2020). “Computing the Oja Median in R: The Package OjaNP.” \emph{Journal of Statistical Software}, \bold{92}(8), pp. 1-36. doi: 10.18637/jss.v092.i08 (URL:
http://doi.org/10.18637/jss.v092.i08).}

\cite{Visuri, S., Koivunen, V., Oja, H. (1999), Sign and rank covariance matrices, 
\emph{J. Stat. Plann. Inference}, \bold{91}, 557--575.} 

\cite{Ollila, E., Croux, C., Oja, H. (2004), Influence function and asymptotic effiency of the affine equivariant rank covariance matrix, \emph{Statistica Sinica}, \bold{14}, 297--316. }
}

\author{Daniel Vogel}

\seealso{ 
\code{\link{ojaRank}}, \code{\link{ojaSCM}} 
}

\examples{
### ----<< Example 1 >>---- : biochem data
data(biochem)
X <- biochem[,1:2]
ojaRCM(X)

# Oja ranks are centered 
# (i.e. they add up to zero), and 
# the following two return the same.
ojaRCM(X)
(1 - 1/nrow(X))*cov(ojaRank(X))



### ----<< Example 2 >>---- : 300 points in R^7 
# The merit of subsampling.
# The following example might take a bit longer:
\dontrun{
A <- matrix(c(1,0.5,1,4,2,0.5,-0.5,1,4), ncol = 3)
B <- A \%x\% A;  Sigma  <- (B \%*\% t(B))[1:7, 1:7]
# Sigma is some arbitrary positive definite matrix.
set.seed(123)
X <- rmvnorm(n = 300, sigma = Sigma) 

cov2cor(Sigma) # the true correlation matrix
cor(X)  # Bravais-Pearson correlation
cov2cor(solve(ojaRCM(X))
# correlation estimate based on Oja ranks 
# The subsampling fraction in this example
# is p = 1.081438e-10.
# Yet it returns a sensible estimate.
}
}
\keyword{multivariate}
\keyword{nonparametric}
