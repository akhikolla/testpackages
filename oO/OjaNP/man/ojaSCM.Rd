\name{ojaSCM}
\alias{ojaSCM}
\title{ Oja Sign Convariance Matrix }
\description{
  The function computes the Oja sign covariance matrix of a data set \code{X}.
}
\usage{
ojaSCM(X, center = "ojaMedian", p = NULL, silent = FALSE, 
       na.action = na.fail, ...)
}


\arguments{
\item{X}{numeric data.frame or matrix containing the data points as rows.}

\item{center}{one of the following three: 
\itemize{
	\item 
	a numeric vector giving the location of the data, 
	\item
		a function that computes a multivariate location
		(see details below) or
	\item
		one of the following strings: 
	\itemize{
	\item	\code{"colMean"} (vector of means, function \code{\link{colMeans}} is called), 
	\item	\code{"ojaMedian"} (function \code{\link{ojaMedian}}),
	\item	\code{"spatialMedian"} (function \code{\link[ICSNP]{spatial.median}} from package \cr \code{ICSNP}),
	\item	\code{"compMedian"}  (marginal median) or
	\item	\code{"HRMedian"}  (Hettmansperger and Randles median, function \cr \code{\link[ICSNP]{HR.Mest}} from package \code{ICSNP}).			
	}	% inneres itemize
}	% au"seres itemize
The default is \code{"ojaMedian"}.
} % \item von \arguments

\item{p}{\code{NULL} or a number between 0 and 1 which specifies the fraction of hyperplanes to be used for subsampling. If \code{p = 1}, no subsampling is done.
If \code{p = NULL}, the value of \code{p} is determined based on the size of the data set. See function \code{\link{ojaSign}} for details.
}

\item{silent}{logical, if subsampling is done or the expected computation time is too long, a warning message will be printed unless \code{silent} is \code{TRUE}. The default is \code{FALSE}.
}

\item{na.action}{a function which indicates what should happen when the data contain 'NA's. Default is to fail.}

\item{\dots}{arguments passed on to the location function.}}


\details{
The function computes the Oja sign covariance matrix of the data set \code{X}, that is (if the Oja signs are centered by the Oja median) the covariance matrix of the Oja signs of the 
data points in \code{X}, taken w.r.t.  \code{X}. 

For a definition of the \dfn{Oja sign covariance matrix} and its properties see references below.
The matrix \code{X} needs to have at least two columns and at least as many rows as columns in order to give sensible results. 
The return value is a quadratic, symmetric matrix having as many columns as \code{X}.

Oja signs (contrary to Oja ranks) require the computation of a centre (\dfn{location}) of the data cloud. The function offers various ways to specify the location. For details on 
location computation see function \code{\link{ojaSign}}.

The function offers a subsampling option in order to speed up computation for large data sets. The subsampling fraction is controlled by the parameter \code{p}. If \code{p} is not 
specified (which defaults to \code{p = NULL}), it is automatically determined based on the dimension of the problem. The function tries to realize a reasonable compromise 
between accuracy and computing time, that is, for sufficiently small data matrices \code{X} the sampling fraction \code{p} is set to 1. Subsampling is applied to hyperplanes, 
not data points. A sample is drawn once, all Oja signs are then computed based on this sample. For further details on subsampling see function \code{\link{ojaSign}}. 
Subsampling is useful. Even for very small \code{p} useable results can be expected, see e.g. Example 2.
  
}
\value{
  a symmetric matrix with \code{ncol(X)} columns and rows.
}  

\references{
\cite{Fischer D, Mosler K, \enc{Möttönen}{Mottonen} J, Nordhausen K, Pokotylo O and Vogel D (2020). “Computing the Oja Median in R: The Package OjaNP.” \emph{Journal of Statistical Software}, \bold{92}(8), pp. 1-36. doi: 10.18637/jss.v092.i08 (URL:
http://doi.org/10.18637/jss.v092.i08).}

\cite{Visuri, S., Koivunen, V., Oja, H. (1999), Sign and rank covariance matrices, 
\emph{J. Stat. Plann. Inference}, \bold{91}, 557--575.}

\cite{Ollila, E., Oja, H., Croux, C. (2003), The affine equivariant sign covariance matrix: Asymptotic behavior and efficiencies, \emph{J. Multivariate Analysis}, \bold{87}, 328--355.}
}

\author{Daniel Vogel}

\seealso{ 
\code{\link{ojaSign}}, 
 \code{\link{ojaRCM}}, 
\code{\link{ojaMedian}}, 
\code{\link[ICSNP]{spatial.median}}, 
\code{\link[ICSNP]{HR.Mest}} }



\examples{

### ----<< Example 1 >>---- : biochem data
data(biochem)
X <- biochem[,1:2]
ojaSCM(X)

# Oja signs are correctly centered 
# (i.e. they add up to zero) when 
# computed w.r.t. the Oja median
# Hence the following return the same,
ojaSCM(X, center = "ojaMedian", alg = "exact")
(1 - 1/nrow(X))*cov(ojaSign(X, alg = "exact"))
# but the following not.
ojaSCM(X, center = "colMean")
(1 - 1/nrow(X))*cov(ojaSign(X, center = "colMean"))



### ----<< Example 2 >>---- : 300 points in R^7 
# The merit of subsampling.
# The following example might take a bit longer:
\dontrun{
A <- matrix(c(1,0.5,1,4,2,0.5,-0.5,1,4), ncol = 3)
B <- A \%x\% A;  Sigma  <- (B \%*\% t(B))[1:7, 1:7]
# Sigma is some arbitrary positive definite matrix.
set.seed(123)
X <- rmvnorm(n=300,sigma=Sigma) 

cov2cor(Sigma) # the true correlation matrix
cor(X)  # Bravais-Pearson correlation
cov2cor(solve(ojaSCM(X, center = "colMean")))
# correlation estimate based on Oja signs 
# The subsampling fraction in this example
# is p = 4.542038e-09.
# Yet it returns a sensible estimate.
}
}

\keyword{multivariate}
\keyword{nonparametric}
