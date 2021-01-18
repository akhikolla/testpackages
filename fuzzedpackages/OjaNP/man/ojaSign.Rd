\name{ojaSign}
\alias{ojaSign}
\title{ Oja Signs -- Affine Equivariant Multivariate Signs }
\description{
  The function computes the Oja sign of a point \code{x} w.r.t. a data set \code{X} or, if no
  point \code{x} is given, the Oja signs of all points in \code{X}.
}
\usage{
ojaSign(X, x = NULL, center = "ojaMedian", p = NULL, silent = FALSE, 
        na.action = na.fail, ...)
}

\arguments{
\item{X}{numeric data.frame or matrix containing the data points as rows.}

\item{x}{\code{NULL} or a numeric vector, the point for which the Oja sign should be computed.}
  
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
	\item	\code{"spatialMedian"} (function \code{\link[ICSNP]{spatial.median}} from \cr package \code{ICSNP}),
	\item	\code{"compMedian"}  (marginal median) or
	\item	\code{"HRMedian"}  (Hettmansperger and Randles median,\cr function \code{\link[ICSNP]{HR.Mest}} from package \code{ICSNP}).			
	}% itemize		
}% itemize
The default is \code{"ojaMedian"}.
} % \item von \arguments

\item{p}{\code{NULL} or a number between 0 and 1 which specifies the fraction of hyperplanes to be used for subsampling. If \code{p = 1}, no subsampling is done.
If \code{p = NULL}, the value of \code{p} is determined based on the size of the data set. See details.
}

\item{silent}{logical, if subsampling is done or the expected computation time is too long, a warning message will be printed unless \code{silent} is \code{TRUE}. The default is \code{FALSE}.
}

\item{na.action}{a function which indicates what should happen when the data contain 'NA's. Default is to fail.}

\item{\dots}{arguments passed on to the location function.}
}% arguments



\details{
The function computes the Oja sign of the point \code{x} w.r.t. to the data set \code{X} or, if no \code{x} is specified, the Oja signs of all data points in \code{X} w.r.t. \code{X}.
For a definition of \dfn{Oja sign} see reference below.

The matrix \code{X} needs to have at least as many rows as columns in order to give sensible results. 
The vector \code{x} has to be of length \code{ncol(X)}. 
If \code{x} is specified, a vector of length \code{ncol(X)} is returned. Otherwise the return value is a matrix of the same dimensions as \code{X} where the 
\eqn{i}-th row contains the Oja sign of the \eqn{i}-th row of \code{X}.
The matrix \code{X} must have at least \bold{two columns}. For univariate signs use \code{\link[base]{sign}}.

Oja signs (contrary to Oja ranks) require the computation of a center (\dfn{location}) of the data cloud. The function offers 
various ways to do this. One can explicitly pass a location as a numeric vector (which has to be of length \code{ncol(X)}), 
or pass a function that computes a multivariate location. %Some care is needed when using this option.

This function must accept data sets of the same type as \code{X} (i.e. row-wise) and return a numeric vector of length \code{ncol(X)}. 
For example \code{center = colMeans} will do, whereas \code{center = mean} will not, since \code{mean(X)} computes the (univariate) mean of 
all elements of the matrix \code{X}. Note that some location functions return a list also containing information on the data or the computation 
rather than the numeric location estimate only. Via \code{...} arguments can be passed on to the location function, see example 2 below. 
The mean and several nonparametric location estimates are implemented and can also be specified by passing a corresponding string. 
The \code{...} option is then available, too. Being based on the same nonparametric concept the Oja median is the natural location for Oja 
signs and is hence the default. But since it is the solution of an optimization problem, it may -- depending on the optimization algorithm -- return different values when run on the same data.

For large data sets the function offers a subsampling option in order to deliver (approximate) results within reasonable time. For \eqn{n} data points in 
\eqn{R^k}{the k-dimenstional space} the computation of the Oja sign necessitates the evaluation of \eqn{N = choose(n,k-1)} hyperplanes in \eqn{R^k}. If 
\eqn{p < 1} is passed to the function, the computation is based on a random sample of only \eqn{p N} of all possible \eqn{N} hyperplanes. If \code{p} is not specified, 
it is automatically determined based on \eqn{n} and \eqn{k} to yield a sensible trade-off between accuracy and computing time. If 
\eqn{N k^3 < 6 \cdot 10^6}{choose(n,k-1)*(k-1)^3 < 6e+06}, the sample fraction \code{p} is set to 1 (no subsampling). Otherwise \code{p} is chosen such that the computation 
(of one sign) \emph{usually} takes around 20 seconds (on a 1.66 GHz CPU and 1 GB RAM). If all Oja signs of \code{X} are requested, 
a hyperplane sample is drawn once, all Oja signs are then computed based on this sample.

Finally, subsampling is feasible. Even for very small \code{p} useable results can be expected, see e.g. the examples for the function \code{\link{ojaSCM}}.

Claudia \enc{Köllmann}{Kollmann} is acknowledged for bug-fixing this function.
}


\value{
Either a numeric vector, the Oja sign of \code{x}, or 
 a matrix of the same dimensions as \code{X} containing the Oja signs of \code{X} as rows.
}
\references{ 
\cite{Fischer D, Mosler K, \enc{Möttönen}{Mottonen} J, Nordhausen K, Pokotylo O and Vogel D (2020). “Computing the Oja Median in R: The Package OjaNP.” \emph{Journal of Statistical Software}, \bold{92}(8), pp. 1-36. doi: 10.18637/jss.v092.i08 (URL:
http://doi.org/10.18637/jss.v092.i08).}

\cite{Oja, H. (1999), Affine invariant multivariate sign and rank tests and corresponding estimates: A review, 
\emph{Scand. J. Statist.}, \bold{26}, 319--343.}
}


\author{Daniel Vogel, Jyrki \enc{Möttönen}{Mottonen}}



\seealso{
\code{\link{ojaRank}},  
\code{\link{ojaSignedRank}},  
\code{\link{ojaMedian}}, 
\code{\link[ICSNP]{spatial.median}}, 
\code{\link[ICSNP]{HR.Mest}}, 
\code{\link{ojaSCM}} 
}


\examples{
### ----<< Example 1 >>---- : 30 points in R^2
set.seed(123)
X <- rmvnorm(n = 30, mean = c(0,0)) # from package 'mvtnorm'
y <- c(100,100)
om <- ojaMedian(X, alg = "exact")

ojaSign(X)
ojaSign(X,y)
# possible ways of specifying the mean as location:
ojaSign(X, center = "colMean")
ojaSign(X, center = colMeans)
ojaSign(X, center = colMeans(X))

# The following two return the same (only in different time),
ojaSign(X, center = colMeans)
t(apply(X, 1, function(y){ojaSign(X, y, center = colMeans)}))

# but the following not (due to different subsampling).
# 1)
set.seed(123)
ojaSign(X, center = colMeans, p = 0.9, silent = TRUE)
# 2)
set.seed(123)
t(apply(X, 1, function(y){ojaSign(X, y, c = colMeans,p = 0.9, s = TRUE)}))
# In 1) one subsample for all signs is drawn, whereas in 2)
# a different sample for each sign is drawn.

### ----<< Example 2 >>---- : Oja median
# The Oja sign of the Oja median is zero:
ojaSign(X, x = om, alg = "exact") 
# The default location function 'ojaMedian()' 
# is called with method "exact",
# which gives the same result as:
ojaSign(X, x = om, center = om) 
# But note: The following is likely to not return zero.  
ojaSign(X, x = ojaMedian(X))
# The default method of 'ojaMedian()' is "evo", 
# which is fast, but returns approximate results.



### ----<< Example 3 >>---- : 400 points in R^6
# Subsampling is done.
# The following example might take a bit longer:
\dontrun{
set.seed(123)
X <- rmvnorm(n = 400, mean = rep(0, 6))
os1 <- ojaSign(X, x = 1:6, c = colMeans)
# Note: the following command may take several minutes
os2 <- ojaSign(X, x = 1:6, p = 0.0000001, c = colMeans)
}
}

\keyword{multivariate}
\keyword{nonparametric}
