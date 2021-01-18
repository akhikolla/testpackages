\name{tMD}
\alias{tMD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Minimum Distance Index of a Kronecker Product
}
\description{
A shortcut function for computing the minimum distance index of a tensorial ICA estimate on the Kronecker product ``scale'' (the vectorized space).
}
\usage{
tMD(W.hat, A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W.hat}{A list of \code{r} unmixing matrix estimates, W_1, W_2, ..., W_r.
}
  \item{A}{A list of \code{r} mixing matrices, A_1, A_2, ..., A_r.
}
}
\details{
The function computes the minimum distance index between \code{W.hat[[r]] \%x\%  ... \%x\% W.hat[[1]]} and \code{A[[r]] \%x\%  ... \%x\% A[[1]]}. The index is useful for comparing the performance of a tensor-valued ICA method to that of a method using first vectorization and then some vector-valued ICA method.
}
\value{
The value of the MD index of the Kronecker product.
}
\references{
\cite{Ilmonen, P., Nordhausen, K., Oja, H. and Ollila, E. (2010), A New Performance Index for ICA: Properties, Computation and Asymptotic Analysis. In Vigneron, V., Zarzoso, V., Moreau, E., Gribonval, R. and Vincent, E. (editors) Latent Variable Analysis and Signal Separation, 229-236, Springer.}

\cite{Virta, J., Li, B., Nordhausen, K. and Oja, H., (2017), Independent component analysis for tensor-valued data, Journal of Multivariate Analysis, \doi{10.1016/j.jmva.2017.09.008}}
}
\author{
Joni Virta
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[JADE]{MD}}
}
\examples{
n <- 1000
S <- t(cbind(rexp(n)-1,
             rnorm(n),
             runif(n, -sqrt(3), sqrt(3)),
             rt(n,5)*sqrt(0.6),
             (rchisq(n,1)-1)/sqrt(2),
             (rchisq(n,2)-2)/sqrt(4)))

dim(S) <- c(3, 2, n)

A1 <- matrix(rnorm(9), 3, 3)
A2 <- matrix(rnorm(4), 2, 2)

X <- tensorTransform(S, A1, 1)
X <- tensorTransform(X, A2, 2)

tfobi <- tFOBI(X)

MD(tfobi$W[[2]] \%x\% tfobi$W[[1]], A2 \%x\% A1)
tMD(list(tfobi$W[[2]]), list(A2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")
\keyword{ array }% __ONLY ONE__ keyword per line
