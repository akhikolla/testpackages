\name{tFOBI}
\alias{tFOBI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{FOBI for Tensor-Valued Observations 
}
\description{
Computes the tensorial FOBI in an independent component model. 
}
\usage{
tFOBI(x, norm = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numeric array of an order at least two. It is assumed that the last dimension corresponds to the sampling units.}
  \item{norm}{A Boolean vector with number of entries equal to the number of modes in a single observation. The elements tell which modes use the ``normed'' version of tensorial FOBI. If \code{NULL} then all modes use the non-normed version.}
}
\details{
It is assumed that \eqn{S} is a tensor (array) of size \eqn{p_1 \times p_2 \times \ldots \times p_r}{p_1 x p_2 x ... x p_r} with mutually independent elements and measured on \eqn{N} units. The tensor independent component model further assumes that the tensors S are mixed from each mode
\eqn{m} by the mixing matrix \eqn{A_m}, \eqn{m = 1, \ldots, r}{m= 1, ..., r}, yielding the observed data \eqn{X}. In R the sample of \eqn{X} is saved as an \code{\link{array}} of dimensions
\eqn{p_1, p_2, \ldots, p_r, N}{p_1, p_2, ..., p_r, N}.

\code{tFOBI} recovers then based on \code{x} the underlying independent components \eqn{S} by estimating the \eqn{r} unmixing matrices 
\eqn{W_1, \ldots, W_r}{W_1, ..., W_r} using fourth joint moments. 
 
The unmixing can in each mode be done in two ways, using a ``non-normed'' or ``normed'' method and this is controlled by the argument \code{norm}. The authors advocate the general use of non-normed version, see the reference below for their comparison.
 
If \code{x} is a matrix, that is, \eqn{r = 1}{r = 1}, the method reduces to FOBI and the function calls \code{\link[JADE]{FOBI}}.
 
For a generalization for tensor-valued time series see \code{\link{tgFOBI}}.
}
\value{
A list with class 'tbss', inheriting from class 'bss', containing the following components: 
\item{S }{Array of the same size as x containing the independent components.}
\item{W }{List containing all the unmixing matrices.}
\item{norm}{The vector indicating which modes used the ``normed'' version.}
\item{Xmu }{The data location.}
\item{datatype}{Character string with value "iid". Relevant for \code{\link{plot.tbss}}.}
}
\references{
\cite{Virta, J., Li, B., Nordhausen, K. and Oja, H., (2017), Independent component analysis for tensor-valued data, Journal of Multivariate Analysis, \doi{10.1016/j.jmva.2017.09.008}}
}
\author{
Joni Virta
}



\seealso{
\code{\link[JADE]{FOBI}}, \code{\link{tgFOBI}}
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

MD(tfobi$W[[1]], A1)
MD(tfobi$W[[2]], A2) 
tMD(tfobi$W, list(A1, A2))

# Digit data example

data(zip.train)
x <- zip.train

rows <- which(x[, 1] == 0 | x[, 1] == 1)
x0 <- x[rows, 2:257]
y0 <- x[rows, 1] + 1

x0 <- t(x0)
dim(x0) <- c(16, 16, 2199)

tfobi <- tFOBI(x0)
plot(tfobi, col=y0)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{array }
