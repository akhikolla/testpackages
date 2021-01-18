\name{k_tJADE}
\alias{k_tJADE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
k-tJADE for Tensor-Valued Observations 
}
\description{
Computes the faster ``k''-version of tensorial JADE in an independent component model. 
}
\usage{
k_tJADE(x, k = NULL, maxiter = 100, eps = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numeric array of an order at least two. It is assumed that the last dimension corresponds to the sampling units.}
  \item{k}{A vector with one less element than dimensions in \code{x}. The elements of \code{k} give upper bounds for cumulant matrix indices we diagonalize in each mode. Lower values mean faster computation times. The default value \code{NULL} puts k equal to 1 in each mode (the fastest choice).}
  \item{maxiter}{Maximum number of iterations. Passed on to \code{\link[JADE]{rjd}}.}
  \item{eps}{Convergence tolerance. Passed on to \code{\link[JADE]{rjd}}.}
}

\details{
It is assumed that \eqn{S} is a tensor (array) of size \eqn{p_1 \times p_2 \times \ldots \times p_r}{p_1 x p_2 x ... x p_r} with mutually independent elements and measured on \eqn{N} units. The tensor independent component model further assumes that the tensors S are mixed from each mode
\eqn{m} by the mixing matrix \eqn{A_m}, \eqn{m = 1, \ldots, r}{m= 1, ..., r}, yielding the observed data \eqn{X}. In R the sample of \eqn{X} is saved as an \code{\link{array}} of dimensions
\eqn{p_1, p_2, \ldots, p_r, N}{p_1, p_2, ..., p_r, N}.

\code{k_tJADE} recovers then based on \code{x} the underlying independent components \eqn{S} by estimating the \eqn{r} unmixing matrices 
 \eqn{W_1, \ldots, W_r}{W_1, ..., W_r} using fourth joint moments at the same time in a more efficient way than \code{\link{tFOBI}} but also in fewer numbers than \code{\link{tJADE}}. \code{k_tJADE} diagonalizes in each mode only those cumulant matrices \eqn{C^{ij}}{C^ij} for which \eqn{|i - j| < k_m}{|i - j| < k_m}.

If \code{x} is a matrix, that is, \eqn{r = 1}{r = 1}, the method reduces to JADE and the function calls \code{\link[JADE]{k_JADE}}.

}
\value{
A list with class 'tbss', inheriting from class 'bss', containing the following components: 
\item{S }{Array of the same size as x containing the independent components.}
\item{W }{List containing all the unmixing matrices}
\item{Xmu }{The data location.}
\item{k}{The used vector of \eqn{k}{k}-values.}
\item{datatype}{Character string with value "iid". Relevant for \code{\link{plot.tbss}}.}
}

\references{
\cite{Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2013), Fast Equivariant JADE, 
In the Proceedings of \emph{38th IEEE International Conference on Acoustics, Speech, and Signal Processing (ICASSP 2013)}, 6153--6157, \doi{10.1109/ICASSP.2013.6638847}}

\cite{Virta J., Li B., Nordhausen K., Oja H. (2018): JADE for tensor-valued observations, Journal of Computational and Graphical Statistics, Volume 27, p. 628 - 637, \doi{10.1080/10618600.2017.1407324}} 

\cite{Virta J., Lietzen N., Ilmonen P., Nordhausen K., (2018), Asymptotically and computationally efficient tensorial JADE, preprint available on arXiv {https://arxiv.org/abs/1808.00791}.}
}

\author{
Joni Virta
}

\seealso{
\code{\link[JADE]{k_JADE}}, \code{\link{tJADE}}, \code{\link[JADE]{JADE}}
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


k_tjade <- k_tJADE(X)

MD(k_tjade$W[[1]], A1)
MD(k_tjade$W[[2]], A2) 
tMD(k_tjade$W, list(A1, A2))


k_tjade <- k_tJADE(X, k = c(2, 1))

MD(k_tjade$W[[1]], A1)
MD(k_tjade$W[[2]], A2) 
tMD(k_tjade$W, list(A1, A2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{array }

