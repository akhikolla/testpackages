\name{tPP}
\alias{tPP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Projection pursuit for Tensor-Valued Observations}
\description{
  Applies mode-wise projection pursuit to tensorial data with respect to the chosen measure of interestingness. 
}
\usage{
  tPP(x, nl = "pow3", eps = 1e-6, maxiter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numeric array of an order at least three. It is assumed that the last dimension corresponds to the sampling units.}
  \item{nl}{The chosen measure of interestingness/objective function. Current choices include \code{pow3} (default) and \code{skew}, see the details below}
  \item{eps}{The convergence tolerance of the iterative algortihm.}
  \item{maxiter}{The maximum number of iterations.}
}
\details{
  The observed tensors (arrays) \eqn{X} of size \eqn{p_1 \times p_2 \times \ldots \times p_r}{p_1 x p_2 x ... x p_r} measured on \eqn{N} units are standardized from each mode and then projected mode-wise onto the directions that maximize the \eqn{L_2}{L_2}-norm of the vector of the values \eqn{E[G(u_k^T X X^T u_k)] - E[G(c^2)]}{E[G(u_k^T X X^T u_k)] - E[G(c^2)]}, where \eqn{G}{G} is the chosen objective function and \eqn{c^2}{c^2} obeys the chi-squared distribution with \eqn{q}{q} degress of freedom. Currently the function allows the choices \eqn{G(x) = x^2}{G(x) = x^2} (\code{pow3}) and \eqn{G(x) = x \sqrt x}{G(x) = x \sqrt{x}} (\code{skew}), which correspond roughly to the maximization of kurtosis and skewness, respectively. The algorithm is the multilinear extension of FastICA, where the names of the objective functions also come from.
}
\value{
A list with class 'tbss', inheriting from class 'bss', containing the following components: 
\item{S }{Array of the same size as x containing the estimated components.}
\item{W }{List containing all the unmixing matrices.}
\item{iter }{The numbers of iteration used per mode.}
\item{Xmu }{The data location.}
\item{datatype}{Character string with value "iid". Relevant for \code{\link{plot.tbss}}.}
}
\references{
  \cite{Nordhausen, K. and Virta, J. (2018), Tensorial projection pursuit, Manuscript in preparation.}
  
  \cite{Hyvarinen, A. (1999) Fast and robust fixed-point algorithms for independent component analysis, IEEE transactions on Neural Networks 10.3: 626-634.}
}
\author{
  Joni Virta
}

\seealso{
\code{\link[fICA]{fICA}}, \code{\link[ICtest]{NGPP}}
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

tpp <- tPP(X)

MD(tpp$W[[1]], A1)
MD(tpp$W[[2]], A2) 
tMD(tpp$W, list(A1, A2))

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{array }
