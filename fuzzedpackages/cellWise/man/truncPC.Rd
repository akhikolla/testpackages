\name{truncPC}
\alias{truncPC}

\title{
Classical Principal Components by truncated SVD.
}
\description{
 Similar usage to robustbase::classPC except for the new argument \code{ncomb} which is the desired number of components. Only this many PC's are computed in order to save computation time. Makes use of \code{propack.svd} of package \pkg{svd}.
}
\usage{
truncPC(X, ncomp = NULL, scale = FALSE, center = TRUE, 
                    signflip = TRUE, via.svd = NULL, scores = FALSE)
}

\arguments{
  \item{X}{
a numeric matrix.
}
  \item{ncomp}{
the desired number of components (if not specified, all components are computed).
}
\item{scale}{
logical, or numeric vector for scaling the columns.
}
  \item{center}{
logical or numeric vector for centering the matrix.
}
  \item{signflip}{
logical indicating if the signs of the loadings should be 
flipped such that the absolutely largest value is always positive.
}
\item{via.svd}{
dummy argument for compatibility with classPC calls, will be ignored.
}
  \item{scores}{
logical indicating whether or not scores should be returned.
}
}

\value{
 A list with components: \cr
  \item{rank}{the (numerical) matrix rank of \code{X}, i.e. an integer number
  between 0 and \code{min(dim(x))}.}
  \item{eigenvalues}{the \code{k} eigenvalues, proportional to the variances, where \code{k} is the rank above.}
  \item{loadings}{the loadings, a \eqn{d \times k} matrix.}
  \item{scores}{if the \code{scores} argument was \code{TRUE},
  the \eqn{n \times k} matrix of scores.}
  \item{center}{a vector of means, unless the center argument was \code{FALSE}.}
 \item{scale}{a vector of column scales, unless the scale argument was false.}
}


\author{
P.J. Rousseeuw
}


\seealso{
\code{\link[robustbase]{classPC}}}



\examples{
library(MASS) 
set.seed(12345) 
n <- 100; d <- 10
A <- diag(d) * 0.1 + 0.9
x <- mvrnorm(n, rep(0,d), A)
truncPCA.out <- truncPC(x, ncomp = 2, scores = TRUE)
plot(truncPCA.out$scores)
}
