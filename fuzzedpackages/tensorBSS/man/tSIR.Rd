\name{tSIR}
\alias{tSIR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{SIR for Tensor-Valued Observations 
}
\description{
  Computes the tensorial SIR. 
}
\usage{
  tSIR(x, y, h = 10, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numeric array of an order at least three. It is assumed that the last dimension corresponds to the sampling units.
}
  \item{y}{
A numeric or factor response vector.
}
  \item{h}{
The number of slices. If \code{y} is a factor the number of factor levels is automatically used as the number of slices.
}
  \item{\dots}{
Arguments passed on to \code{\link{quantile}}.
}
}
\details{
 Computes the mode-wise sliced inverse regression (SIR) estimators for a tensor-valued data set and a univariate response variable.

}
\value{
  A list with class 'tbss', inheriting from class 'bss', containing the following components: 
  \item{S }{Array of the same size as x containing the predictors.}
  \item{W }{List containing all the unmixing matrices.}
  \item{Xmu }{The data location.}
  \item{datatype}{Character string with value "iid". Relevant for \code{\link{plot.tbss}}.}
}

\author{
  Joni Virta, Klaus Nordhausen
}



\examples{
data(zip.train)
x <- zip.train

rows <- which(x[, 1] == 0 | x[, 1] == 3)
x0 <- x[rows, 2:257]
y0 <- as.factor(x[rows, 1])

x0 <- t(x0)
dim(x0) <- c(16, 16, length(y0))


res <- tSIR(x0, y0)
plot(res$S[1, 1, ], res$S[1, 2, ], col = y0)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{array }
