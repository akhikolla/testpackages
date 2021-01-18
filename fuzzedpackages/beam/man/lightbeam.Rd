\name{lightbeam}
\alias{lightbeam}

\title{Fast inference of a conditional independence graph}

\usage{
lightbeam(X, thres = 0.1, method = "holm", verbose=TRUE)
}

\arguments{
  \item{X}{n by p data matrix}
  \item{thres}{numeric. Significance threshold to be applied on adjusted tail probabilities.}
  \item{method}{character. Method to use for multiple comparison adjustment of tail probabilities.}
  \item{verbose}{logical. Whether information on progress should be be printed.}
}

\value{
An object of class \code{\link{dgCMatrix-class}}
}

\description{
Fast and memory efficient reconstruction of large conditional independence networks.
}

\details{
The function returns a sparse upper triangular matrix containing the (adjusted) tails probabilities that are below the given threshold (i.e. significant).
}

\references{
Leday, G.G.R. and Richardson, S. (2019). Fast Bayesian inference in large Gaussian graphical models. \emph{Biometrics}.
}

\author{
Gwenael G.R. Leday
}

\examples{

# Load data
data(TCPAprad)

# beam
res <- lightbeam(X = TCPAprad, thres=0.1)

}
