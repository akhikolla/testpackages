\name{prior_unif_dirichlet}
\alias{prior_unif_dirichlet}

\title{
Uniform prior counts for the hyper-Dirichlet distribution.
}
\description{
The function generates a \eqn{(r,r,p,p)}-array filled with uniform counts for the hyper-Dirichlet distribution used as prior in \link{lweights_multinomial} when there are \eqn{p} variables with \eqn{r} levels.
\code{Neq} is the equivalent prior sample size.
}
\usage{
prior_unif_dirichlet(p, r, Neq = 0.5 * r^2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
Number of variables.
}
  \item{r}{
Number of levels.
}
  \item{Neq}{
Equivalent sample size.
}
}

\value{
\item{prior}{A \eqn{(r,r,p,p)}-array containings counts.}
}
\references{
This package implements the method described in the paper "Bayesian Inference of Graphical Model Structures Using Trees" by L. Schwaller, S. Robin, M. Stumpf, 2015 (submitted and availavable on arXiv).
}
\author{
Loïc Schwaller
}

\examples{
library('saturnin')
p <- 100
r <- 10
prior <- prior_unif_dirichlet(p,r)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hyper dirichlet }

