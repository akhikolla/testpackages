\name{estPI}
\alias{estPI}

\title{
  Calculate Probabilistic Index
}

\description{
  This function calculates the probabilistic indices \eqn{\hat{P}_t} , \eqn{\hat{P}_{tt'}} and \eqn{\hat{P}_{tt't''}} 
  to compare the groups of observations.
}
\usage{
  estPI(X, g, type="pair", goi=NULL, mc=1, order=TRUE, alg="Cnaive")
}

\arguments{
  \item{X}{Matrix or vector with observations. The rows refer to individuals, the columns to variables.}
  \item{g}{Vector of group numbers for observations in \code{X}. Its length has to be the number of observations in \code{X}.}
  \item{type}{Type of probabilistic index, see details.}
  \item{goi}{Groups Of Interest, see details.}
  \item{mc}{Multiple Cores, set the amount of cores to use for parallel calculation (only available for Linux).}
  \item{order}{Boolean, calculate probabilistic index only for a specified order of groups or for all possible orders.}
  \item{alg}{Internal function, which algorithm should be used to calculate the probabilistic index.}
}

\details{
  The matrix \code{X} contains the data. Each column refers to a variable, each row to an observation. 
  The group memberships of the observations are given in \code{g}. In the case of one dimensional data,  
 \code{X} is a vector.

  There is also an option to calculate the probabilistic indices (PI) only for a subset of the groups. 
  In that case those group labels can be specified with the \code{goi} option by giving the corresponding
  group labels. 

  Different types of PI can be calculated: \code{"single"} calculates the probability \eqn{\hat{P}_t} for each group,
  \code{"pair"} produces the probabilistic indices \eqn{\hat{P}_{tt'}} for all pairs of groups \eqn{t<t'}, and
  \code{"triple"} provides the probabilities \eqn{\hat{P}_{tt't''}} for all triples of groups \eqn{t<t'<t''}.
  See Fischer et al. (2013) for more details.

  Note that the PIs are calculated using the group numbering as given in \code{g}. See also the function
  \code{\link[gMWT]{createGroups}} for renumbering the group labels. By specifying the option \code{order=FALSE}
  the PIs for all possible group orders will be calculated. The default is that the PI is only calculated for the
  natural order given by \code{g}.

  In case this code is executed on a Linux machine and \code{X} is a data matrix, the calculation can be parallelized
  using the option \code{mc} to specify the amount of used calculation cores. 

  Different algorithms are available for the  calculation of the PIs. The default is the fastest
  possible option and the user should not change this option. Different algorithms are only provided for validation and testing
  purposes. Options here are currently \code{Cnaive}, \code{Rnaive}, \code{Rgrid}, \code{Rsub}, \code{Csub} but not all combinations
  of \code{type}/\code{alg} are available.
}
\value{
  A list with class 'estPI' containing the following components:
  \item{probs}{Matrix or vector of the  PIs.}
  \item{type}{String, storing the type of PI.}
  \item{goi}{Vector, the Groups Of Interest, as given in the function call.}
  \item{order}{Boolean, PI just for the specified order or for all orders.}
  \item{obs}{Matrix, the original data matrix.}
  \item{g}{Vector, the original group vector.}
  \item{alg}{String, the requested algorithm.}
}

\references{
  Fischer, D., Oja, H., Schleutker, Sen, P.K., J., Wahlfors, T. (2013): Generalized Mann-Whitney Type Tests for Microarray Experiments, Scandinavian Journal of Statistics, to appear.\cr

  Fischer, D., Oja, H. (2013): Mann-Whitney Type Tests for Microarray Experiments: The R Package gMWT, submitted article.
}

\author{
  Daniel Fischer
}

\seealso{
\code{\link[gMWT]{createGroups}}
}

\examples{

X <- c(sample(15))
g <- c(1,1,1,2,2,2,2,3,3,3,4,4,4,4,4)
estPI(X,g,type="single")

X <- matrix(c(rnorm(5000,1.5,2),rnorm(6000,2,2),rnorm(4000,3.5,1)),byrow=TRUE, ncol=10)
colnames(X) <- letters[1:10]
g <- c(rep(1,500),rep(2,600),rep(3,400))

estPI(X,g,type="single",mc=1)

}
\keyword{methods}
 
