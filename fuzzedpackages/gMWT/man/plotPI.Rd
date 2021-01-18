\name{plotPI}
\alias{plotPI}

\title{
  Creates Scatterplot Matrices for Probabilistic Indices.
}

\description{
  This function creates scatterplot matrices for the different probabilistic indices (PI) \eqn{P_t} , \eqn{P_{tt'}} and \eqn{P_{tt't''}}
  computed for different variables.
}
\usage{
  plotPI(X,g,type="pair",goi=NULL,mc=1,alg="Cnaive",col="black",
           highlight=NULL,hlCol="red",pch=20,zoom=FALSE,order=NULL,...)
}

\arguments{
  \item{X}{Matrix or vector with observations. Each column is a variable, each row an individual.}
  \item{g}{Vector of group labels for observations in \code{X}. Has to be the same length as \code{X} has observations.}
  \item{type}{Type of probabilistic index, see details.}
  \item{goi}{Groups of Interest, see details.}
  \item{mc}{Set the amount of cores to use for parallel calculation (only available for Linux).}
  \item{order}{Boolean, calculate PI only for natural order or for all combinations.}
  \item{alg}{Internal function, which implementation should be used to calculate the PI.}
  \item{col}{Vector of colors of the scatterplot.}
  \item{highlight}{Vector with positions, which are marked in special color.}
  \item{hlCol}{Color of highlighted spots.}
  \item{pch}{Dot type of the plot.}
  \item{zoom}{Logical, shall the plots be zoomed to interesting areas?}
  \item{...}{Additional plotting arguments.}
}

\details{
  This function creates the scatterplot matrices for the PI, in case that they haven't
  been calculated previoulsy. This means that all arguments of the \code{\link{estPI}} are valid here, since
  this function is called first and the results will then be passed to the plot function of the \code{estPI}
  object. 

  See also \code{\link{plot.estPI}} for further details on the specific plot parameters.
}

\references{
  Fischer, D., Oja, H., Schleutker, J., Sen, P.K., Wahlfors, T. (2013): Generalized Mann-Whitney Type Tests for Microarray Experiments, Scandinavian Journal of Statistic, to appear.\cr

  Fischer, D., Oja, H. (2013): Mann-Whitney Type Tests for Microarray Experiments: The R Package gMWT, submitted article.
}

\author{
  Daniel Fischer
}

\seealso{
\code{\link{estPI}}, \code{\link{plot.estPI}}
}

\examples{

X <- c(sample(15))
g <- c(1,1,1,2,2,2,2,3,3,3,4,4,4,4,4)
estPI(X,g,type="single")

X <- matrix(c(rnorm(5000,1.5,2),rnorm(6000,2,2),rnorm(4000,3.5,1)),byrow=TRUE, ncol=10)
colnames(X) <- letters[1:10]
g <- c(rep(1,500),rep(2,600),rep(3,400))

plotPI(X,g,type="single",mc=1)

}
\keyword{methods}
\keyword{hplot}
 
