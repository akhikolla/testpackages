\name{rejectionPlot}
\alias{rejectionPlot}

\title{
 Creates a Rejection Plot. 
}

\description{
  This function plots the observed and expected numbers of rejected tests and their ratio
  as a function of alpha.
}
\usage{
  rejectionPlot(X,lCol="red",xlim=NULL,crit=NULL,rejLine = NULL,
                alpha = 0.01, rej.lty = c("dotted"), rej.col = "black",
		incLegend=NULL)
}

\arguments{
  \item{X}{Matrix or vector of p-values, as obtained from \code{gmw} function.}
  \item{lCol}{Vector of of line colors}
  \item{xlim}{An interval for alpha values (for the plotting). Default is [0,1].}
  \item{crit}{Shall the distance or the ratio be plotted, see details.}
  \item{rejLine}{Type of rejection line. Options are \code{"bh","bonferroni"} and \code{"simes"}.}
  \item{alpha}{Alpha for the multiple testing adjusted p-value- determines the slope of the rejection line.}
  \item{rej.lty}{Line type of the rejection line.}
  \item{rej.col}{Colour of the rejection line.}
  \item{incLegend}{Position of a possible legend.}
}

\details{
 The matrix \code{X} contains p-values as given from the \code{gmw} function. The rows correspond to different tests,
 and the columns to different variables. 

 The function plots, for each test separately, the observed and expected numbers of rejected tests and their ratio (\code{crit="ratio"})
 or distance (\code{crit="distance"})  as a function of alpha.

}

\references{
  Fischer, D., Oja, H., Schleutker, J., Sen, P.K., Wahlfors, T. (2013): Generalized Mann-Whitney Type Tests for Microarray Experiments, Scandinavian Journal of Statistic, to appear.\cr
 
  Fischer, D., Oja, H. (2013): Mann-Whitney Type Tests for Microarray Experiments: The R Package gMWT, submitted article.
}

\author{
  Daniel Fischer
}

\seealso{
\code{\link{gmw}}
}

\examples{
X <- matrix(c(rnorm(50,2,1),rnorm(60,2,1),rnorm(40,2.5,1)),byrow=TRUE, ncol=10)
colnames(X) <- letters[1:10]
g <- c(rep(1,5),rep(2,6),rep(3,4))


test1 <- gmw(X,g,test="mw",type="external",alternative="greater",goi=c(2,3))
test2 <- gmw(X,g,test="kw",type="external")
test3 <- gmw(X,g,test="jt",type="external",alternative="greater")

tests <- rbind(test1$p.values,test2$p.values,test3$p.values)

rejectionPlot(tests,lCol=c("green","red","blue"))
}
\keyword{methods}
 
