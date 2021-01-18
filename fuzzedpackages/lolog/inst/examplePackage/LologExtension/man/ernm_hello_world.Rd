\name{ernm_hello_world}
\alias{ernm_hello_world}
\title{An basic example of a function in C++ using ernm}
\usage{
  ernm_hello_world()
}
\value{
  a list of a character vector, a numeric vector, and an
  lolog DirectedNet
}
\description{
  An basic example of a function in C++ using ernm
}
\examples{
lolog_hello_world()

#MinDegree is an new statistic defined by this package
if(require(network)){
	data(flo)
	nflo<-as.BinaryNet(network(flo,directed=FALSE) )
	fit <- lolog(nflo ~ edges() + minDegree(3),verbose=0)
	summary(fit)
}
}

