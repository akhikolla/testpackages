\name{genecont}
\Rdversion{1.1}
\alias{genecont}
\title{Calculates Genetic Contributions using Pedigrees. 
}
\description{
Calculates the genetic contributions each individual has from specified ancestors.
}
\usage{
genecont(Pedig, from=NULL, to=NULL)
}


\arguments{
\item{Pedig}{Data frame containing the pedigree, where the first columns are \code{Indiv} (Individual ID),  \code{Sire}, and \code{Dam}.}
\item{from}{Vector with ancestors whose contributions to the individuals should be calculated. By default, the contributions from all individuals will be calculated.}
\item{to}{Vector with individuals for which the contributions from ancestors should be calculated. By default, the contributions are calculated for all individuals.}
}



\details{
This function calculates genetic contributions of specified ancestors to each individual.
}

\value{
Lower triangular matrix with genetic contributions for each pair of individuals. Column i contains the genetic contribution of ancestor i to all individuals.
}


\examples{

data(ExamplePed)
Pedig <- prePed(ExamplePed)
cont  <- genecont(Pedig)

plot(Pedig$Born, cont[,"276000803611144"], pch=18, ylim=c(0,1))
Pedig["276000803611144",]

#faster:
cont  <- genecont(Pedig, from="276000803611144")
head(cont)
plot(Pedig$Born, cont[,"276000803611144"], pch=18, ylim=c(0,1))
}



\author{Robin Wellmann}
