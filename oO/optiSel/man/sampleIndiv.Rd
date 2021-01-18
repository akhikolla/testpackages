\name{sampleIndiv}
\Rdversion{1.1}
\alias{sampleIndiv}
\title{Sample Individuals from Pedigree}

\description{
Sampling Individuals from a Pedigree.
}

\usage{
sampleIndiv(Pedig, from="Born", each=100)
}

\arguments{
\item{Pedig}{Pedigree with column \code{Indiv} and the column specified in parameter \code{from}.}
\item{from}{Column name. From each cohort specified in this column (e.g. year of birth), the number of individuals specified in parameter \code{each} is sampled. If a cohort contains less individuals, then all individuals are chosen.}
\item{each}{Number of individuals to be sampled from each cohort.}
}

\details{
From each cohort, a specified number of individuals will be sampled. If a cohort contains less individuals, then all individuals are sampled.
This may be needed for estimating population specific parameters from a subset of a large pedigree to reduce computation time.
}

\value{
Character vector containing the IDs of the individuals.
}

\examples{
data("PedigWithErrors")
set.seed(1)
Pedig <- prePed(PedigWithErrors)
use   <- Pedig$Breed=="Hinterwaelder"
keep  <- sampleIndiv(Pedig[use, ], from="Born", each=5)
keep
}

\author{Robin Wellmann}
