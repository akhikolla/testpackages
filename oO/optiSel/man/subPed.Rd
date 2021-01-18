\name{subPed}
\Rdversion{1.1}
\alias{subPed}
\title{Creates a Subset of a Large Pedigree}
\description{
Creates a subset of a large pedigree that includes only individuals related with specified individuals in a predefined way.
}
\usage{
subPed(Pedig, keep, prevGen=3, succGen=0)}

\arguments{
\item{Pedig}{Data frame containing the pedigree where the first 3 columns correspond to: Individual ID,  Sire, and Dam. More columns can be passed in the \code{Pedig} argument including columns named \code{Sex}, \code{Breed} (with breed names), and \code{Born} (with years of birth). Missing parents are coded as NA, 0, or "0".}
\item{keep}{Vector with IDs of individuals. Only these individuals and individuals related with them in a predefined way will be kept in the pedigree.}
\item{prevGen}{Number of previous (ancestral) generations to be included in the pedigree.}
\item{succGen}{Number of succeeding (descendant) generations to be included in the pedigree.}
}

\details{

This function creates a subset of a large pedigree that includes only individuals related with the individuals specified in the vector \code{keep} in a predefined way.
}

\value{
A data frame containing the reduced pedigree. A column \code{keep} is appended indicating which individuals were included in parameter \code{keep}.
}


\examples{

data(PedigWithErrors)

sPed <- subPed(PedigWithErrors, keep="276000891974272", prevGen=3, succGen=2)
sPed

label <- c("Indiv", "Born", "Breed")
pedplot(sPed, mar=c(2,4,2,4), label=label, cex=0.7)

}



\author{Robin Wellmann}
