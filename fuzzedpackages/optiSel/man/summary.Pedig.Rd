\name{summary.Pedig}
\Rdversion{1.1}
\alias{summary.Pedig}
\title{Calculates Summary Statistics for Pedigrees.
}
\description{
Calculates summary statistics for pedigrees.
}
\usage{
\method{summary}{Pedig}(object, keep.only=NULL, maxd=50, d=4, ...)}

\arguments{
\item{object}{An object from class \code{Pedig}, which is usually created with function \link{prePed}.}
\item{keep.only}{The individuals to be included in the summary.}
\item{maxd}{Maximum pedigree depth.}
\item{d}{Number of generations taken into account for computing the PCI.}
\item{...}{further arguments passed to or from other methods}
}

\details{
Computes summary statistics for pedigrees, including the numbers of equivalent complete generations, numbers of fully traced generations, numbers of maximum generations traced, indexes of pedigree completeness (MacCluer et al, 1983), and the inbreeding coefficients.
}

\value{
A data frame with the following columns:
\tabular{ll}{
\code{Indiv} \tab IDs of the individuals, \cr
\code{equiGen} \tab Number of equivalent complete generations, \cr
\code{fullGen} \tab Number of fully traced generations, \cr
\code{maxGen} \tab Number of maximum generations traced, \cr
\code{PCI} \tab Index of pedigree completeness (MacCluer et al, 1983) in generation \code{d}. \cr
\code{Inbreeding} \tab Inbreeding coefficient. \cr
}
}


\examples{
data(PedigWithErrors)
Pedig <- prePed(PedigWithErrors)
Summary <- summary(Pedig, keep.only=Pedig$Born \%in\% (2006:2007))
head(Summary)

hist(Summary$PCI,        xlim=c(0,1),  main="Pedigree Completeness")
hist(Summary$Inbreeding, xlim=c(0,1),  main="Inbreeding")
hist(Summary$equiGen,    xlim=c(0,20), main="Number of Equivalent Complete Generations")
hist(Summary$fullGen,    xlim=c(0,20), main="Number of Fully Traced Generations")
hist(Summary$maxGen,     xlim=c(0,20), main="Number of Maximum Generations Traced")
}

\references{
MacCluer J W, Boyce A J, Dyke B, Weitkamp L R, Pfenning D W, Parsons C J (1983). Inbreeding and pedigree structure in Standardbred horses. J Hered 74 (6): 394-399. 
}

\author{Robin Wellmann}
