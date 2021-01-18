\name{completeness}
\Rdversion{1.1}
\alias{completeness}
\title{Calculates Pedigree Completeness}

\description{
Calculates completeness of the pedigree for individuals and for groups of individuals in each ancestral generation.
}
\usage{
completeness(Pedig, keep=NULL, maxd=50, by="Indiv")}

\arguments{
\item{Pedig}{Data frame containing the pedigree, where the first columns are \code{Indiv} (Individual ID),  \code{Sire}, and \code{Dam}. More columns can be passed in the \code{Pedig} argument, in particular a column for grouping with the name defined by argument \code{by}.}
\item{keep}{Vector with IDs of the individuals for which the completeness will be calculated, or a logical vector indicating the individuals. By default, all individuals are used.}
\item{maxd}{Number of generations for which completeness should be calculated.}
\item{by}{Name of a column in data frame \code{Pedig}. The completeness will be computed separately for each group defined by the column.}
}

\details{
The function computes the completeness of the pedigree for the specified individuals and for groups of individuals. It is the proportion of known ancestors in each generation. Generation 0 corresponds to the individual itself, so the completeness is always 1 in generation 0.
}

\value{

Data frame with the following columns

\item{Indiv (or 'by')}{ID of the individual or level of the grouping factor,}
\item{Generation}{Generation number,}
\item{Completeness}{Completeness of the pedigree in the respective generation.}
}

\seealso{
Another function for characterizing pedigree completeness is \code{\link{summary.Pedig}}.
}

\examples{

#Computes the pedigree completeness of Hinterwald cattle 
#born between 2006 and 2007 in each ancestral generation.

data(PedigWithErrors)
Pedig <- prePed(PedigWithErrors)
compl <- completeness(Pedig, keep=Pedig$Born \%in\% (2006:2007), maxd=50, by="Indiv")
head(compl)

#Summary statistics can be computed directly from the pedigree:
Summary <- summary(Pedig, keep=Pedig$Born \%in\% (2006:2007))
head(Summary)

hist(Summary$PCI,        xlim=c(0,1),  main="Pedigree Completeness")
hist(Summary$Inbreeding, xlim=c(0,1),  main="Inbreeding")
hist(Summary$equiGen,    xlim=c(0,20), main="Number of Equivalent Complete Generations")
hist(Summary$fullGen,    xlim=c(0,20), main="Number of Fully Traced Generations")
hist(Summary$maxGen,     xlim=c(0,20), main="Number of Maximum Generations Traced")

compl <- completeness(Pedig, keep=Pedig$Born \%in\% (2006:2007), maxd=50, by="Sex")
head(compl)

\dontrun{
library("ggplot2")
ggplot(compl, aes(Generation, Completeness, col=Sex))+geom_line()
}
}

\references{Cazes P, Cazes MH. (1996) Comment mesurer la profondeur genealogique d'une ascendance? Population (French Ed) 51:117-140.
}

\author{Robin Wellmann}
