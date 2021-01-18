\name{conttac}
\Rdversion{1.1}
\alias{conttac}
\title{Calculates \bold{Cont}ributions \bold{T}o \bold{A}ge \bold{C}ohorts
}
\description{
Calculates genetic contributions of other breeds to age cohorts
}
\usage{
conttac(cont, cohort, use=rep(TRUE,length(cohort)), mincont=0.05, long=TRUE)}

\arguments{
\item{cont}{Data frame containing the genetic contributions of several ancestors or breeds to all individuals. This is typically the output of function \link{pedBreedComp}.}
\item{cohort}{Numeric vector indicating for every individual the age cohort to which it belongs (typically year of birth).}
\item{use}{Logical vector indicating for every individual whether it should be included in an age cohort (typically \code{TRUE} for individuals belonging to the breed of interest).}
\item{mincont}{Contributions of breeeds with average contribution smaller than \code{mincont} will be summarized in one row}
\item{long}{Should the resutling data frame be melted for easy plotting?}
}

\details{
The genetic contributions from other breeds to all age cohorts are computed. The genetic contribution from a breed is the fraction of genes in the gene pool originating from the respective breed.
}

\value{
Data frame containing the genetic contribution from every breed to every age cohort.
}


\examples{
data(ExamplePed)
Pedig      <- prePed(ExamplePed, thisBreed="Hinterwaelder", lastNative=1970)
cont       <- pedBreedComp(Pedig, thisBreed="Hinterwaelder")
contByYear <- conttac(cont, Pedig$Born, use=Pedig$Breed=="Hinterwaelder", mincont=0.04, long=FALSE)
round(contByYear,2)

barplot(contByYear, ylim=c(0,1), col=1:10, ylab="genetic contribution",
        legend=TRUE, args.legend=list(x="bottomleft",cex=0.5))

}



\author{Robin Wellmann}
