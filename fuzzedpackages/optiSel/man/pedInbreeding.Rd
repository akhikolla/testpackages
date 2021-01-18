\name{pedInbreeding}
\Rdversion{1.0}
\alias{pedInbreeding}
\title{Calculates Pedigree Based Inbreeding 
}
\description{                                                 
Calculates Pedigree Based Inbreeding 
}
\usage{
pedInbreeding(Pedig)}

\arguments{
\item{Pedig}{Data frame containing the Pedigree with the first 3 columns being \code{Indiv} (individual ID), \code{Sire}, and \code{Dam}, which is usually obtained with function \link{prePed}. Missing parents are coded as NA.}
}

\details{
Computation of pedigree based inbreeding.
This function is a wrapper function for \code{\link[pedigree:calcInbreeding]{pedigree}} from package  \code{\link[pedigree:pedigree-package]{pedigree}}. 
}

\value{
A data frame with column \code{Indiv} containing the individual IDs and column \code{Inbr} containing the inbreeding coefficients.
}


\examples{
data(PedigWithErrors)
data(Phen)
keep  <- Phen$Indiv 
Pedig <- prePed(PedigWithErrors, keep=keep)
Res   <- pedInbreeding(Pedig)
mean(Res$Inbr[Res$Indiv \%in\% keep])
#[1] 0.01943394
}



\author{Robin Wellmann}
