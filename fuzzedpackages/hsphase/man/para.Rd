\name{para}
\alias{para}

\title{
Parallel Analysis of Data
}
\description{
This function uses the list of matrices (the output of \code{\link{cs}}) and runs one of the options, on each element of the list, in parallel. 
}
\usage{
para(halfsibs, cpus = 1, option = "bmh", type = "SOCK", bmh_forwardVectorSize = 30,
bmh_excludeFP = TRUE, bmh_nsap = 3, pmMethod = "constant")
}
\keyword{opposing-homozygote}

\arguments{
  \item{halfsibs}{
\code{list} list of matrices of half-sibs (can be generated with \code{\link{hss}} and \code{\link{cs}} functions)
}
  \item{cpus}{
 \code{numeric} number of CPUs (thread)
}
  \item{option}{
 \code{character}  type of analysis
}
  \item{type}{\code{character}
  type of cluster for parallel analysis

}

\item{bmh_forwardVectorSize}{
  \code{integer} number of heterozygous sites used to validate recombination events or check for genotyping errors}


\item{bmh_excludeFP}{
  \code{logical} exclude SNPs that may cause heterozygous sites in the sire due to genotyping errors or map errors}


\item{bmh_nsap}{
  \code{integer} number of SNPs per block}

\item{pmMethod}{
  \code{character} method for creating the recombination matrix}
}
\details{
Type of analysis  can be  \code{\link{bmh}}, \code{\link{ssp}}, \code{\link{aio}}, \code{\link{pm}}, or rec (refer to \code{\link{pm}}, \code{\link{rplot}} and vignette for more information about rec).

}
\value{
Returns a list of matrices with the results (formats specific to the option selected).
}
\examples{
# Please run demo(hsphase)

}

