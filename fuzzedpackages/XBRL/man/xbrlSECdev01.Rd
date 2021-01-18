\name{xbrlSECdev01}
\alias{xbrlSECdev01}
\title{"Mutable state" function that exposes "methods" that analyze a
  US SEC submission.}
\description{
  \code{xbrlSECtmp} is a "mutable state" function (see
  \code{\link{XBRL}}). It exposes the method
  \code{showPresentationHierarchy}. It is a work in progress.
}
\usage{
xbrlSECdev01(xbrl.vars)
}
\arguments{
  \item{xbrl.vars}{list of variables returned by \code{\link{xbrlDoAll}}}
}
\value{a list of "methods"}
\details{
  \itemize{
    \item{showPresentationHierarchy(showLabels=TRUE, showFacts=FALSE,
      file="")}{ Shows a tree of the presentations provided.}
  }
}
\author{
Roberto Bertolusso and Marek Kimmel
}
\seealso{
  \code{\link{xbrlDoAll}}, \code{\link{xbrlParse}}
}
\examples{
\dontrun{
## Setting stringsAsFactors = FALSE is highly recommended
## to avoid data frames to create factors from character vectors.
options(stringsAsFactors = FALSE)

## XBRL instance file to be analyzed, accessed
## directly from SEC website:
inst <- "https://www.sec.gov/Archives/edgar/data/21344/000002134413000050/ko-20130927.xml"

xbrl.vars <- xbrlDoAll(inst)
xbrl.sec <- xbrlSECdev01(xbrl.vars)
xbrl.sec$showPresentationHierarchy(showLabels=TRUE, showFacts=TRUE)
}
}
\keyword{XBRL}
\keyword{mid-level}

