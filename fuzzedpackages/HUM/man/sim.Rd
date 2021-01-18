\encoding{UTF-8}
\name{sim}
\docType{data}
\alias{sim}
\title{
 desease data
}
\description{
  The data set corresponds to absolute (cells/mm2) or relative (percentage of the cell type in question of the entire inflammatory cell population) densities of 5 major inflammatory cell types in synovial tissue specimens from normal human joints (“Normal”) and from patients with osteoarthritis (“OA”), non-inflammatory orthopedic arthropathies (“Orth.A”), early unclassified arthritis (“EA”), rheumatoid arthritis (“RA”), and chronic septic arthritis (“SeA”). An analysis of this data set with binary and multicategory ROC analysis has been published in Della Beffa PLOS One 2013, which also contains additional details about the data set. The dataset consists of 92 cases with 11 features and disease code.
}
\usage{sim}

\format{A data.frame containing 92 observations of 11 variables.}

\source{
Cristina Della Beffa, Elisabeth Slansky, Claudia Pommerenke, Frank Klawonn, Jialiang Li, Lie Dai, H. Ralph Schumacher Jr., Frank Pessler (2013). The Relative Composition of the Inflammatory Infiltrate as an Additional Tool for Synovial Tissue Classification. PLoS ONE. 8(8): e72494.
}

\seealso{
 \code{\link{dataset}}
}

\examples{
# load the dataset
data(sim)
# CD15
with(sim, by(CD15,Disease,mean))

# CD20
with(sim,tapply(CD20, Disease, FUN = mean))
with(sim, table(CD20=ifelse(CD20<=mean(CD20), "1", "2"), Disease))
}

\keyword{datasets}
