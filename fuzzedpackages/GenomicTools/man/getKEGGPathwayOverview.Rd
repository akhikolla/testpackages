\name{getKEGGPathwayOverview}
\alias{getKEGGPathwayOverview}

\title{
  Get a list of pathways in a specific organism.
}

\description{
  This function gets a list of pathways in a specific organism.
}
\usage{
  getKEGGPathwayOverview(code)
}

\arguments{
  \item{code}{Species code}
}

\details{
  This function delievers a vector with KEGG pathway codes and names for a specific species. For a list of genes in a specific pathway, use the function getKEGGPathway.
}
\value{
  A vector with pathways IDs and names.
}

\examples{
  \dontrun{
  getKEGGPathwayOverview(code = "bta")
  }
}

\author{
  Daniel Fischer
}

 
