\name{getKEGGPathway}
\alias{getKEGGPathway}

\title{
  Get a list of genes in a specific KEGG pathway.
}

\description{
  This function gets a list of organisms in the KEGG database.
}
\usage{
  getKEGGPathway(pathway)
}

\arguments{
  \item{pathway}{Pathway name}
}

\details{
  This function delievers a vector with KEGG pathway information.
}
\value{
  A vector with pathways genes and entries. For a overview of all availabl pathway names in a specific organism use the getKEGGPathwayOverview function.
}

\examples{
  \dontrun{
  getKEGGPathway(pathway = "bta00620")
  }
}

\author{
  Daniel Fischer
}

 
