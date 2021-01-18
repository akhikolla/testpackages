\name{getKEGGOrganisms}
\alias{getKEGGOrganisms}

\title{
  Get a list of organisms in the KEGG database
}

\description{
  This function gets a list of organisms in the KEGG database.
}
\usage{
  getKEGGOrganisms(url="http://rest.kegg.jp/list/organism")
}

\arguments{
  \item{url}{Link to KEGG database.}
}

\details{
  This function extracts a list of KEGG database organism names.
}
\value{
  A data.frame with organism names.
}

\author{
  Daniel Fischer
}

 
