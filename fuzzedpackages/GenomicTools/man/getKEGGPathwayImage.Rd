\name{getKEGGPathwayImage}
\alias{getKEGGPathwayImage}

\title{
  Get a visualization of a specific KEGG pathway.
}

\description{
  This function gets a visualizatin of a specific KEGG pathway.
}
\usage{
  getKEGGPathwayImage(pathway, folder=NULL)
}

\arguments{
  \item{pathway}{Pathway name}
  \item{folder}{Path to where the image is stored}
}

\details{
  thi function downloads a png figure of a specifi pathway and stores it in a given folder
}
\value{
  An image.
}

\examples{
 \dontrun{
  getKEGGPathwayImage(pathway = "bta00620")
  }
}

\author{
  Daniel Fischer
}

 
