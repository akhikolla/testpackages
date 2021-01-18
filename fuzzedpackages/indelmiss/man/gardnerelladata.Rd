\name{gardnerelladata}
\alias{gardnerelladata}
\docType{data}
\title{
Gardnerella vaginalis phyletic data (gene family memberships) and phylogenetic tree
}
\description{
Gene family memberships for thirty five species from the genus \emph{Gardnerella vaginalis}.
}
\usage{data("gardnerelladata")}
\format{
Contains a list that comprises a tree (called "tree") and phyletic gene family membership data (called "phyl") as its components. The tree is in the ape package phylo format. The data component consists of a matrix of 0/1 patterns with the different patterns as the rows and the 35 taxa as the columns. An entry of 1 (0) describe gene presence (absence) in the taxa.
  }
\examples{
data(gardnerelladata)
}
\references{  
Devault, A. (2014). Genomics of Ancient Pathogenic Bacteria: Novel Techniques & Extraordinary Substrates. Ph. D. thesis, McMaster University, Hamilton.

Devault, A. M., T. D. Mortimer, H. Kiesewetter, T. Smith, M. Kuch, A. Hussain, J. M. Enk, G. B. Golding, W. Aylward, H. N. Poinar, and C. S. Pepperell (2016). Ancient DNA analysis of a unique calcified urogenital infection from medieval Troy reveals insights into past female health. Submitted.
}
\keyword{datasets}