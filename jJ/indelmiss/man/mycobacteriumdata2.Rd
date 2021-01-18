\name{mycobacteriumdata2}
\alias{mycobacteriumdata2}
\docType{data}
\title{
Alternate Mycobacterium tree.
}
\description{
Gene family memberships and a phylogenetic tree for ten species from the genus \emph{Mycobacterium}.
}
\usage{data("mycobacteriumdata2")}
\format{
  Contains a list that comprises a tree (called "tree") and phyletic gene family membership data (called "phyl") as its components. The tree is in the ape package phylo format. The data component consists of a data frame of 0/1 patterns with the different patterns as the rows and the 10 taxa as the columns. An entry of 1 (0) describe gene presence (absence) in the taxa.
}
\examples{
data(mycobacteriumdata2)
}
\keyword{datasets}