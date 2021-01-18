\name{createGroups}
\alias{createGroups}

\title{
 Renumbering the Group Labels 
}

\description{
  This function transforms the k group labels to integers 1,...,k in a specified order.
}
\usage{
  createGroups(g, desOrder)
}

\arguments{
  \item{g}{Vector with group labels.}
  \item{desOrder}{Specified order for the groups.}
}

\details{
 This function helps to renumber the group labels in a specified way. 
 The vector \code{desOrder} gives the old numbers as the function of the new numbers.
 The first element of \code{desOrder} gives the old number that is transformed to be one, and so on.
 The renumbering is used to specify alternative hypotheses.

 The output is then a transformed vector  \code{g}  with new group numbers.

 For non-numeric old group labels, the old number is the rank in the lexicographical order.
}
\value{
  A vector containing the new group numbers.
}

\author{
  Daniel Fischer
}

\examples{
g <- c(rep(1,5),rep(2,6),rep(3,4))
createGroups(g,c(2,3,1))

# 'Smallest' letter should become 3, second smallest 2 and largest 1:
g <- sample(c("A","B","C"),10,replace=TRUE)
createGroups(g,c(3,2,1))
}

\keyword{methods}
 
