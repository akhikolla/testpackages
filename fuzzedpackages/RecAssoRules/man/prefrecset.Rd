\name{prefrecset}
\alias{prefrecset}
\title{Mining for frequent patterns. }
\usage{
prefrecset (Bitmax , varnames , relativeSup )
}
\arguments{
\item{Bitmax}{A binary DataFrame}

\item{varnames}{A vector of character or a Matrix character.}

\item{relativeSup}{The wanted relative min Support value.}
}

\description{
Provides a R DataFrame with 3 columns , the set-names , their support and relative support. When lauching the function, The BinaryMatrix is ordered in ascendant colsums order. For each column, the frequent itemSets containing the current tested item are extracted with the traversal function.
}

\examples{

data("Transacfruits")
## Import and build Transafruits as binary matrix.
binaryfruits = transactiontoBitmax (as.matrix(Transacfruits),",")
itemsets = prefrecset (binaryfruits, names(binaryfruits), 0.02)
## Mine frequents itemSet from Transacfruits with a relative min support of 0.02.
}
