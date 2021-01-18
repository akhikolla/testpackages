\name{prefrecrules}
\alias{prefrecrules}
\title{Mining for frequent patterns and confident rules.}
\usage{
prefrecrules (Bitmax  , varnames , relativeSup , Minconf )
}
\arguments{
\item{Bitmax}{A binary DataFrame.}

\item{varnames}{A vector of character or a Matrix character.}

\item{relativeSup}{The wanted relative min Support value. }

\item{Minconf}{The wanted min confident value. }
}

\description{
Returns a list with two objets : freqset a R DataFrame with 3 columns , the set-names , their support and relative support and rulesset a R Dataframe with 3 colums , the antecedant names , the consequent names and their confidence. When lauching the function, the BinaryMatrix is ordered in ascendant colsums order. For each column, the frequent itemSets and confident Rules containing the current tested item are extracted with the ARM traversal function.
}

\examples{

data("Transacfruits")
## Import and build Transafruits as binary matrix.
binaryfruits = transactiontoBitmax (as.matrix(Transacfruits),",")
outsetrules = prefrecrules (binaryfruits, names(binaryfruits),0.02, 0.4)
## Mine frequents itemsets and confidentes rules with :
## relative min support of 0.02 and min confidence of 0.4.
}
