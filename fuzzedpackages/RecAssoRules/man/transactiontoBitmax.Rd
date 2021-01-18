\name{transactiontoBitmax}
\alias{transactiontoBitmax}
\title{Creation of a Binary DataFRame from transactions type}
\usage{
transactiontoBitmax(transac , deli)
}
\arguments{
\item{transac}{a character vector or matrix.}

\item{deli}{The character used to delimit the objects present in a transaction.}

}

\description{
Returns the corresponding Binary R DataFrame. If you imported your transaction as R DataFrame, please use the R function as.matrix(yourdataframe).
}

\examples{

data("Transacfruits")
binaryfruits = transactiontoBitmax (as.matrix(Transacfruits),",")
## Import transaction data Transacfruits and create the correspondant binary matrix.
}
