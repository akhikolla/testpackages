\name{eigenMapMatMult}
\alias{eigenMapMatMult}
\docType{package}
\title{
Matrix multiplication using RcppEigen
}
\description{
Matrix multiplication using RcppEigen.
}

\usage{
eigenMapMatMult(A,B)
}
\arguments{
  \item{A, B}{
    numeric (double) complex matrices or vectors.
  }

}

\value{
  The matrix product. The value is the same as A \%*\% B
}

\examples{
library(SMUT)
A=matrix(1:9,3,3)
A=A+0
B=as.matrix(c(5.0, 2.0, 0.0))
eigenMapMatMult(A,B)
# the result is the same as A \%*\% B

# Thanks for using our R package SMUT 
}

