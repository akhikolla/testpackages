\name{nroRcppMatrix}
\alias{nroRcppMatrix}

\title{Safety check for Rcpp calls}

\description{Forces all values to numeric to be passed to C++ functions.}

\usage{
nroRcppMatrix(data, trim)
}

\arguments{
  \item{data}{
A matrix or a data frame.
  }
  
  \item{trim}{
if TRUE, empty rows and columns are removed.
  }
}

\details{
Converts all columns to values that have a numeric representation.
Detects columns that can be represented as 0s and 1s.
}

\value{
A matrix or a data frame with the attribute 'binary' that contains the
names of binary columns and 'excl.rows' and 'excl.columns' contain the names
of rows and columns that were excluded.
}

\author{
Ville-Petteri Makinen
}

\examples{
# Fully numeric data frame.
x <- data.frame(A=c(1,2,3,4), B=c(0,1,0,NA), C=c(2,3,4,5))
print(nroRcppMatrix(data=x, trim=TRUE))

# Matrix of characters, some of which can be converted to numbers.
x <- matrix(c("1","2","b","4","","6","7","8"), nrow=4, ncol=2)
print(nroRcppMatrix(data=x, trim=TRUE))

# Object that can be converted to numbers.
x <- list(text="abc", value="123")
print(nroRcppMatrix(data=x, trim=TRUE))

# Unusable object.
x <- list(text="abc", value="123", multiple=c("a","b","c"))
print(nroRcppMatrix(data=x, trim=TRUE))
}
