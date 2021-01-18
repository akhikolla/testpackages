\name{nroRcppVector}
\alias{nroRcppVector}

\title{Safety check for Rcpp calls}

\description{Ensures vectors can be passed safely to C++ functions.}

\usage{
nroRcppVector(data, default, numeric = TRUE, empty = TRUE)
}

\arguments{
  \item{data}{
A vector.
  }
  
  \item{default}{
Default output if input is empty.
  }
  
  \item{numeric}{
If TRUE, output is converted into a numeric vector, otherwise a
character vector is returned.
  }

  \item{empty}{
If TRUE, empty output is allowed.
  }
}

\details{
Checks the input for size (0 is allowed if empty flag is set) and that it
is a vector (or can be converted to a vector).
}

\value{
A numeric or a character vector.
}

\author{
Ville-Petteri Makinen
}

\examples{
# Empty input reverts to default.
x <- nroRcppVector(data=NULL, default=NA, empty=TRUE)
print(x)

# Empty input reverts to default, then to specified type.
x <- nroRcppVector(data=NULL, default=123, empty=TRUE, numeric=FALSE)
print(x)

# Convert a logical vector to numbers.
x <- c(TRUE, TRUE, FALSE, TRUE)
names(x) <- c("a","b","c","d")
y <- nroRcppVector(data=x, numeric=TRUE)
print(y)
}
