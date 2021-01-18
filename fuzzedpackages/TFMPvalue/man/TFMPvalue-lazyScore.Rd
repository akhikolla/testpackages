\name{TFMLazyScore}
\alias{TFMLazyScore}

\title{
  Compute the score from P-value.
}
\description{
  Computes the score threshold associated with P-value p using the algorithm of 
  Beckstette 2006.
}
\usage{
TFMLazyScore(mat, pvalue, bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
             type=c("PFM", "PWM"), granularity=1e-5)
}

\arguments{
  \item{mat}{
    The input matrix. It can be a Position Frequency Matrix (PFM) or
    Position Weight Matrix (PWM) in log ratio.
    The matrix must have row names with "A", "C", "G", "T".
  }
  \item{pvalue}{
    The required P-value.
  }
  \item{bg}{
    The background frequency of the sequences. 
    A numeric vector with names "A", "C", "G", "T".
  }
  \item{type}{
    The type of input matrix. Can be "PFM" or "PWM".
  }
  \item{granularity}{
    The granularity used in the computation.
  }
}

\value{
  The score is returned based on the matrix, given P-value and granularity.
}

\author{
  Ge Tan
}

\examples{
  \donttest{
  ## This example is not tested due to running time > 5s
  pfm <- matrix(c(3, 5, 4, 2, 7, 0, 3, 4, 9, 1, 1, 3, 3, 6, 4, 1, 11,
                  0, 3, 0, 11, 0, 2, 1, 11, 0, 2, 1, 3, 3, 2, 6, 4, 1,
                  8, 1, 3, 4, 6, 1, 8, 5, 1, 0, 8, 1, 4, 1, 9, 0, 2, 3,
                  9, 5, 0, 0, 11, 0, 3, 0, 2, 7, 0, 5),
                nrow = 4, dimnames = list(c("A","C","G","T"))
                )
  bg <- c(A=0.25, C=0.25, G=0.25, T=0.25)
  pvalue <- 1e-5
  type <- "PFM"
  granularity <- 1e-5
  TFMLazyScore(pfm, pvalue, bg, type, granularity)
  }
}
