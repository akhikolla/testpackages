\name{TFMpv2sc}
\alias{TFMpv2sc}

\title{
  Compute score from P-value.
}
\description{
  Computes the score threshold associated with a P-value.
}
\usage{
TFMpv2sc(mat, pvalue, bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
         type=c("PFM", "PWM"))
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
}

\value{
  The score is returned based on the matrix, given P-value.
}

\references{
  Touzet, H., and Varre, J.-S. (2007). Efficient and accurate P-value computation for Position Weight Matrices. Algorithms Mol Biol 2, 15.

}

\author{
  Ge Tan
}

\examples{
  pfm <- matrix(c(3, 5, 4, 2, 7, 0, 3, 4, 9, 1, 1, 3, 3, 6, 4, 1, 11,
                  0, 3, 0, 11, 0, 2, 1, 11, 0, 2, 1, 3, 3, 2, 6, 4, 1,
                  8, 1, 3, 4, 6, 1, 8, 5, 1, 0, 8, 1, 4, 1, 9, 0, 2, 3,
                  9, 5, 0, 0, 11, 0, 3, 0, 2, 7, 0, 5),
                nrow = 4, dimnames = list(c("A","C","G","T"))
                )
  bg <- c(A=0.25, C=0.25, G=0.25, T=0.25)
  pvalue <- 1e-5
  type <- "PFM"
  score <- TFMpv2sc(pfm, pvalue, bg, type)
}

