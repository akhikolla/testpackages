\name{q1LowRank}
\alias{q1LowRank}
\encoding{UTF8}
\title{

  Qualitative Correlation or Covariance Kernel with one Input and
  Low-Rank Correlation

}
\usage{

q1LowRank(factor, rank = 2L, input = "x",
          cov = c("corr", "homo", "hete"), intAsChar = TRUE)

}
\arguments{
  \item{factor}{

    A factor with the wanted levels for the covariance kernel object.

  }
  \item{rank}{

    The wanted rank, which must be \eqn{\geq 2}{>= 2} and \eqn{< m}
    where \eqn{m} is the number of levels.
    
  }
  \item{input}{

    Name of (qualitative) input for the kernel.

  }
  \item{cov}{

    Character telling what variance structure will be chosen:
    \emph{correlation} with no variance parameter, \emph{homoscedastic}
    with one variance parameter or \emph{heteroscedastic} with \eqn{m}
    variance parameters.
    
  }
  \item{intAsChar}{
    
    Logical. If \code{TRUE} (default), an integer-valued input will be
    coerced into a character. Otherwise, it will be coerced into a factor.
    
  }
}
\value{

  An object with class \code{"covQual"} with \code{d = 1} qualitative
  input.

}
\description{
  
  Qualitative correlation or covariance kernel with one input and
  low-rank correlation.
  
}
\details{

  The correlation structure involves \eqn{(r - 1)(m - r /2)} parameters.
  The parameterization of Rapisarda et al is used: the correlation
  parameters are angles \eqn{\theta_{i,j}}{Theta[i, j]} corresponding to
  \eqn{1 < i \leq r}{1 < i <= r} and \eqn{1 \leq j < i}{1 < j < i}
  or to \eqn{r < i \leq m}{r < i <= m} and
  \eqn{1 \leq j < r}{1 <)= j < r}. The
  correlation matrix \eqn{\mathbf{C}}{C} for the levels, with size
  \eqn{m}, factors as
  \eqn{\mathbf{C} = \mathbf{L}\mathbf{L}^\top}{C = L \%*\% t(L)}
  where \eqn{\mathbf{L}}{L} is a lower-triangular
  matrix with dimension \eqn{m \times r}{c(m, r)} with all its rows
  having unit Euclidean norm. Note that the diagonal elements of
  \eqn{\mathbf{L}}{L} can be negative and correspondingly the angles
  \eqn{\theta_{i,1}}{Theta[i, 1]} are taken in the interval
  \eqn{[0, 2\pi)}{[0, 2*pi]} for \eqn{1 < i \leq r}{1 < i <= r}. The
  matrix \eqn{\mathbf{L}}{L} is not unique. As explained in Grubišić and
  Pietersz, the parameterization is surjective: any correlation with
  rank \eqn{\leq r}{<= r} is obtained by choosing a suitable vector
  of parameters, but this vector is not unique.
  
  Correlation kernels are needed in tensor products because the tensor
  product of two covariance kernels each with unknown variance would not
  be identifiable.

}

\references{

  Francesco Rapisarda, Damanio Brigo, Fabio Mercurio (2007).
  "Parameterizing Correlations a Geometric Interpretation".
  \emph{IMA Journal of Management Mathematics}, \bold{18}(1):
  55-73.
  
  Igor Grubišić, Raoul Pietersz
  (2007). "Efficient Rank Reduction of Correlation Matrices". \emph{Linear
  Algebra and its Applications}, \bold{422}: 629-653.

}
\seealso{

  The \code{\link{q1Symm}} function to create a kernel object for the
  full-rank case and \code{\link{corLevLowRank}} for the correlation
  function.

}
\examples{
myFact <- factor(letters[1:8])
myCov <- q1LowRank(factor = myFact, rank = 3)
## corrplot
plot(myCov)
## find the rank using a pivoted Cholesky
chol(covMat(myCov), pivot = TRUE)
}