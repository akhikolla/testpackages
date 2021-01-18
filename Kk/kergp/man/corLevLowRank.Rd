\name{corLevLowRank}
\alias{corLevLowRank}
\encoding{UTF8}
\title{

  Correlation Matrix for a Low-Rank Structure

}
\usage{

corLevLowRank(par, nlevels, rank, levels,
              lowerSQRT = FALSE, compGrad = TRUE,
              cov = 0, impl = c("C", "R"))

}
\arguments{
  \item{par}{

    A numeric vector with length \code{npCor + npVar} where
    \code{npCor = (rank - 1) *} \code{(nlevels - rank / 2)} is the number of
    correlation parameters, and \code{npVar} is the number of variance
    parameters, which depends on the value of \code{cov}. The value of
    \code{npVar} is \code{0}, \code{1} or \code{nlevels} corresponding
    to the values of \code{cov}: \code{0}, \code{1} and \code{2}.  The
    correlation parameters are assumed to be located at the head of
    \code{par} i.e. at indices \code{1} to \code{npCor}. The variance
    parameter(s) are assumed to be at the tail, i.e. at indices
    \code{npCor +1 } to \code{npCor + npVar}.

  }
  \item{nlevels}{

    Number of levels \eqn{m}.

  }
  \item{rank}{

    The rank, which must be \code{>1} and \code{< nlevels}.

  }
  \item{levels}{

    Character representing the levels.

  }
  \item{lowerSQRT}{

    Logical. When \code{TRUE} a lower-triangular root
    \eqn{\mathbf{L}}{L} of the correlation or covariance matrix
    \eqn{\mathbf{C}}{C} is returned instead of the correlation
    matrix. Note that this matrix can have negative diagonal elements
    hence is not a (pivoted) Cholesky root.
    
  }
  
  \item{compGrad}{
    
    Logical. Should the gradient be computed? This is only possible for
    the C implementation.
    
  }
  \item{cov}{
    
    Integer \code{0}, \code{1} or \code{2}. If \code{cov} is \code{0},
    the matrix is a \emph{correlation} matrix (or its root). If
    \code{cov} is \code{1} or \code{2}, the matrix is a
    \emph{covariance} (or its root) with constant variance vector for
    \code{code = 1} and with arbitrary variance for \code{code = 2}. The
    variance parameters \code{par} are located at the tail of the
    \code{par} vector, so at locations \code{npCor + 1} to \code{npCor +
    nlevels} when \code{code = 2} where \code{npCor} is the number of
    correlation parameters.
  
  }
  \item{impl}{
    
    A character telling which of the C and R implementations should be
    chosen. The R implementation is only for checks and should not be
    used.

  }
  
} 
\value{
  
  A correlation matrix (or its root) with the optional \code{gradient}
  attribute.

}
\description{

  Compute the correlation matrix for a low-rank structure.

}
\details{

  The correlation matrix with size \eqn{m} is the general symmetric
  correlation matrix with rank \eqn{\leq r}{<= r} where \eqn{r} is
  given, as described by Rapisarda et al.  It depends on \eqn{(r - 1)
  \times (m - r / 2) / 2}{(r - 1) * (m - r / 2)} parameters
  \eqn{\theta_{ij}}{theta[i, j]} where the indices \eqn{i} and \eqn{j}
  are such that \eqn{1 \leq j < i}{1 <= j < i} for \eqn{i \leq r}{i <=
    r} or such that \eqn{1 \leq j < r}{1 <= j < r} for
  \eqn{r < i \leq n}{r < i <= m}. The parameters \eqn{\theta_{ij}}{theta[i, j]} are angles
  and are to be taken to be in \eqn{[0, 2\pi)}{[0, 2*pi)} if \eqn{j =
    1}{j = 1} and in \eqn{[0, \pi)}{[0, pi)} otherwise.

}
\note{

  This function is essentially for internal use and the corresponding
  correlation or covariance kernels are created as \code{covQual}
  objects by using the \code{\link{q1LowRank}} creator. 
  
  Here the parameters \eqn{\theta_{ij}}{theta[i, j]} are used \emph{in
  row order} rather than in the column order. This order simplifies the
  computation of the gradient.

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

  The \code{\link{q1LowRank}} creator of a correponding kernel object
  with class \code{"covQual"}, and the similar \code{\link{corLevSymm}}
  function for the full-rank case.

}
%%\author{Yves Deville}
%% \examples{XXX add examples here}