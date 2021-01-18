\name{check_irred}
\alias{check_irred}
\title{
Routine for Checking Irreducibility 
}

\description{
This routine checks the irreducibility of a set of edges, which are supposed to form the columns of a matrix. If a column is a positive linear combination of other columns, then it can be removed without affecting the problem; if there is a positive linear combination of columns of the matrix that equals the zero vector, then there is an implicit equality constraint in the matrix. In the former case, this routine delete the redundant columns and return a set of irreducible edges, while in the latter case, this routine will give the number of equality constraints in the matrix, and will leave this issue to the user to fix.
}

\usage{check_irred(mat)}

\arguments{
  \item{mat}{A matrix whose columns are edges.}
}

\value{
\item{edge}{The edges kept after being checked about irreducibility.}
\item{reducible}{A vector of the indice of the edges that are redundant in the original set of edges.}
\item{equal}{A vector showing the number of equality constraints in the original set of edges.}
}

\references{
Meyer, M. C. (1999) An extension of the mixed primal-dual bases algorithm to the case of more constraints than 
dimensions. \emph{Journal of Statistical Planning and Inference \bold{81}}, 13--31.

Meyer, M. C. (2013b) A simple new algorithm for quadratic programming with applications in statistics. \emph{Communications
in Statistics \bold{42(5)}}, 1126--1139.  

Liao, X. and M. C. Meyer (2014) coneproj: An R package for the primal or dual cone projections with routines for constrained regression. \emph{Journal of Statistical Software \bold{61(12)}}, 1--22.
}



\author{
Mary C. Meyer and Xiyue Liao  
}

\examples{
\dontrun{
  data(TwoDamat)
  dim(TwoDamat)
  ans <- check_irred(t(TwoDamat))
  }
}

