\name{Mtxprod1}
\alias{Mtxprod1}
\title{Alternative Polynomial Matrix Product
}
\description{Compute the product of two polynomial matrices
}
\usage{
Mtxprod1(Mtx, sMtx, p, P)
}
\arguments{
  \item{Mtx}{The coefficient matrix of a regular polynomial matrix
}
  \item{sMtx}{The coefficient matrix of a seasonal polynomial matrix
}
  \item{p}{Degree of the regular polynomial matrix. p is less than P.
}
  \item{P}{Degree of the seasonal polynomial matrix
}
}

\details{This polynomial product is used in seasonal VARMA modeling to check the multiplicative nature between the regular and seasonal polynomial matrices
}

\value{Coefficient matrix of the product. The product matrix is in the form 
sAR * reg-AR, etc.
}
\author{Ruey S. Tsay
}




