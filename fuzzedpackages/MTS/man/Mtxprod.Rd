\name{Mtxprod}
\alias{Mtxprod}
\title{Polynomial Matrix Product
}
\description{Compute the product of two polynomial matrices
}
\usage{
Mtxprod(Mtx, sMtx, p, P)
}
\arguments{
  \item{Mtx}{The coefficient matrix of a regular polynomial matrix
}
  \item{sMtx}{The coefficient matrix of a seasonal polynomial matrix
}
  \item{p}{Degree of the regular polynomial matrix
}
  \item{P}{Degree of the seasonal polynomial matrix
}
}
\value{Coefficient matrix of the product. The product is in the form 
reg-AR * sAR, etc.
}
\author{Ruey S. Tsay
}




