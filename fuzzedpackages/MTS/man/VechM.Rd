\name{VechM}
\alias{VechM}
\title{Matrix constructed from output of the Vech Command. 
In other words, restore the original symmetric matrix from its 
half-stacking vector.
}
\description{Restores the symmetric matrix from the Vech command
}
\usage{
VechM(vec)
}
\arguments{
  \item{vec}{A vector representing the half-stacking of a symmetric matrix
}
}
\details{This command re-construct a symmetric matrix from output 
of the Vech command
}
\value{A symmetric matrix
}
\references{Tsay (2014, Appendix A)
}
\author{Ruey S. Tsay
}
\seealso{Vech
}
\examples{
v1=c(2,1,3)
m1=VechM(v1)
m1
}




