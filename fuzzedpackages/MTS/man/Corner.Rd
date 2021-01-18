\name{Corner}
\alias{Corner}
\title{Compute the Corner table for transfer function model specification
}
\description{For a given dependent variable and an input variable, the program computes the Corner table for 
specifying the order (r,s,d) of a transfer function
}
\usage{
Corner(y,x,Nrow=11,Ncol=7)
}
\arguments{
  \item{y}{A pre-whitened dependent (or output) variable
}
  \item{x}{A pre-whitened independent (or input) variable. It should be a white noise series
}
  \item{Nrow}{The number of rows of the Corner table. Default is 11.
}
 \item{Ncol}{The number of columns of the Corner table. Default is 7.
 }
}
\details{For the pair of pre-whitened output and input variables, the program compute the 
Corner table and its simplified version for specifying the order of a transfer function.
}
\value{
\item{corner }{The Corner table}
}
\author{Ruey S. Tsay
}
