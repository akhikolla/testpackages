\name{simulPar}
\alias{simulPar}
\title{

  Generic function: Draw Random Values for the Parameters of a
  Covariance Kernel

}
\description{

 Generic function to draw random values for the parameters of a
 covariance kernel object.

}

\usage{
  simulPar(object, nsim = 1L, seed = NULL, ...)
}
\arguments{
  \item{object}{

    A covariance kernel.

  }
  \item{nsim}{

    Number of drawings.

  }
  \item{seed}{

    Seed for the random generator.

  }
  \item{...}{

    Other arguments for methods.

  }
}
\value{

  A matrix with \code{nsim} rows and \code{npar(object)} columns.

}
\details{
  
  Draw random values for the parameters of a covariance kernel object
  using the informations \code{coefLower} and \code{coefUpper}.

}
%% \author{}
