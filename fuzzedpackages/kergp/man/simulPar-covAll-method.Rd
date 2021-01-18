\docType{methods}
\name{simulPar,covAll-method}
\alias{simulPar,covAll-method}
\title{

  Draw Random Values for the Parameters of a Covariance Kernel

}
\usage{

\S4method{simulPar}{covAll}(object, nsim = 1L, seed = NULL)

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
}
\value{

  A matrix with \code{nsim} rows and \code{npar(object)} columns.

}
\description{

  Draw random values for the parameters of a covariance kernel

  object.

}
\details{

  Draw random values for the parameters of a covariance kernel
  object using the informations \code{coefLower} and
  \code{coefUpper}.

}
%% \author{ODY}
