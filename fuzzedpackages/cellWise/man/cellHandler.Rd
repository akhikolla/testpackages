\name{cellHandler}
\alias{cellHandler}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cellHandler algorithm
}
\description{
This function flags cellwise outliers in \code{X} and imputes them, if robust estimates of the center \code{mu} and scatter matrix \code{Sigma} are given. When the latter are not known, as is typically the case, one can use the function \code{\link{DDC}} which only requires the data matrix \code{X}. Alternatively, the unknown center mu and scatter matrix Sigma can be estimated robustly from \code{X} by the function \code{\link{DI}}. 
}
\usage{
cellHandler(X, mu, Sigma, quant = 0.99)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{\code{X} is the input data, and must be an \eqn{n} by \eqn{d} matrix or a data frame.
  }
  \item{mu}{An estimate of the center of the data
  }
  \item{Sigma}{An estimate of the covariance matrix of the data
}
  \item{quant}{Cutoff used in the detection of cellwise outliers. Defaults to \code{0.99}
  }
  }
  
\value{
 A list with components: \cr
  \itemize{
    \item{\code{Ximp} \cr
    The imputed data matrix.
    }
      \item{\code{indcells} \cr
    Indices of the cells which were flagged in the analysis.
    }
    \item{\code{indNAs} \cr
    Indices of the NAs in the data. 
    }
    \item{\code{Zres} \cr
    Matrix with standardized cellwise residuals of the flagged cells. Contains zeroes in the unflagged cells.
    }
    \item{\code{Zres_denom} \cr
    Denominator of the standardized cellwise residuals.
    }
    \item{\code{cellPaths} \cr
    Matrix with the same dimensions as X, in which each row contains the path of least angle regression through the cells of that row, i.e. the order of the coordinates in the path (1=first, 2=second,...)
    }
    
    }
}
\references{
J. Raymaekers and P.J. Rousseeuw (2020). Handling cellwise outliers by sparse
regression and robust covariance. \emph{Arxiv: 1912.12446}. \href{https://arxiv.org/abs/1912.12446}{(link to open access pdf)}
}
\author{
J. Raymaekers and P.J. Rousseeuw
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DI}}
}
\examples{
mu <- rep(0, 3)
Sigma <- diag(3) * 0.1 + 0.9
X <- rbind(c(0.5, 1.0, 5.0), c(-3.0, 0.0, 1.0))
n <- nrow(X); d <- ncol(X)
out <- cellHandler(X, mu, Sigma)
Xres <- X - out$Ximp # unstandardized residual
mean(abs(as.vector(Xres - out$Zres*out$Zres_denom))) # 0
W <- matrix(rep(0,n*d),nrow=n) # weight matrix 
W[out$Zres != 0] <- 1 # 1 indicates cells that were flagged
# For more examples, we refer to the vignette:
vignette("DI_examples")
}