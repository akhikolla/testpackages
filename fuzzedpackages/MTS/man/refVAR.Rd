\name{refVAR}
\alias{refVAR}
\title{Refining a VAR Model
}
\description{Refine a fitted VAR model by removing simultaneously insignificant parameters
}
\usage{
refVAR(model, fixed = NULL, thres = 1)
}
\arguments{
  \item{model}{An output object of the command VAR or the refVAR command
}
  \item{fixed}{A logical matrix for VAR polynomial structure
}
  \item{thres}{Threshold used to set parameter to zero. Default is 1.
}
}
\details{Refine a VAR fitting by setting all estimates with t-ratio 
less than the threshold (in modulus) to zero.
}
\value{The same as those of the command VAR
}
\references{Tsay (2014, Chapter 2)
}
\author{Ruey S. Tsay
}
\seealso{VAR
}
\examples{
data("mts-examples",package="MTS")
gdp=log(qgdp[,3:5])
zt=diffM(gdp)
m1=VAR(zt,3)
m2=refVAR(m1,thres=1.0)
names(m2)
}




