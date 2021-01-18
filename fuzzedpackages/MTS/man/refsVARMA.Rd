\name{refsVARMA}
\alias{refsVARMA}
\title{Refining a Seasonal VARMA Model
}
\description{Refines a fitted seasonal VARMA model by setting insignificant 
estimates to zero
}
\usage{
refsVARMA(model, thres = 0.8)
}
\arguments{
  \item{model}{An output object of the sVARMA command or the refsVARMA 
command
}
  \item{thres}{Threshold for individual t-ratio. Default is 0.8. 
}
}
\details{The command removes simultaneously all parameters with t-ratio 
less than the threshold in modulus.
}
\value{The same as those of the command sVARMA
}
\references{Tsay (2014, Chapter 6)
}
\author{Ruey S. Tsay
}
\seealso{sVARMA
}




