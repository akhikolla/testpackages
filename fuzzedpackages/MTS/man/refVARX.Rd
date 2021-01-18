\name{refVARX}
\alias{refVARX}
\title{Refining a VARX Model
}
\description{Refine a fitted VARX model by setting insignificant parameters 
to zero
}
\usage{
refVARX(m1, thres = 1)
}
\arguments{
  \item{m1}{An output object of the VARX command or the refVARX command
}
  \item{thres}{A threshold for the individual t-ratio. Default is 1.
}
}
\details{The program sets simultaneously all estimates with t-ratio less than 
threshold (in modulus) to zero and re-estimate the VARX model.
}
\value{The same as those of the command VARX.
}
\references{Tsay (2014, Chapter 6). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{VARX
}




