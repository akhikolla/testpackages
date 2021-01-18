\name{refSCMfit}
\alias{refSCMfit}
\title{Refining Estimation of VARMA Model via SCM Approach
}
\description{Refine estimation of a VARMA model specified via the SCM 
approach by removing insignificant parameters
}
\usage{
refSCMfit(model, thres = 1)
}
\arguments{
  \item{model}{Name of the model from the SCMfit command or the refSCMfit command
}
  \item{thres}{Threshold for the t-ratio of individual coefficient. 
Default is 1.
}
}
\value{The same as those of the command SCMfit.
}
\references{Tsay (2014, Chapter 4)
}
\author{Ruey S. Tsay
}
\seealso{SCMfit
}




