\name{summary.tfCox}
\alias{summary.tfCox}
\title{
Summarize \code{tfCox} object
}
\description{
This function summarizes \code{tfCox} object
}
\usage{
\method{summary}{tfCox}(object, \dots)
}
\arguments{
  \item{object}{an object of class "tfCox"}
    \item{\dots}{
additional arguments to be passed. These are ignored in this function.
}
}

\details{
Summarize the fit by the number of knots and percent sparsity achieved. Percent sparsity is the percentage of features estimated to have no relationship with the outcome.
}

\author{
Jiacheng Wu
}

\seealso{
\code{\link{tfCox}}, \code{\link{plot.tfCox}}
}
\examples{
#generate data
set.seed(1234)
dat = sim_dat(n=100, zerof=0, scenario=1)

#fit piecewise constant for alpha=1 and a range of lambda
fit = tfCox(dat, ord=0, alpha=1)

#summarize the fit by the number of knots and percent sparsity achieved.
#Percent sparsity is the percentage of features estimated to have
#no relationship with outcome
summary(fit)
}
