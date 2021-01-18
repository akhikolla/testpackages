\name{summary.starma}
\alias{summary.starma}
\alias{summary.starma}
\alias{print.summary.starma}
\title{Summary method for space-time series fitted models}
\description{
\code{summary} method for class "\code{\link{starma}}".
}
\usage{
\method{summary}{starma}(object, ...)
\method{print}{summary.starma}(x, ...)
}
\arguments{
\item{object}{ a \code{starma} class object. }
\item{x}{a \code{summary.starma} class object. }
\item{...}{ unused }
}
\details{
\code{print.summary.starma} formats the coefficients, standard errors, etc. and additionally gives 'significance stars'.
}
\value{
An object of class \code{summary.starma} containing the following elements:
\item{call}{ An object of mode "\code{call}": a symbolic description of the fitted model }
\item{coefficients}{ A data frame containing the estimates, standard errors, etc. of the coefficients of the fitted model }
}
\author{Felix Cheysson}
\examples{
data(nb_mat)	# Get neighbourhood matrices

# Simulate a STARMA model
eps <- matrix(rnorm(94*200), 200, 94)
sim <- eps
for (t in 3:200) {
	sim[t,] <- (.4*diag(94) + .25*blist[[2]]) \%*\% sim[t-1,] +
		(.25*diag(94)                ) \%*\% sim[t-2,] +
		(            - .3*blist[[2]]) \%*\% eps[t-1,] +
		eps[t, ]
}

sim <- sim[101:200,]
sim <- stcenter(sim)	# Center and scale the dataset

# Select parameters to estimate
ar <- matrix(0, 2, 2)
ar[ ,1] <- 1	# phi10 and phi20
ar[1,2] <- 1	# phi11
ma <- matrix(c(0,1), 1, 2)	# theta11

# Run the Kalman filter algorithm
model <- starma(sim, blist, ar, ma)

# Get summary
summary(model)
}
\keyword{starma}
\keyword{summary}
\keyword{summary.starma}