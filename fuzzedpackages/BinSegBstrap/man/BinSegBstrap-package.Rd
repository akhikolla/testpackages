\name{BinSegBstrap-package}
\alias{BinSegBstrap-package}
\docType{package}

\title{Piecewise smooth regression by bootstrapped binary segmentation}

\description{
Provides methods for piecewise smooth regression. The main function \code{\link{BinSegBstrap}} estimates a piecewise smooth signal by applying a bootstrapped test recursively (binary segmentation approach). A single bootstrapped test for the hypothesis that the underlying signal is smooth versus the alternative that the underlying signal contains at least one change-point can be performed by the function \code{\link{BstrapTest}}. A single change-point is estimated by the function \code{\link{estimateSingleCp}}. More details can be found in the vignette. Parts of this work were inspired by Gijbels and Goderniaux (2004).}


\section{Acknowledgement}{
This work results from a summer research project at the University of Cambridge in 2019. Kate McDaid was supported by a bursary from the summer research programme of the Centre of Mathematics at the University of Cambridge. Florian Pein's position is funded by the EPSRC programme grant 'StatScale: Statistical Scalability for Streaming Data'.
}


\encoding{utf8}
\references{
Gijbels, I., Goderniaux, A-C. (2004) Bootstrap test for change-points in nonparametric regression. \emph{Journal of Nonparametric Statistics} \bold{16}(3-4), 591--611.
}

\keyword{package}
\keyword{nonparametric}

\seealso{\code{\link{BinSegBstrap}}, \code{\link{BstrapTest}}, \code{\link{estimateSingleCp}}}

\examples{
n <- 200
signal <- sin(2 * pi * 1:n / n)
signal[51:100] <- signal[51:100] + 5
signal[151:200] <- signal[151:200] + 5

y <- rnorm(n) + signal

est <- BinSegBstrap(y = y)

plot(y)
lines(signal)
lines(est$est, col = "red")


n <- 100
signal <- sin(2 * pi * 1:n / n)
signal[51:100] <- signal[51:100] + 5

y <- rnorm(n) + signal

test <- BstrapTest(y = y)
est <- estimateSingleCp(y = y)

plot(y)
lines(signal)
lines(est$est, col = "red")
}
