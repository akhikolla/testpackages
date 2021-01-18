\name{summary.reg}
\alias{summary.reg}

\title{Summarizes Outcomes from \code{reg} Object.}

\description{The function summarizes selected outcomes obtained from \code{\link{roll.reg}} and \code{\link{rec.reg}}.}

\usage{\method{summary}{reg}(object, ...)}

\arguments{
\item{object}{an object of \code{reg} class}
\item{...}{not used}
}

\details{
The function produces the outcomes as \code{\link{print.reg}}.

Additionally, it provides how often p-values for t-test of statistical significance for each independent variable in the model is below 1\%, 5\% and 10\%, respectively.
}

\examples{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
roll <- roll.reg(y=ld.wti,x=ld.drivers,window=100)
rec <- rec.reg(y=ld.wti,x=ld.drivers)
summary(roll)
summary(rec)
}
