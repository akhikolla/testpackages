\name{plot.etm}
\alias{plot.etm}

\title{Plot method for an etm object}
\description{
  Plot method for an object of class 'etm'. It draws the estimated
  transition probabilities in a basic scatterplot.
}
\usage{
\S3method{plot}{etm}(x, tr.choice, xlab = "Time",
     ylab = "Transition Probability", col = 1, lty, xlim, ylim,
     conf.int = FALSE, level = 0.95, ci.fun = "linear",
     ci.col = col, ci.lty = 3,
     legend = TRUE, legend.pos, curvlab, legend.bty = "n", ...)
}
\arguments{
  \item{x}{An object of class 'etm'}
  \item{tr.choice}{ character vector of the form 'c("from to","from
    to")' specifying which transitions should be plotted. Default, all
    the transition probabilities are plotted}
  \item{xlab}{x-axis label. Default is "Time"}
  \item{ylab}{y-axis label. Default is "Transition Probability"}
  \item{col}{Vector of colour. Default is black}
  \item{lty}{Vector of line type. Default is 1:number of transitions}
  \item{xlim}{Limits of x-axis for the plot}
  \item{ylim}{Limits of y-axis for the plot}
  \item{conf.int}{Logical. Whether to display pointwise confidence
    intervals. Default is FALSE.}
  \item{level}{Level of the conficence intervals. Default is 0.95.}
  \item{ci.fun}{Transformation applied to the confidence intervals. It
    could be different for all transition probabilities, though if
    \code{length(ci.fun) != number of transitions}, only \code{ci.fun[1]}
    will be used. Possible choices are "linear", "log", "log-log" and
    "cloglog". Default is "linear".}
  \item{ci.col}{Colour of the confidence intervals. Default is
    \code{col}.}
  \item{ci.lty}{Line type of the confidence intervals. Default is 3.}
  \item{legend}{A logical specifying if a legend should be added}
  \item{legend.pos}{A vector giving the legend's position. See
    \code{\link{legend}} for further details}
  \item{curvlab}{A character or expression vector to appear in the
    legend. Default is the name of the transitions}
  \item{legend.bty}{Box type for the legend}
  \item{\dots}{Further arguments for plot}
}
\details{
  By default, if the argument \code{strata} was used for creating the \code{etm}
  object, the first transition probability for all strata will be
  plotted. If there is no strata, all transition probabilities are
  plotted by default.

  In any case, a legend with be created by the labels are likely to be
  ugly. Please use the \code{curvlab} argument to control the test or
  use \code{legend = FALSE} and build your own legend. 
}
\value{
  No value returned
}
\author{Arthur Allignol, \email{arthur.allignol@gmail.com}}
\seealso{\code{\link{plot.default}}, \code{\link{legend}}, \code{\link{etm}} }
\examples{
data(sir.cont)

# Modification for patients entering and leaving a state
# at the same date
sir.cont <- sir.cont[order(sir.cont$id, sir.cont$time), ]
for (i in 2:nrow(sir.cont)) {
  if (sir.cont$id[i]==sir.cont$id[i-1]) {
    if (sir.cont$time[i]==sir.cont$time[i-1]) {
      sir.cont$time[i-1] <- sir.cont$time[i-1] - 0.5
    }
  }
}

tra <- matrix(ncol=3,nrow=3,FALSE)
tra[1, 2:3] <- TRUE
tra[2, c(1, 3)] <- TRUE

my.etm <- etm(sir.cont,c("0","1","2"),tra,"cens", s = 0)

plot(my.etm, tr.choice = c("0 0"))
}
\keyword{hplot}
