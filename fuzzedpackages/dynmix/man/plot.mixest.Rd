\name{plot.mixest}
\alias{plot.mixest}
\alias{plot}

\title{Plots Selected Outcomes from \code{mixest} Object.}

\description{The function plots selected outcomes from \code{mixest} object.}

\usage{\method{plot}{mixest}(x, ...)}

\arguments{
\item{x}{an object of \code{mixest} class}
\item{...}{not used}
}

\details{
The function plots a few outcomes from \code{mixest} object. First, the estimated regression coefficients are plotted separately for each variable. Credible intervals of 90\% are added. Next, if averaging was chosen for forecasting, then relative variable importances are plotted, i.e., sum of weights of models containing the given variable. If selection procedure was chosen for forecasting, it is plotted whether the given variable is included in the selected model at the given time. Finally weights from all component models are presented in one plot.
}

\examples{
\donttest{
data(oil) 
m1 <- mixest1(y=oil[,1,drop=FALSE],x=oil[,-1,drop=FALSE],ftype=2,V=100,W=100)
plot(m1)
}
}

\seealso{\code{\link[dynmix]{convts}}}
