\name{plot.tvpreg}
\alias{plot.tvpreg}

\title{Plots Selected Outcomes from \code{tvpreg} Object.}

\description{The function plots selected outcomes from \code{tvpreg} object.}

\usage{\method{plot}{tvpreg}(x, ...)}

\arguments{
\item{x}{an object of \code{tvpreg} class}
\item{...}{not used}
}

\details{
The function plots the estimated regression coefficients, separately for each variable. 90\% credible intervals are added. 
}

\examples{
data(oil)
t1<- tvp.reg(y=oil[,1,drop=FALSE],x=oil[,-1,drop=FALSE],lambda=0.99,V=100,W=100)
plot(t1)
}

\seealso{\code{\link[dynmix]{convts}}}
