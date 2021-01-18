\name{plot3D.splmm}
\alias{plot3D.splmm}
\title{3D Plot the tuning results of a \code{'splmm.tuning'} object when tuning over both lambda 1 and lambda 2 grids}
\description{This function inputs an \code{'splmm.tuning'} object and plot the model selection criterion values in a 3D plot over the lambda 1 and lambda 2 tuning parameters grid.}
\usage{
\method{plot3D}{splmm}(x, criteria=c("BIC","AIC","BBIC","EBIC"),type=c("line","surface"),...)
}

\arguments{
  \item{x}{a \code{'splmm.tuning'} object with both \code{lam1.tuning=TRUE} and \code{lam2.tuning=TRUE}}
  \item{criteria}{A parameter specifying whether the criteria value the user want to plot is \code{BIC}, \code{AIC}, \code{BBIC} or \code{EBIC}. The default is \code{BIC}}
  \item{type}{A parameter specifying which type of 3D plot to use for plotting. Currently the available options include \code{line} plot and \code{surface} plot. The default is \code{surface} plot. }
  \item{\dots}{not used}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{plot3D}
}
\examples{

\donttest{
data(cognitive)

x <- model.matrix(ravens ~schoolid+treatment+year+sex+age_at_time0
                  +height+weight+head_circ+ses+mom_read+mom_write
                  +mom_edu, cognitive)
z <- x

## Tuning over lambda1 grid and lambda2 grid
lam1 = seq(0.1,0.5,0.1)
lam2 = seq(0.1,0.5,0.1)
fit1 <-splmmTuning(x=x,y=cognitive$ravens,z=z,grp=cognitive$id,lam1=lam1,
lam2=lam2,penalty.b="scad", penalty.L="scad")
plot3D(fit1)
}

}

