\name{plot.altf3}
\alias{plot.altf3}

\title{Plots Selected Outcomes from \code{altf3} Object.}

\description{The function plots selected outcomes from \code{altf3} object.}

\usage{\method{plot}{altf3}(x,non.interactive=NULL, ...)}

\arguments{
\item{x}{an object of \code{altf3} class}
\item{non.interactive}{optional, \code{\link[base]{logical}}, indicate whether plots should be made in non-interactive mode, by default \code{non.interactive=FALSE}, i.e., the user specifies in the interactive menu which plots will be made}
\item{...}{not used}
}

\details{
After executing the command, the User is asked to choose

1 - for plotting expected coefficients in separate png files, saved in the current working directory, and moreover, to paste them into one big plot (also saved as a png file in the current working directory),

2 - for plotting p-values (averaged over selected window sizes) for t-test of statistical significance for coefficients from a rolling regression, in separate png files, saved in the current working directory, and moreover, to paste them into one big plot (also saved as a png file in the current working directory),

3 - for plotting weights of all models used in averaging,

4 - for plotting expected window size.


Chosing 0 exits the \code{plot} command.

If \code{non.interactive=TRUE} all the above plots are made.
}

\note{
It is suggested to execute \code{\link[grDevices]{graphics.off}} before exectuing \code{plot} command for \code{altf3} object. However, the User should take care to save all other plots before executing this command, as they can be lost.

If \code{\link[grDevices]{graphics.off}} is not executed before plotting \code{altf3} object, sometimes a legend might cover the important parts of the plot.
}

\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
a <- altf3(y=wti,x=drivers,windows=c(36,100,150))

plot(a)
}
}
