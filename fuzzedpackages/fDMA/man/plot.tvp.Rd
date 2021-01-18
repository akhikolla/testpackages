\name{plot.tvp}
\alias{plot.tvp}

\title{Plots Selected Outcomes from \code{tvp} Object.}

\description{The function plots selected outcomes from \code{tvp} object.}

\usage{\method{plot}{tvp}(x,non.interactive=NULL, ...)}

\arguments{
\item{x}{an object of \code{tvp} class}
\item{non.interactive}{optional, \code{\link[base]{logical}}, indicate whether plots should be made in non-interactive mode, by default \code{non.interactive=FALSE}, i.e., the user specifies in the interactive menu which plots will be made}
\item{...}{not used}
}

\details{
After executing the command, the User is asked to choose

1 - for plotting actual and predicted values,

2 - for plotting residuals,

3 - for plotting regression coefficients on one plot,

4 - for plotting regression coefficients in separate png files, saved in the current working directory, and moreover, to paste them into one big plot (also saved as a png file in the current working directory).

Chosing 0 exits the \code{plot} command.

If \code{non.interactive=TRUE} all the above plots are made.
}

\note{
It is suggested to execute \code{\link[grDevices]{graphics.off}} before exectuing \code{plot} command for \code{tvp} object. However, the User should take care to save all other plots before executing this command, as they can be lost.

If \code{\link[grDevices]{graphics.off}} is not executed before plotting \code{tvp} object, sometimes a legend might cover the important parts of the plot.
}

\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]

tvp <- tvp(y=ld.wti,x=ld.drivers,V=1,lambda=0.99)

# graphics.off()
plot(tvp)
}
}
