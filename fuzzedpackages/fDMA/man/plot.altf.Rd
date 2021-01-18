\name{plot.altf}
\alias{plot.altf}

\title{Plots Selected Outcomes from \code{altf} Object.}

\description{The function plots selected outcomes from \code{altf} object.}

\usage{\method{plot}{altf}(x,non.interactive=NULL, ...)}

\arguments{
\item{x}{an object of \code{altf} class}
\item{non.interactive}{optional, \code{\link[base]{logical}}, indicate whether plots should be made in non-interactive mode, by default \code{non.interactive=FALSE}, i.e., the user specifies in the interactive menu which plots will be made}
\item{...}{not used}
}

\details{
After executing the command, the User is asked to choose

1 - for plotting regression coefficients in separate png files, saved in the current working directory, and moreover, to paste them into one big plot (also saved as a png file in the current working directory),

2 - for plotting p-values for t-test of statistical significance for regression coefficients from applied models, in separate png files, saved in the current working directory, and moreover, to paste them into one big plot (also saved as a png file in the current working directory).

Chosing 0 exits the \code{plot} command.

If \code{non.interactive=TRUE} all the above plots are made.
}

\note{
Coefficients are plotted only for rec. OLS, roll. OLS, TVP, TVP-AR(1) and TVP-AR(2) models. P-values -- for rec. OLS and roll. OLS. 

It is suggested to execute \code{\link[grDevices]{graphics.off}} before exectuing \code{plot} command for \code{altf} object. However, the User should take care to save all other plots before executing this command, as they can be lost.

If \code{\link[grDevices]{graphics.off}} is not executed before plotting \code{altf} object, sometimes a legend might cover the important parts of the plot.
}

\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
a <- altf(y=wti,x=drivers)

plot(a)
}
}
