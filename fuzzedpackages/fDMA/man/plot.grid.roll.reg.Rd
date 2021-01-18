\name{plot.grid.roll.reg}
\alias{plot.grid.roll.reg}

\title{Plots Selected Outcomes from \code{\link{grid.roll.reg}} Function.}

\description{The function plots selected outcomes from \code{\link{grid.roll.reg}}.}

\usage{\method{plot}{grid.roll.reg}(x,non.interactive=NULL, ...)}

\arguments{
\item{x}{an object of \code{grid.roll.reg} class}
\item{non.interactive}{optional, \code{\link[base]{logical}}, indicate whether plots should be made in non-interactive mode, by default \code{non.interactive=FALSE}, i.e., the user specifies in the interactive menu which plots will be made}
\item{...}{not used}
}

\details{
After executing the command, the User is asked to choose

1 - for plotting Root Mean Squared Error (RMSE) for all estimated models,

2 - for plotting Mean Absolute Error (MAE) for all estimated models,

3 - for plotting coefficients (including constant) for all estimated models, the outcomes are saved in separate png files in the current working directory, and additionally, plots for different variables are collected into one big plot (also saved as a png file in the current working directory),

4 - for plotting p-values for t-test of statistical significance for regression coefficients for all estimated models, the outcomes are saved in separate png files in the current working directory, and additionally, plots for different variables are collected into one big plot (also saved as a png file in the current working directory),

Chosing 0 exits the \code{plot} command.

If \code{non.interactive=TRUE} all the above plots are made.
}

\note{
It is suggested to execute \code{\link[grDevices]{graphics.off}} before exectuing \code{plot} command for \code{grid.roll.reg} object. However, the User should take care to save all other plots before executing this command, as they can be lost.

If \code{\link[grDevices]{graphics.off}} is not executed before plotting \code{grid.roll.reg} object, sometimes a legend might cover the important parts of the plot.
}

\examples{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
\dontrun{
grw <- c(50,100,150)
g <- grid.roll.reg(y=ld.wti,x=ld.drivers,grid.window=grw)

plot(g)
}
}
