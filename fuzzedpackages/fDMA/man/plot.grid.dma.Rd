\name{plot.grid.dma}
\alias{plot.grid.dma}

\title{Plots Selected Outcomes from \code{\link{grid.DMA}} Function.}

\description{The function plots selected outcomes from \code{\link{grid.DMA}}.}

\usage{\method{plot}{grid.dma}(x,non.interactive=NULL, ...)}

\arguments{
\item{x}{an object of \code{grid.dma} class}
\item{non.interactive}{optional, \code{\link[base]{logical}}, indicate whether plots should be made in non-interactive mode, by default \code{non.interactive=FALSE}, i.e., the user specifies in the interactive menu which plots will be made}
\item{...}{not used}
}

\details{
If \code{x} comes from estimation of Dynamic Model Averaging (DMA), after executing the command, the User is asked to choose

1 - for plotting Root Mean Squared Error (RMSE) for all estimated models,

2 - for plotting Mean Absolute Error (MAE) for all estimated models,

3 - for plotting posterior inclusion probabilities (including constant) for all estimated models, the outcomes are saved in separate png files in the current working directory, and additionally, plots for different variables are collected into one big plot (also saved as a png file in the current working directory),

4 - for plotting expected coefficients (including constant) for all estimated models, the outcomes are saved in separate png files in the current working directory, and additionally, plots for different variables are collected into one big plot (also saved as a png file in the current working directory).

Chosing 0 exits the \code{plot} command.

---------------------------------------

If \code{x} comes from estimation of Dynamic Model Selection (DMS) or Median Probability Model (MED), after executing the command, the User is asked to choose

1 - for plotting Root Mean Squared Error (RMSE) for all estimated models,

2 - for plotting Mean Absolute Error (MAE) for all estimated models,

3 - for plotting expected coefficients (including constant) for all estimated models, the outcomes are saved in separate png files in the current working directory, and additionally, plots for different variables are collected into one big plot (also saved as a png file in the current working directory).

Chosing 0 exits the \code{plot} command.

If \code{non.interactive=TRUE} all the above plots are made.
}

\note{
It is suggested to execute \code{\link[grDevices]{graphics.off}} before exectuing \code{plot} command for \code{grid.dma} object. However, the User should take care to save all other plots before executing this command, as they can be lost.

If \code{\link[grDevices]{graphics.off}} is not executed before plotting \code{grid.dma} object, sometimes a legend might cover the important parts of the plot.

If any of the models comes from using multiple lambda (see \code{\link{fDMA}}), then RMSE and MAE are not plotted.

Also, if \code{\link[base]{length}(grid.alpha)} or \code{\link[base]{length}(grid.lambda)} is less than 2, then RMSE and MAE are not plotted.
}

\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
gra <- c(0.99,0.98,0.97)
grl <- c(0.99,0.95)

g1 <- grid.DMA(y=ld.wti,x=ld.drivers,grid.alpha=gra,grid.lambda=grl,initvar=1)
g2 <- grid.DMA(y=ld.wti,x=ld.drivers,grid.alpha=gra,grid.lambda=grl,initvar=1,model="dms")

# graphics.off()
plot(g1)
# graphics.off()
plot(g2)
}
}
