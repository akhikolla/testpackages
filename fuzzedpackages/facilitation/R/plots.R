#' Plots 
#' 
#' Plotting functions.
#' 
#' The \code{stackplot} function produces a stacked plot of the population over time.
#' Notice that the population should have at least two stages for this function to work.
#' 
#' @param mat A population matrix, as produced by \code{\link{abundance.matrix}} or something that
#' can be coerced to matrix
#' @param col Optional. A color vector
#' @param legend Optional. An array of names
#' @param log.y Logical. Should the y-axis be plotted in a logarithmic scale?
#' @param perc Logical. If set to true, will output the y-axis as a percentage instead of the absolute numbers
#' @param qt Optional. For distributions, show only up to quantile qt (percentage)
#' @param \dots Further parameters to be passed to the lower level plot function
#' @examples
#' data(twospecies)
#' ab <- abundance.matrix(twospecies,seq(0,twospecies$maxtime,by=1))
#' # species 1
#' stackplot(ab[,1:3])
#' # species 2
#' stackplot(ab[,4:5])
#' @export
#' @import grDevices 
#' @import graphics
stackplot <- function(mat, col, legend, log.y = FALSE, perc=F, qt=100, ...) {
	dots <- list(...)
	if(missing(col))
		#col <- terrain.colors(dim(mat)[2])
		col <- colorRampPalette(c("darkred","pink"))(dim(mat)[2])
	if (log.y) {
		minp <- 1
		log <- "y"
	} else {
		minp <- 0
		log <- ""
	}

    mat<-as.matrix(mat)

	N <- dim(mat)[2]
	time <- as.numeric(rownames(mat))
    if(N>1){
        for (i in (N-1):1) # sums populations IN REVERSE
            mat[,i] = mat[,i] + mat[,i+1]
    }
	mat <- cbind(mat, rep(minp, length(time)))
	# maximo da escala do grafico
	maxp <-max(mat[,1])

    # percentage
    if(perc){
        mat <- mat*100.0/maxp
        maxp <- 100
        minp <- 100.0*minp/maxp
        if (! "ylab" %in% names(dots)) dots$ylab = "Population (%)"
    }

    # cap at quantile
    if(qt<100){
        quant <- maxp*(100.0-qt)/100.0
        linemax <- max(which(mat[,1]>=quant))
        mat <- mat[1:linemax,]
        time <- time[1:linemax]
    }

	if (! "ylim" %in% names(dots)) dots$ylim = c(minp, maxp)
	if (! "xlim" %in% names(dots)) dots$xlim = c(min(time),max(time))
	if (! "main" %in% names(dots)) dots$main = "Population dynamics"
	if (! "ylab" %in% names(dots)) dots$ylab = "Population"
	if (! "xlab" %in% names(dots)) dots$xlab = "Time"

	do.call(plot, c(list(1, type='n', log=log), dots))
	x <- c(time, rev(time))
	for (i in 1:(N)) {
		y <- c(mat[,i], rev(mat[,i+1]))
		polygon(x,y, col=col[i])
	}
    if(N>1){ # legend is unnecessary if N==1
        if (missing(legend)) { 
            if(N == 2) legend <- c("Juveniles", "Adults")
            if(N == 3) legend <- c("Seeds", "Juveniles", "Adults")
            if(N > 3) legend <- c(1:N)
        }
        legend("topleft", legend=legend, fill=col, bg="white")
    }
}

#' Function for ploting simulation as a gif
#'
#' The spatialanimation function plots the individuals of the selected stages over time. Use plotsnapshot
#' for plotting a single instant.
#'
#' @author Alexandre Adalardo de Oliveira - 16/03/2016
#' @author M. Salles
#' @param data	result of a simulation, created by \code{\link{community}}
#' @param times	array of times at which to plot
#' @param interval a time length to wait between frames
#' @param draw an array of stages id, to be drawn bottom to top. Absent stages will not be
#' drawn.
#' @param xlim	Optional. Limits to the x-axis
#' @param ylim	Optional. Limits to the y-axis
#' @param color 	Optional. A color vector
#' @param radius Optional. Array representing the sizes in which the individuals will be drawn. 
#' Defaults to interaction radius.
#' @param movie.name The filename of the gif that will be saved.
#' @examples
#' data(twospecies)
#' spatialanimation(twospecies,draw=c(5,3),times=seq(0,10,1),movie.name="ts.gif")
#' @export
#' @import grDevices 
#' @import graphics
#' @import animation
#' @import grid
spatialanimation = function(data, times=seq(0,data$maxtime,length.out=50), interval=0.1,
                            draw=data$num.total:1,
                            radius=data$param$radius[draw],
                            color=colorRampPalette(c("darkred","lightgreen"))(length(draw)),
                            movie.name="facilitationmovie.gif",
                            xlim=c(0,data$w), ylim=c(0,data$h)
                            )
{
    # creates list of dataframes, one for each time
    d<-data$data
	dtlist <- lapply(times,function(t){subset(d,d$begintime <= t & (d$endtime >= t | is.na(d$endtime)),select=c(1,3,4))})
    maxst <- data$num.total
    # set minimum radius for stages with rad=0
    for(i in 1:length(radius)) if(radius[i] == 0) radius[i] = 0.05
    saveGIF(spatialplot(dtlist,times=times,xlim=xlim,ylim=ylim,sp=draw,color,radius),interval=interval,movie.name=movie.name)
}


# function for ploting simulation frames
#
# @author Alexandre Adalardo de Oliveira - 16/03/2016
# @author M. Salles
# @param dtlist	A list of data.frames, one for each time in times, containing columns sp,x,y 
# @param times	Array of times, corresponding to the data.frames
# @param xlim	Limits to the x-axis
# @param ylim	Limits to the y-axis
# @param sp    Array of species/stages id, in order of plotting bottom to top
# @param col 	A color array, one for each species/stage
# @param radius Array or radiuses, one for each species/stage
spatialplot = function(dtlist, times, xlim, ylim, sp, 
                       col, radius)
{
    n <- length(sp)
    # init viewport
    vp <- viewport(width = 0.8, height = 0.8, xscale=xlim, yscale=ylim)
    # loop through times
    for (i in 1:length(times))
    {
        dt = dtlist[[i]]
        if(dim(dt)[1] > 0) {# interrupt if population is zero

            grid.newpage()
            pushViewport(vp)
            grid.rect(gp = gpar(col = "gray"))
            for (j in 1:n){
                dtsp <- dt[dt$sp==sp[j],]
                if(dim(dtsp)[1] > 0){
                    grid.circle(x = dtsp$x, y=dtsp$y, r=radius[j],default.units="native", gp=gpar(fill=col[j],col=col[j]))
                }
            }
            grid.text(paste("t =",round(times[i],digits=4)), y=1.06)
            grid.xaxis(at=round(seq(xlim[1],xlim[2], len=5)))
            grid.yaxis(at=round(seq(ylim[1],ylim[2], len=5)))
        }
    }
}

#' @export
#' @param t a single time at which to plot
#' @param \dots additional parameters to be passed to spatialanimation
#' @rdname spatialanimation
#' @examples
#' data(twospecies)
#' plotsnapshot(twospecies,t=10)
plotsnapshot <- function(data,t,...) {
    spatialanimation(data,c(t,t),...)
}
