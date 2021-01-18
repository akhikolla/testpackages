#' @title plot the output from transformPhylo.MCMC
#' @description Plots a histogram of the estimated parameter and a trace of the results
#' @param mcmc.input an object of class "mcmc.motmot" output from \code{\link{transformPhylo.MCMC}}
#' @param y.limit the limits for the y axes for the plots
#' @param x.limit the limits for the x axes for the plots
#' @param label.text the labels for the two plots defaults to '(a)' etc., for the histogram and '(b)' etc., for the trace plot
#' @param cex.axis character expansion for the plot axis labels
#' @param cex.labels character expansion for the plot axis names
#' @param col.hist colour for the histogram bars
#' @param col.trace colour for the trace plot
#' @return Two plots showing the histogram of the estimated parameter value and a trace of the MCMC estimation
#' @author Mark Puttick
#' @examples
#' library(motmot)
#' data(anolis.tree)
#' data(anolis.data)
#' attach(anolis.data)
#' male.length <- matrix(Male_SVL, dimnames=list(rownames(anolis.data)))
#' sortedData <- sortTraitData(anolis.tree, male.length)
#' phy <- sortedData$phy
#' male.length <- sortedData$trait
#' phy.clade <- extract.clade(phy, 182)
#' male.length.clade <- as.matrix(male.length[match(phy.clade$tip.label, rownames(male.length)),])
#' ## not run
#' # please note, this model will be need to run for longer to achieve convergence
#' # lambda.mcmc <- transformPhylo.MCMC(y=male.length.clade, phy=phy.clade, 
#' # model="lambda", mcmc.iteration=100, burn.in=0.1)
#' # mcmc.plot(lambda.mcmc)
#' @export


mcmc.plot <- function(mcmc.input, y.limit=NULL, x.limit=NULL, label.text=NULL, cex.axis=1, cex.labels=0.7, col.hist="green4", col.trace="navy") {
	if(!class(mcmc.input) == "motmot.mcmc") stop("please supply object of class motmot.mcmc")
	n.param <- ncol(mcmc.input$mcmc)
	names.param <- names(mcmc.input[[1]])
	
	hist.l <- lapply(1:n.param, function(x) suppressWarnings(hist(mcmc.input$mcmc[,x], xlim=x.limit, ylim=y.limit, xaxs="i", yaxs="i", las=1, ylab="density", xlab=names.param[x], plot=FALSE)))
	
	if(is.null(y.limit)) y.limit <- sapply(hist.l, function(x) c(0, max(x$counts) * 1.05))
	if(is.null(x.limit)) x.limit <- sapply(hist.l, function(x) {
		if(x$breaks[1] > 0) {
			diffs <- diff(range(x$breaks)) / 20
			out <- c(0, max(x$breaks) + diffs)
			} else {
			diffs <- diff(range(x$breaks)) / 20
			out <- c(min(x$breaks) - diffs, max(x$breaks) + diffs)
			}
			out
		}
	)
	
	if(is.null(label.text))  label.text <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)")
	par(mfrow=c(n.param,2))
	count <- 1
	for(u in 1:n.param) {
		plot(hist.l[[u]], xaxs="i", yaxs="i", las=1, main="", xlab=names.param[u], col=col.hist, border="white", xlim=x.limit[,u], ylim=y.limit[,u], cex.axis=cex.axis)
		mtext(label.text[count], 3, at=0, cex=cex.labels, font=2)
		count <- count + 1
		box()
		plot(mcmc.input$mcmc[,u], ylim=x.limit[,u], xaxs="i", yaxs="i", las=1, ylab=paste0(names.param[u]), xlab="generations", type="l", cex.axis=cex.axis, col=col.trace)
		mtext(label.text[count], 3, at=0, cex=cex.labels, font=2)
		count <- count + 1
		}
}






