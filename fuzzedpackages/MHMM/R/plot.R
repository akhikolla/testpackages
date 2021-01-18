#' Plots of an instance of \code{\linkS4class{mhmmresults}}
#' 
#' @param x instance of  \code{\linkS4class{mhmmresults}}.
#' @param y numeric index of the subject to visualize.
#' @param col numeric indicates the latent state at each time (length must be equal to the length of x)
#' @param xlab character label of the x-axis
#' @param ylab character label of the y-axis
#' @param ylim numeric range of the y-axis
#'
#' @name plot
#' @rdname plot-methods
#' @docType methods
#' @exportMethod plot 
#' @aliases plot plot,mhmmresults-method plot,mhmmresults,numeric-method plot,mhmmresults,ANY-method 
#' @aliases plot plot,mhmmresults,numeric-method
#' @aliases plot plot,mhmmresults,ANY-method 
##' @examples
##' data(accelero)
##' # To make the estimation <5
##' res <- mhmm(accelero, K = 2, M = 4, nbcores = 1, nbinit = 5, iterSmall = 2)
##' plot(res, 1)
##' 
##'  \donttest{
##' data(accelero)
##' # It is better to increase the number of random initializations
##' res <- mhmm(accelero, K = 2, M = 4, nbcores = 1)
##' plot(res, 1)
##' }
setMethod(
  f="plot",
  signature = c("mhmmresults", "numeric"),
  definition = function(x, y, col=x@partitions$states[[y]], 
                        xlab="Time", ylab= "Activity", ylim=range(na.omit(x@data@yi[[y]]))){
    MHMM.plot(x, y, col = col, xlab = xlab, ylab = ylab, yrange = ylim)
})

ggcol <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

MHMM.plot <- function(x, y, col, xlab, ylab, yrange){
  cols <- ggcol(x@param@M)
  df <- data.frame(Time = 1:length(x@data@yi[[y]]))
  dfp <- data.frame(as.factor(df$Time), x@probabilities$states[[y]])
  colnames(dfp) <- c("Time", as.character(1:x@param@M))
  p1 <- ggplot(df, aes(x = 1:length(x@data@yi[[y]]),
                       y = x@data@yi[[y]],
                       color = factor(x@partitions$states[[y]], levels = as.character(1:x@param@M)))) +
    geom_point(pch = 21, size = 1.5, stroke = 1) +
    scale_x_continuous(xlab,
                       labels = df$Time[seq(1, nrow(dfp), length.out = 5)],
                       breaks = df$Time[seq(1, nrow(dfp), length.out = 5)],
                       expand = c(0.02,0)) +
    scale_y_continuous(ylab, limits = yrange) +
    scale_colour_manual(values=cols[1:x@param@M]) +
    theme(legend.position = "none",
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"))
  
  tabplot <- melt(dfp, id.vars="Time")
  p2 <- ggplot(tabplot, aes(x=tabplot$Time, tabplot$value, fill=tabplot$variable)) +
    geom_bar(stat="identity", position="fill", width=0.7) +
    scale_x_discrete(xlab,
                     labels = dfp$Time[seq(1, nrow(dfp), length.out = 5)],
                     breaks = dfp$Time[seq(1, nrow(dfp), length.out = 5)],
                     expand = c(0.02,0)) +
    scale_y_continuous("Probability",expand = c(0,0))+
    scale_colour_manual(values=cols[1:x@param@M])+
    scale_fill_discrete(name="State")+
    theme(legend.position="bottom", 
          legend.text=element_text(size=8),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          legend.key.size = unit(0.8,"line"),
          legend.margin = margin(0.01,0.01,0.01,0.01),
          plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"))
  return(grid.arrange(p1, p2, nrow=2, heights=c(0.6, 0.4)))
}


