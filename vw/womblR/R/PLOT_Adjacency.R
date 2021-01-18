###Function for plotting posterior adjacencies Wij(alpha_t)
#'
#' PlotAdjacency
#'
#' Plots a heat map of the differential light sensitivity on the Humphrey Field
#' Analyzer-II visual field.
#'
#' @param Wij a \code{\link{PosteriorAdj}} object.
#'
#' @param Visit either an integer \code{(1,...,Nu)} indicating the visit number for which
#' you want to get the adjacencies to plot or NA. If NA, then the plot will produce the
#' dissimilarity metric at each adjacency.
#'
#' @param stat either "mean" or "sd" (only used for Visit != NA).
#'
#' @param main an overall title for the plot.
#'
#' @param color.scheme a vector of colors to be used to show the adjacencies changing.
#'
#' @param edgewidth a scalar indicating the width of the edges.
#'
#' @param cornerwidth a scalar indicating the width of the corners.
#'
#' @param lwd.border a scalar indicating width of the visual field border.
#'
#' @param color.bs one color specifying the blind spot.
#'
#' @param zlim the limits used for the legend (default are c(0,1)).
#'
#' @param legend logical, indicating whether the legend should be present (default = TRUE).
#'
#' @param DM a dissimilarity metric to be plotted at each location on the visual field (default = NULL).
#'
#' @param W an adjacency matrix that specifies the visual field, required if Wij is not provided (default = NULL).
#'
#' @details \code{PlotAdjacency} is used in the application of glaucoma progression to
#'  plot the posterior mean and standard deviation neighborhood adjacencies across the
#'  visual field.
#'
#' @examples
#' ###Define blind spot locations on the HFA-II
#' blind_spot <- c(26, 35)
#'
#' ###Load visual field adjacency matrix
#' W <- HFAII_Queen[ -blind_spot, -blind_spot]
#'
#' ###Load Garway-Heath angles for dissimiliarity metric
#' DM <- GarwayHeath[-blind_spot] #Uses Garway-Heath angles object "GarwayHeath"
#'
#' ###Adjacency plots
#' PlotAdjacency(W = W, DM = DM, zlim = c(0, 180), Visit = NA,
#'               main = "Garway-Heath dissimilarity metric\n across the visual field")
#'
#' @author Samuel I. Berchuck
#'
#' @export
PlotAdjacency <- function(Wij,
                          Visit = 1,
                          stat = "mean",
                          main = "Estimated Adjacencies",
                          color.scheme = c("Black","White"),
                          edgewidth = 2,
                          cornerwidth = 1 / 4,
                          lwd.border = 3,
                          color.bs = "gray",
                          zlim = c(0, 1),
                          legend = TRUE,
                          DM = NULL,
                          W = NULL) {

  ##Note: Depends on library classInt
  # You need the suggested package for this function
  if (!requireNamespace("classInt", quietly = TRUE)) {
    stop("classInt needed for this function to work. Please install it.",
          call. = FALSE)
  }

  ###Logical function to check for colors
  areColors <- function(x) {
    sapply(x, function(X) {
      tryCatch(is.matrix(col2rgb(X)),
               error = function(e) FALSE)
    })
  }

  ###Check inputs
  if (missing(Wij)) {

    if (is.null(W)) stop('"W" is required when "Wij" is missing')
    if (!is.na(Visit)) stop('"Visit" must be NA when Wij is missing')
    if (!any(areColors(color.scheme))) stop('"color.scheme" can only include colors')

    ###Border information
    AdjacentEdgesBoolean <- (W == 1) & (!lower.tri(W))
    Dist <- function(x, y) pmin(abs(x - y), (360 - pmax(x, y) + pmin(x, y))) #arc length of optic nerve
    DM_Grid <- expand.grid(DM, DM)
    Z_Vector <- Dist(DM_Grid[ , 1], DM_Grid[ , 2])
    Z_Matrix <- matrix(Z_Vector, nrow = dim(W)[1], ncol = dim(W)[1], byrow = TRUE)
    Z <- matrix(Z_Matrix[AdjacentEdgesBoolean], ncol = 1)
    Boundary <- cbind(which(AdjacentEdgesBoolean, arr.ind = TRUE), Z)
    Boundary[ ,3] <- 180 - Boundary[ ,3]
    NAdjacency <- dim(Boundary)[1]
    Degree <- DM

    ###Set color pallete
    col.breaks <- seq(zlim[1],zlim[2],length.out=(length(unique(Boundary[,3]))+1))
    col.br <- colorRampPalette(color.scheme)
    col.pal <- col.br(NAdjacency)
    # suppressWarnings(fixed_obs<-classIntervals(Boundary[,3],n=length(unique(Boundary[,3]))))
    suppressWarnings(fixed_obs<-classInt::classIntervals(Boundary[,3],style="fixed",fixedBreaks=col.breaks))
    color.adj<-classInt::findColours(fixed_obs,col.pal)

    ###Plotting functions and Parameters
    lwd<-edgewidth
    lend<-2
    format2<-function(x) format(round(x,2),nsmall=2)
    format0<-function(x) format(round(x,0),nsmall=0)
    tri<-cornerwidth
    point<-function(x,y,ulbr=TRUE, col=colorAdj) {
      if (!ulbr) {
        polygon(c(x,x+tri,x),c(y,y,y+tri),col=col,border=col)
        polygon(c(x,x-tri,x),c(y,y,y-tri),col=col,border=col)
      }
      if (ulbr) {
        polygon(c(x,x-tri,x),c(y,y,y+tri),col=col,border=col)
        polygon(c(x,x+tri,x),c(y,y,y-tri),col=col,border=col)
      }
    }

    ###Create edge indicator
    EdgesInd<-numeric(length=NAdjacency)
    WhichEdges<-c(1,2,3,5,7,9,11,13,15,17,18,20,22,24,26,28,30,32,34,36,38,40,42,43,45,47,49,51,53,55,57,59,61,63,65,67, 69,71,72,75,77,79,81,83,85,87,89,91,93,95,96,97,99,102,104,106,108,110,112,114,116,118,119,122,123,124,126,129,131,133,135,137,139,141,143,145,147,149,152,154,156,158,160,162)
    EdgesInd[WhichEdges]<-rep(1,length(WhichEdges))
    WhichCorners<-which(EdgesInd==0)

    ###Create plot
    pardefault <- suppressWarnings(par(no.readonly = T))
    par(mfcol = c(1, 1), pty = "m", mai = c(0, 0, 0.75, 0))
    plot(1,1,xlim=c(1,13),ylim=c(1.25,8.75),type="n",xaxt="n",yaxt="n",bty="n",ylab="",xlab="",main=main, asp = 1, cex.main = 1.5)

    ###Plot edges
    for (i in WhichEdges) {
      Indeces<-as.numeric(Boundary[i,1:2])
      colorAdj<-color.adj[i]

      if (identical(Indeces,c(1,2))) segments(5,8,5,9,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(2,3))) segments(6,8,6,9,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(3,4))) segments(7,8,7,9,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(1,5))) point(4,8,ulbr=FALSE, colorAdj)
      if (identical(Indeces,c(1,6))) segments(4,8,5,8,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(2,6))) point(5,8,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(5,6))) segments(4,7,4,8,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(1,7))) point(5,8,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(2,7))) segments(5,8,6,8,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(3,7))) point(6,8,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(6,7))) segments(5,7,5,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(2,8))) point(6,8,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(3,8))) segments(6,8,7,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(4,8))) point(7,8,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(7,8))) segments(6,7,6,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(3,9))) point(7,8,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(4,9))) segments(7,8,8,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(8,9))) segments(7,7,7,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(4,10))) point(8,8,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(9,10))) segments(8,7,8,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(5,11))) point(3,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(5,12))) segments(3,7,4,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(6,12))) point(4,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(11,12))) segments(3,6,3,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(5,13))) point(4,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(6,13))) segments(4,7,5,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(7,13))) point(5,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(12,13))) segments(4,6,4,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(6,14))) point(5,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(7,14))) segments(5,7,6,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(8,14))) point(6,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(13,14))) segments(5,6,5,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(7,15))) point(6,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(8,15))) segments(6,7,7,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(9,15))) point(7,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(14,15))) segments(6,6,6,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(8,16))) point(7,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(9,16))) segments(7,7,8,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(10,16))) point(8,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(15,16))) segments(7,6,7,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(9,17))) point(8,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(10,17))) segments(8,7,9,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(16,17))) segments(8,6,8,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(10,18))) point(9,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(17,18))) segments(9,6,9,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(11,19))) point(2,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(11,20))) segments(2,6,3,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(12,20))) point(3,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(19,20))) segments(2,5,2,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(11,21))) point(3,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(12,21))) segments(3,6,4,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(13,21))) point(4,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(20,21))) segments(3,5,3,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(12,22))) point(4,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(13,22))) segments(4,6,5,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(14,22))) point(5,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(21,22))) segments(4,5,4,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(13,23))) point(5,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(14,23))) segments(5,6,6,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(15,23))) point(6,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(22,23))) segments(5,5,5,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(14,24))) point(6,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(15,24))) segments(6,6,7,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(16,24))) point(7,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(23,24))) segments(6,5,6,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(15,25))) point(7,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(16,25))) segments(7,6,8,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(17,25))) point(8,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(24,25))) segments(7,5,7,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(17,26))) point(9,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(18,26))) segments(9,6,10,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(19,27))) segments(1,5,2,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(20,27))) point(2,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(19,28))) point(2,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(20,28))) segments(2,5,3,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(21,28))) point(3,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(27,28))) segments(2,4,2,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(20,29))) point(3,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(21,29))) segments(3,5,4,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(22,29))) point(4,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(28,29))) segments(3,4,3,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(21,30))) point(4,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(22,30))) segments(4,5,5,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(23,30))) point(5,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(29,30))) segments(4,4,4,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(22,31))) point(5,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(23,31))) segments(5,5,6,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(24,31))) point(6,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(30,31))) segments(5,4,5,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(23,32))) point(6,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(24,32))) segments(6,5,7,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(25,32))) point(7,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(31,32))) segments(6,4,6,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(24,33))) point(7,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(25,33))) segments(7,5,8,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(32,33))) segments(7,4,7,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(26,34))) segments(9,5,10,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(27,35))) point(2,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(28,35))) segments(2,4,3,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(29,35))) point(3,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(28,36))) point(3,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(29,36))) segments(3,4,4,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(30,36))) point(4,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(35,36))) segments(3,3,3,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(29,37))) point(4,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(30,37))) segments(4,4,5,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(31,37))) point(5,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(36,37))) segments(4,3,4,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(30,38))) point(5,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(31,38))) segments(5,4,6,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(32,38))) point(6,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(37,38))) segments(5,3,5,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(31,39))) point(6,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(32,39))) segments(6,4,7,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(33,39))) point(7,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(38,39))) segments(6,3,6,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(32,40))) point(7,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(33,40))) segments(7,4,8,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(39,40))) segments(7,3,7,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(33,41))) point(8,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(34,41))) point(9,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(40,41))) segments(8,3,8,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(34,42))) segments(9,4,10,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(41,42))) segments(9,3,9,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(35,43))) point(3,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(36,43))) segments(3,3,4,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(37,43))) point(4,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(36,44))) point(4,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(37,44))) segments(4,3,5,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(38,44))) point(5,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(43,44))) segments(4,2,4,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(37,45))) point(5,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(38,45))) segments(5,3,6,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(39,45))) point(6,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(44,45))) segments(5,2,5,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(38,46))) point(6,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(39,46))) segments(6,3,7,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(40,46))) point(7,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(45,46))) segments(6,2,6,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(39,47))) point(7,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(40,47))) segments(7,3,8,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(41,47))) point(8,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(46,47))) segments(7,2,7,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(40,48))) point(8,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(41,48))) segments(8,3,9,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(42,48))) point(9,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(47,48))) segments(8,2,8,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(43,49))) point(4,2,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(44,49))) segments(4,2,5,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(45,49))) point(5,2,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(44,50))) point(5,2,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(45,50))) segments(5,2,6,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(46,50))) point(6,2,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(49,50))) segments(5,1,5,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(45,51))) point(6,2,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(46,51))) segments(6,2,7,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(47,51))) point(7,2,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(50,51))) segments(6,1,6,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(46,52))) point(7,2,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(47,52))) segments(7,2,8,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(48,52))) point(8,2,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(51,52))) segments(7,1,7,2, lwd=lwd,lend=lend,col=colorAdj)

    }

    ###Plot corners
    for (i in WhichCorners) {
      Indeces<-as.numeric(Boundary[i,1:2])
      colorAdj<-color.adj[i]

      if (identical(Indeces,c(1,2))) segments(5,8,5,9,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(2,3))) segments(6,8,6,9,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(3,4))) segments(7,8,7,9,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(1,5))) point(4,8,ulbr=FALSE, colorAdj)
      if (identical(Indeces,c(1,6))) segments(4,8,5,8,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(2,6))) point(5,8,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(5,6))) segments(4,7,4,8,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(1,7))) point(5,8,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(2,7))) segments(5,8,6,8,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(3,7))) point(6,8,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(6,7))) segments(5,7,5,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(2,8))) point(6,8,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(3,8))) segments(6,8,7,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(4,8))) point(7,8,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(7,8))) segments(6,7,6,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(3,9))) point(7,8,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(4,9))) segments(7,8,8,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(8,9))) segments(7,7,7,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(4,10))) point(8,8,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(9,10))) segments(8,7,8,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(5,11))) point(3,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(5,12))) segments(3,7,4,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(6,12))) point(4,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(11,12))) segments(3,6,3,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(5,13))) point(4,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(6,13))) segments(4,7,5,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(7,13))) point(5,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(12,13))) segments(4,6,4,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(6,14))) point(5,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(7,14))) segments(5,7,6,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(8,14))) point(6,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(13,14))) segments(5,6,5,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(7,15))) point(6,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(8,15))) segments(6,7,7,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(9,15))) point(7,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(14,15))) segments(6,6,6,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(8,16))) point(7,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(9,16))) segments(7,7,8,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(10,16))) point(8,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(15,16))) segments(7,6,7,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(9,17))) point(8,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(10,17))) segments(8,7,9,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(16,17))) segments(8,6,8,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(10,18))) point(9,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(17,18))) segments(9,6,9,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(11,19))) point(2,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(11,20))) segments(2,6,3,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(12,20))) point(3,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(19,20))) segments(2,5,2,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(11,21))) point(3,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(12,21))) segments(3,6,4,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(13,21))) point(4,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(20,21))) segments(3,5,3,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(12,22))) point(4,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(13,22))) segments(4,6,5,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(14,22))) point(5,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(21,22))) segments(4,5,4,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(13,23))) point(5,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(14,23))) segments(5,6,6,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(15,23))) point(6,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(22,23))) segments(5,5,5,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(14,24))) point(6,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(15,24))) segments(6,6,7,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(16,24))) point(7,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(23,24))) segments(6,5,6,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(15,25))) point(7,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(16,25))) segments(7,6,8,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(17,25))) point(8,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(24,25))) segments(7,5,7,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(17,26))) point(9,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(18,26))) segments(9,6,10,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(19,27))) segments(1,5,2,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(20,27))) point(2,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(19,28))) point(2,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(20,28))) segments(2,5,3,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(21,28))) point(3,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(27,28))) segments(2,4,2,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(20,29))) point(3,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(21,29))) segments(3,5,4,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(22,29))) point(4,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(28,29))) segments(3,4,3,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(21,30))) point(4,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(22,30))) segments(4,5,5,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(23,30))) point(5,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(29,30))) segments(4,4,4,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(22,31))) point(5,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(23,31))) segments(5,5,6,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(24,31))) point(6,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(30,31))) segments(5,4,5,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(23,32))) point(6,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(24,32))) segments(6,5,7,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(25,32))) point(7,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(31,32))) segments(6,4,6,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(24,33))) point(7,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(25,33))) segments(7,5,8,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(32,33))) segments(7,4,7,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(26,34))) segments(9,5,10,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(27,35))) point(2,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(28,35))) segments(2,4,3,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(29,35))) point(3,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(28,36))) point(3,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(29,36))) segments(3,4,4,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(30,36))) point(4,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(35,36))) segments(3,3,3,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(29,37))) point(4,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(30,37))) segments(4,4,5,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(31,37))) point(5,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(36,37))) segments(4,3,4,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(30,38))) point(5,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(31,38))) segments(5,4,6,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(32,38))) point(6,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(37,38))) segments(5,3,5,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(31,39))) point(6,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(32,39))) segments(6,4,7,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(33,39))) point(7,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(38,39))) segments(6,3,6,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(32,40))) point(7,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(33,40))) segments(7,4,8,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(39,40))) segments(7,3,7,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(33,41))) point(8,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(34,41))) point(9,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(40,41))) segments(8,3,8,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(34,42))) segments(9,4,10,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(41,42))) segments(9,3,9,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(35,43))) point(3,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(36,43))) segments(3,3,4,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(37,43))) point(4,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(36,44))) point(4,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(37,44))) segments(4,3,5,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(38,44))) point(5,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(43,44))) segments(4,2,4,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(37,45))) point(5,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(38,45))) segments(5,3,6,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(39,45))) point(6,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(44,45))) segments(5,2,5,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(38,46))) point(6,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(39,46))) segments(6,3,7,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(40,46))) point(7,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(45,46))) segments(6,2,6,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(39,47))) point(7,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(40,47))) segments(7,3,8,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(41,47))) point(8,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(46,47))) segments(7,2,7,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(40,48))) point(8,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(41,48))) segments(8,3,9,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(42,48))) point(9,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(47,48))) segments(8,2,8,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(43,49))) point(4,2,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(44,49))) segments(4,2,5,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(45,49))) point(5,2,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(44,50))) point(5,2,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(45,50))) segments(5,2,6,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(46,50))) point(6,2,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(49,50))) segments(5,1,5,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(45,51))) point(6,2,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(46,51))) segments(6,2,7,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(47,51))) point(7,2,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(50,51))) segments(6,1,6,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(46,52))) point(7,2,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(47,52))) segments(7,2,8,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(48,52))) point(8,2,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(51,52))) segments(7,1,7,2, lwd=lwd,lend=lend,col=colorAdj)

    }

    ###Add border
    hloop<-list(4:7,c(3,8),c(2,9),1,NULL,1,c(2,9),c(3,8),4:7)
    vloop<-list(4:5,c(3,6),c(2,7),c(1,8),NULL,NULL,NULL,c(1,8),c(2,7),3:6)
    for (j in 1:9) {
      for (i in hloop[[j]]) {
        segments(i,j,i+1,j,lwd=lwd.border)
      }
    }
    for (i in 1:10) {
      for (j in vloop[[i]]) {
        segments(i,j,i,j+1,lwd=lwd.border)
      }
    }

    ###Add blind spot
    rect(8,4,9,6,col=color.bs,border=color.bs)

    ###Add legend
    if (legend) {

      NColors<-length(col.pal)
      Vertical<-seq(3,7,length.out=NColors)
      if (stat == "mean") for (i in 1:NColors) segments(11,Vertical[i],11.75,Vertical[i],col=rev(col.pal)[i],lwd=1.5)
      if (stat == "sd") for (i in 1:NColors) segments(11,Vertical[i],11.75,Vertical[i],col=(col.pal)[i],lwd=1.5)
      minx<-zlim[1]
      maxx<-zlim[2]
      LegendPV<-seq(minx,maxx,length.out=5)
      segments(11.75,3,11.75,7,lwd=2)
      segments(11,3,11,7,lwd=2)
      segments(11,7,11.75,7,lwd=2)
      segments(11,3,11.75,3,lwd=2)
      for (i in 1:length(LegendPV)) {
        if ((stat == "mean") & (!is.na(Visit))) {
          if (is.na(Visit)) text(12.75,(3:7)[i],format0(rev(LegendPV)[i]))
          if (!is.na(Visit)) text(12.75,(3:7)[i],format2(rev(LegendPV)[i]))
        }
        if ((stat == "sd")& ((!is.na(Visit)))) {
          if (is.na(Visit)) text(12.75,(3:7)[i],format0((LegendPV)[i]))
          if (!is.na(Visit)) text(12.75,(3:7)[i],format2((LegendPV)[i]))
        }
        if (is.na(Visit)) {
          if (is.na(Visit)) text(12.75,(3:7)[i],format0((LegendPV)[i]))
          if (!is.na(Visit)) text(12.75,(3:7)[i],format2((LegendPV)[i]))
        }
        segments(11.75,(3:7)[i],12,(3:7)[i],lwd=2)
      }
      if (is.na(Visit)) text(11.5, 7.5, expression(paste("Degree (",degree,")")))
      if (!is.na(Visit) & stat == "mean") text(11.4, 7.5, expression(paste("E[",w[ij],"(",alpha[t],")]")))
      if (!is.na(Visit) & stat == "sd") text(11.4, 7.5, expression(paste("sd[",w[ij],"(",alpha[t],")]")))
    }

    if (!is.null(Degree)) {

      ###Add garway angles

      text(4.5,8.5,Degree[1])
      text(5.5,8.5,Degree[2])
      text(6.5,8.5,Degree[3])
      text(7.5,8.5,Degree[4])

      text(3.5,7.5,Degree[5])
      text(4.5,7.5,Degree[6])
      text(5.5,7.5,Degree[7])
      text(6.5,7.5,Degree[8])
      text(7.5,7.5,Degree[9])
      text(8.5,7.5,Degree[10])

      text(2.5,6.5,Degree[11])
      text(3.5,6.5,Degree[12])
      text(4.5,6.5,Degree[13])
      text(5.5,6.5,Degree[14])
      text(6.5,6.5,Degree[15])
      text(7.5,6.5,Degree[16])
      text(8.5,6.5,Degree[17])
      text(9.5,6.5,Degree[18])

      text(1.5,5.5,Degree[19])
      text(2.5,5.5,Degree[20])
      text(3.5,5.5,Degree[21])
      text(4.5,5.5,Degree[22])
      text(5.5,5.5,Degree[23])
      text(6.5,5.5,Degree[24])
      text(7.5,5.5,Degree[25])
      text(9.5,5.5,Degree[26])

      text(1.5,4.5,Degree[27])
      text(2.5,4.5,Degree[28])
      text(3.5,4.5,Degree[29])
      text(4.5,4.5,Degree[30])
      text(5.5,4.5,Degree[31])
      text(6.5,4.5,Degree[32])
      text(7.5,4.5,Degree[33])
      text(9.5,4.5,Degree[34])

      text(2.5,3.5,Degree[35])
      text(3.5,3.5,Degree[36])
      text(4.5,3.5,Degree[37])
      text(5.5,3.5,Degree[38])
      text(6.5,3.5,Degree[39])
      text(7.5,3.5,Degree[40])
      text(8.5,3.5,Degree[41])
      text(9.5,3.5,Degree[42])

      text(3.5,2.5,Degree[43])
      text(4.5,2.5,Degree[44])
      text(5.5,2.5,Degree[45])
      text(6.5,2.5,Degree[46])
      text(7.5,2.5,Degree[47])
      text(8.5,2.5,Degree[48])

      text(4.5,1.5,Degree[49])
      text(5.5,1.5,Degree[50])
      text(6.5,1.5,Degree[51])
      text(7.5,1.5,Degree[52])

    }

    ###Return par to default
    par(pardefault)

  }

  ###Check inputs
  if (!missing(Wij)) {

    if (!is.PosteriorAdj(Wij)) stop('"Wij" is not a PosteriorAdj object')
    Nu <- as.numeric(unlist(strsplit(colnames(Wij)[dim(Wij)[2]], "sd"))[2])
    if (!(Visit %in% c(NA,1:Nu))) stop('"Visit" must be either NA or an integer between 1 and Nu')
    if (!(stat %in% c("mean", "sd"))) stop('"stat" must be one of "mean" or "sd"')
    if (!any(areColors(color.scheme))) stop('"color.scheme" can only include colors')

    ###Border information
    if (is.na(Visit)) {
      Boundary <- Wij[ , 1 : 3]
      Boundary[ ,3] <- 180 - Boundary[ ,3]
    }
    if (!is.na(Visit) & stat == "mean") Boundary <- Wij[ , c(1, 2, (2 * Visit) + 2)]
    if (!is.na(Visit) & stat == "sd") Boundary <- Wij[ , c(1, 2, (2 * Visit) + 3)]
    NAdjacency <- dim(Boundary)[1]
    Degree <- DM

    ###Set color pallete
    col.breaks <- seq(zlim[1],zlim[2],length.out=(length(unique(Boundary[,3]))+1))
    col.br <- colorRampPalette(color.scheme)
    col.pal <- col.br(NAdjacency)
    # suppressWarnings(fixed_obs<-classIntervals(Boundary[,3],n=length(unique(Boundary[,3]))))
    suppressWarnings(fixed_obs<-classInt::classIntervals(Boundary[,3],style="fixed",fixedBreaks=col.breaks))
    color.adj<-classInt::findColours(fixed_obs,col.pal)

    ###Plotting functions and Parameters
    lwd<-edgewidth
    lend<-2
    format2<-function(x) format(round(x,2),nsmall=2)
    format0<-function(x) format(round(x,0),nsmall=0)
    tri<-cornerwidth
    point<-function(x,y,ulbr=TRUE, col=colorAdj) {
      if (!ulbr) {
        polygon(c(x,x+tri,x),c(y,y,y+tri),col=col,border=col)
        polygon(c(x,x-tri,x),c(y,y,y-tri),col=col,border=col)
      }
      if (ulbr) {
        polygon(c(x,x-tri,x),c(y,y,y+tri),col=col,border=col)
        polygon(c(x,x+tri,x),c(y,y,y-tri),col=col,border=col)
      }
    }

    ###Create edge indicator
    EdgesInd<-numeric(length=NAdjacency)
    WhichEdges<-c(1,2,3,5,7,9,11,13,15,17,18,20,22,24,26,28,30,32,34,36,38,40,42,43,45,47,49,51,53,55,57,59,61,63,65,67, 69,71,72,75,77,79,81,83,85,87,89,91,93,95,96,97,99,102,104,106,108,110,112,114,116,118,119,122,123,124,126,129,131,133,135,137,139,141,143,145,147,149,152,154,156,158,160,162)
    EdgesInd[WhichEdges]<-rep(1,length(WhichEdges))
    WhichCorners<-which(EdgesInd==0)

    ###Create plot
    pardefault <- suppressWarnings(par(no.readonly = T))
    par(mfcol = c(1, 1), pty = "m", mai = c(0, 0, 0.75, 0))
    plot(1,1,xlim=c(1,13),ylim=c(1.25,8.75),type="n",xaxt="n",yaxt="n",bty="n",ylab="",xlab="",main=main, asp = 1, cex.main = 1.5)

    ###Plot edges
    for (i in WhichEdges) {
      Indeces<-as.numeric(Boundary[i,1:2])
      colorAdj<-color.adj[i]

      if (identical(Indeces,c(1,2))) segments(5,8,5,9,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(2,3))) segments(6,8,6,9,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(3,4))) segments(7,8,7,9,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(1,5))) point(4,8,ulbr=FALSE, colorAdj)
      if (identical(Indeces,c(1,6))) segments(4,8,5,8,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(2,6))) point(5,8,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(5,6))) segments(4,7,4,8,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(1,7))) point(5,8,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(2,7))) segments(5,8,6,8,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(3,7))) point(6,8,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(6,7))) segments(5,7,5,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(2,8))) point(6,8,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(3,8))) segments(6,8,7,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(4,8))) point(7,8,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(7,8))) segments(6,7,6,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(3,9))) point(7,8,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(4,9))) segments(7,8,8,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(8,9))) segments(7,7,7,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(4,10))) point(8,8,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(9,10))) segments(8,7,8,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(5,11))) point(3,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(5,12))) segments(3,7,4,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(6,12))) point(4,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(11,12))) segments(3,6,3,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(5,13))) point(4,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(6,13))) segments(4,7,5,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(7,13))) point(5,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(12,13))) segments(4,6,4,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(6,14))) point(5,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(7,14))) segments(5,7,6,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(8,14))) point(6,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(13,14))) segments(5,6,5,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(7,15))) point(6,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(8,15))) segments(6,7,7,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(9,15))) point(7,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(14,15))) segments(6,6,6,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(8,16))) point(7,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(9,16))) segments(7,7,8,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(10,16))) point(8,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(15,16))) segments(7,6,7,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(9,17))) point(8,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(10,17))) segments(8,7,9,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(16,17))) segments(8,6,8,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(10,18))) point(9,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(17,18))) segments(9,6,9,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(11,19))) point(2,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(11,20))) segments(2,6,3,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(12,20))) point(3,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(19,20))) segments(2,5,2,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(11,21))) point(3,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(12,21))) segments(3,6,4,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(13,21))) point(4,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(20,21))) segments(3,5,3,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(12,22))) point(4,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(13,22))) segments(4,6,5,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(14,22))) point(5,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(21,22))) segments(4,5,4,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(13,23))) point(5,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(14,23))) segments(5,6,6,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(15,23))) point(6,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(22,23))) segments(5,5,5,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(14,24))) point(6,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(15,24))) segments(6,6,7,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(16,24))) point(7,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(23,24))) segments(6,5,6,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(15,25))) point(7,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(16,25))) segments(7,6,8,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(17,25))) point(8,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(24,25))) segments(7,5,7,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(17,26))) point(9,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(18,26))) segments(9,6,10,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(19,27))) segments(1,5,2,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(20,27))) point(2,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(19,28))) point(2,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(20,28))) segments(2,5,3,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(21,28))) point(3,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(27,28))) segments(2,4,2,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(20,29))) point(3,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(21,29))) segments(3,5,4,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(22,29))) point(4,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(28,29))) segments(3,4,3,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(21,30))) point(4,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(22,30))) segments(4,5,5,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(23,30))) point(5,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(29,30))) segments(4,4,4,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(22,31))) point(5,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(23,31))) segments(5,5,6,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(24,31))) point(6,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(30,31))) segments(5,4,5,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(23,32))) point(6,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(24,32))) segments(6,5,7,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(25,32))) point(7,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(31,32))) segments(6,4,6,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(24,33))) point(7,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(25,33))) segments(7,5,8,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(32,33))) segments(7,4,7,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(26,34))) segments(9,5,10,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(27,35))) point(2,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(28,35))) segments(2,4,3,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(29,35))) point(3,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(28,36))) point(3,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(29,36))) segments(3,4,4,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(30,36))) point(4,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(35,36))) segments(3,3,3,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(29,37))) point(4,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(30,37))) segments(4,4,5,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(31,37))) point(5,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(36,37))) segments(4,3,4,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(30,38))) point(5,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(31,38))) segments(5,4,6,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(32,38))) point(6,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(37,38))) segments(5,3,5,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(31,39))) point(6,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(32,39))) segments(6,4,7,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(33,39))) point(7,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(38,39))) segments(6,3,6,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(32,40))) point(7,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(33,40))) segments(7,4,8,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(39,40))) segments(7,3,7,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(33,41))) point(8,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(34,41))) point(9,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(40,41))) segments(8,3,8,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(34,42))) segments(9,4,10,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(41,42))) segments(9,3,9,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(35,43))) point(3,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(36,43))) segments(3,3,4,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(37,43))) point(4,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(36,44))) point(4,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(37,44))) segments(4,3,5,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(38,44))) point(5,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(43,44))) segments(4,2,4,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(37,45))) point(5,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(38,45))) segments(5,3,6,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(39,45))) point(6,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(44,45))) segments(5,2,5,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(38,46))) point(6,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(39,46))) segments(6,3,7,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(40,46))) point(7,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(45,46))) segments(6,2,6,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(39,47))) point(7,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(40,47))) segments(7,3,8,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(41,47))) point(8,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(46,47))) segments(7,2,7,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(40,48))) point(8,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(41,48))) segments(8,3,9,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(42,48))) point(9,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(47,48))) segments(8,2,8,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(43,49))) point(4,2,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(44,49))) segments(4,2,5,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(45,49))) point(5,2,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(44,50))) point(5,2,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(45,50))) segments(5,2,6,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(46,50))) point(6,2,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(49,50))) segments(5,1,5,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(45,51))) point(6,2,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(46,51))) segments(6,2,7,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(47,51))) point(7,2,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(50,51))) segments(6,1,6,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(46,52))) point(7,2,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(47,52))) segments(7,2,8,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(48,52))) point(8,2,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(51,52))) segments(7,1,7,2, lwd=lwd,lend=lend,col=colorAdj)

    }

    ###Plot corners
    for (i in WhichCorners) {
      Indeces<-as.numeric(Boundary[i,1:2])
      colorAdj<-color.adj[i]

      if (identical(Indeces,c(1,2))) segments(5,8,5,9,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(2,3))) segments(6,8,6,9,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(3,4))) segments(7,8,7,9,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(1,5))) point(4,8,ulbr=FALSE, colorAdj)
      if (identical(Indeces,c(1,6))) segments(4,8,5,8,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(2,6))) point(5,8,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(5,6))) segments(4,7,4,8,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(1,7))) point(5,8,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(2,7))) segments(5,8,6,8,lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(3,7))) point(6,8,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(6,7))) segments(5,7,5,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(2,8))) point(6,8,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(3,8))) segments(6,8,7,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(4,8))) point(7,8,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(7,8))) segments(6,7,6,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(3,9))) point(7,8,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(4,9))) segments(7,8,8,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(8,9))) segments(7,7,7,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(4,10))) point(8,8,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(9,10))) segments(8,7,8,8, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(5,11))) point(3,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(5,12))) segments(3,7,4,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(6,12))) point(4,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(11,12))) segments(3,6,3,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(5,13))) point(4,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(6,13))) segments(4,7,5,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(7,13))) point(5,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(12,13))) segments(4,6,4,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(6,14))) point(5,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(7,14))) segments(5,7,6,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(8,14))) point(6,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(13,14))) segments(5,6,5,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(7,15))) point(6,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(8,15))) segments(6,7,7,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(9,15))) point(7,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(14,15))) segments(6,6,6,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(8,16))) point(7,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(9,16))) segments(7,7,8,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(10,16))) point(8,7,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(15,16))) segments(7,6,7,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(9,17))) point(8,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(10,17))) segments(8,7,9,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(16,17))) segments(8,6,8,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(10,18))) point(9,7,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(17,18))) segments(9,6,9,7, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(11,19))) point(2,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(11,20))) segments(2,6,3,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(12,20))) point(3,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(19,20))) segments(2,5,2,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(11,21))) point(3,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(12,21))) segments(3,6,4,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(13,21))) point(4,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(20,21))) segments(3,5,3,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(12,22))) point(4,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(13,22))) segments(4,6,5,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(14,22))) point(5,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(21,22))) segments(4,5,4,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(13,23))) point(5,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(14,23))) segments(5,6,6,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(15,23))) point(6,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(22,23))) segments(5,5,5,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(14,24))) point(6,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(15,24))) segments(6,6,7,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(16,24))) point(7,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(23,24))) segments(6,5,6,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(15,25))) point(7,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(16,25))) segments(7,6,8,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(17,25))) point(8,6,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(24,25))) segments(7,5,7,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(17,26))) point(9,6,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(18,26))) segments(9,6,10,6, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(19,27))) segments(1,5,2,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(20,27))) point(2,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(19,28))) point(2,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(20,28))) segments(2,5,3,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(21,28))) point(3,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(27,28))) segments(2,4,2,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(20,29))) point(3,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(21,29))) segments(3,5,4,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(22,29))) point(4,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(28,29))) segments(3,4,3,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(21,30))) point(4,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(22,30))) segments(4,5,5,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(23,30))) point(5,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(29,30))) segments(4,4,4,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(22,31))) point(5,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(23,31))) segments(5,5,6,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(24,31))) point(6,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(30,31))) segments(5,4,5,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(23,32))) point(6,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(24,32))) segments(6,5,7,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(25,32))) point(7,5,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(31,32))) segments(6,4,6,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(24,33))) point(7,5,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(25,33))) segments(7,5,8,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(32,33))) segments(7,4,7,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(26,34))) segments(9,5,10,5, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(27,35))) point(2,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(28,35))) segments(2,4,3,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(29,35))) point(3,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(28,36))) point(3,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(29,36))) segments(3,4,4,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(30,36))) point(4,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(35,36))) segments(3,3,3,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(29,37))) point(4,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(30,37))) segments(4,4,5,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(31,37))) point(5,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(36,37))) segments(4,3,4,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(30,38))) point(5,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(31,38))) segments(5,4,6,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(32,38))) point(6,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(37,38))) segments(5,3,5,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(31,39))) point(6,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(32,39))) segments(6,4,7,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(33,39))) point(7,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(38,39))) segments(6,3,6,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(32,40))) point(7,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(33,40))) segments(7,4,8,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(39,40))) segments(7,3,7,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(33,41))) point(8,4,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(34,41))) point(9,4,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(40,41))) segments(8,3,8,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(34,42))) segments(9,4,10,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(41,42))) segments(9,3,9,4, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(35,43))) point(3,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(36,43))) segments(3,3,4,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(37,43))) point(4,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(36,44))) point(4,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(37,44))) segments(4,3,5,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(38,44))) point(5,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(43,44))) segments(4,2,4,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(37,45))) point(5,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(38,45))) segments(5,3,6,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(39,45))) point(6,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(44,45))) segments(5,2,5,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(38,46))) point(6,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(39,46))) segments(6,3,7,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(40,46))) point(7,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(45,46))) segments(6,2,6,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(39,47))) point(7,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(40,47))) segments(7,3,8,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(41,47))) point(8,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(46,47))) segments(7,2,7,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(40,48))) point(8,3,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(41,48))) segments(8,3,9,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(42,48))) point(9,3,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(47,48))) segments(8,2,8,3, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(43,49))) point(4,2,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(44,49))) segments(4,2,5,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(45,49))) point(5,2,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(44,50))) point(5,2,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(45,50))) segments(5,2,6,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(46,50))) point(6,2,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(49,50))) segments(5,1,5,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(45,51))) point(6,2,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(46,51))) segments(6,2,7,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(47,51))) point(7,2,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(50,51))) segments(6,1,6,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(46,52))) point(7,2,ulbr=TRUE,  colorAdj)
      if (identical(Indeces,c(47,52))) segments(7,2,8,2, lwd=lwd,lend=lend,col=colorAdj)
      if (identical(Indeces,c(48,52))) point(8,2,ulbr=FALSE,  colorAdj)
      if (identical(Indeces,c(51,52))) segments(7,1,7,2, lwd=lwd,lend=lend,col=colorAdj)

    }

    ###Add border
    hloop<-list(4:7,c(3,8),c(2,9),1,NULL,1,c(2,9),c(3,8),4:7)
    vloop<-list(4:5,c(3,6),c(2,7),c(1,8),NULL,NULL,NULL,c(1,8),c(2,7),3:6)
    for (j in 1:9) {
      for (i in hloop[[j]]) {
        segments(i,j,i+1,j,lwd=lwd.border)
      }
    }
    for (i in 1:10) {
      for (j in vloop[[i]]) {
        segments(i,j,i,j+1,lwd=lwd.border)
      }
    }

    ###Add blind spot
    rect(8,4,9,6,col=color.bs,border=color.bs)

    ###Add legend
    if (legend) {

      NColors<-length(col.pal)
      Vertical<-seq(3,7,length.out=NColors)
      if (stat == "mean") for (i in 1:NColors) segments(11,Vertical[i],11.75,Vertical[i],col=rev(col.pal)[i],lwd=1.5)
      if (stat == "sd") for (i in 1:NColors) segments(11,Vertical[i],11.75,Vertical[i],col=(col.pal)[i],lwd=1.5)
      minx<-zlim[1]
      maxx<-zlim[2]
      LegendPV<-seq(minx,maxx,length.out=5)
      segments(11.75,3,11.75,7,lwd=2)
      segments(11,3,11,7,lwd=2)
      segments(11,7,11.75,7,lwd=2)
      segments(11,3,11.75,3,lwd=2)
      for (i in 1:length(LegendPV)) {
        if ((stat == "mean") & (!is.na(Visit))) {
          if (is.na(Visit)) text(12.75,(3:7)[i],format0(rev(LegendPV)[i]))
          if (!is.na(Visit)) text(12.75,(3:7)[i],format2(rev(LegendPV)[i]))
        }
        if ((stat == "sd")& ((!is.na(Visit)))) {
          if (is.na(Visit)) text(12.75,(3:7)[i],format0((LegendPV)[i]))
          if (!is.na(Visit)) text(12.75,(3:7)[i],format2((LegendPV)[i]))
        }
        if (is.na(Visit)) {
          if (is.na(Visit)) text(12.75,(3:7)[i],format0((LegendPV)[i]))
          if (!is.na(Visit)) text(12.75,(3:7)[i],format2((LegendPV)[i]))
        }
        segments(11.75,(3:7)[i],12,(3:7)[i],lwd=2)
      }
      if (is.na(Visit)) text(11.5, 7.5, expression(paste("Degree (",degree,")")))
      if (!is.na(Visit) & stat == "mean") text(11.4, 7.5, expression(paste("E[",w[ij],"(",alpha[t],")]")))
      if (!is.na(Visit) & stat == "sd") text(11.4, 7.5, expression(paste("sd[",w[ij],"(",alpha[t],")]")))
    }

    if (!is.null(Degree)) {

      ###Add garway angles

      text(4.5,8.5,Degree[1])
      text(5.5,8.5,Degree[2])
      text(6.5,8.5,Degree[3])
      text(7.5,8.5,Degree[4])

      text(3.5,7.5,Degree[5])
      text(4.5,7.5,Degree[6])
      text(5.5,7.5,Degree[7])
      text(6.5,7.5,Degree[8])
      text(7.5,7.5,Degree[9])
      text(8.5,7.5,Degree[10])

      text(2.5,6.5,Degree[11])
      text(3.5,6.5,Degree[12])
      text(4.5,6.5,Degree[13])
      text(5.5,6.5,Degree[14])
      text(6.5,6.5,Degree[15])
      text(7.5,6.5,Degree[16])
      text(8.5,6.5,Degree[17])
      text(9.5,6.5,Degree[18])

      text(1.5,5.5,Degree[19])
      text(2.5,5.5,Degree[20])
      text(3.5,5.5,Degree[21])
      text(4.5,5.5,Degree[22])
      text(5.5,5.5,Degree[23])
      text(6.5,5.5,Degree[24])
      text(7.5,5.5,Degree[25])
      text(9.5,5.5,Degree[26])

      text(1.5,4.5,Degree[27])
      text(2.5,4.5,Degree[28])
      text(3.5,4.5,Degree[29])
      text(4.5,4.5,Degree[30])
      text(5.5,4.5,Degree[31])
      text(6.5,4.5,Degree[32])
      text(7.5,4.5,Degree[33])
      text(9.5,4.5,Degree[34])

      text(2.5,3.5,Degree[35])
      text(3.5,3.5,Degree[36])
      text(4.5,3.5,Degree[37])
      text(5.5,3.5,Degree[38])
      text(6.5,3.5,Degree[39])
      text(7.5,3.5,Degree[40])
      text(8.5,3.5,Degree[41])
      text(9.5,3.5,Degree[42])

      text(3.5,2.5,Degree[43])
      text(4.5,2.5,Degree[44])
      text(5.5,2.5,Degree[45])
      text(6.5,2.5,Degree[46])
      text(7.5,2.5,Degree[47])
      text(8.5,2.5,Degree[48])

      text(4.5,1.5,Degree[49])
      text(5.5,1.5,Degree[50])
      text(6.5,1.5,Degree[51])
      text(7.5,1.5,Degree[52])

    }

    ###Return par to default
    par(pardefault)

  }

###End Function
}
