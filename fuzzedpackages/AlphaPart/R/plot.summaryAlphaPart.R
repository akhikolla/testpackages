#' plot.summaryAlphaPart.R
#'
#' A function to plot summary of partitioned breeding values.
#'
#' @details
#' Information in summaries of partitions of breeding values can be
#' overhelming due to a large volume of numbers. Plot method can be used to
#' visualise this data in eye pleasing way using ggplot2
#' graphics.
#'
#' @param x summaryAlphaPart, object from the \code{AlphaPart(...)} or \code{summary(AlphaPart(...), ...)} call.
#' @param sortValue Logical, affect legend attributes via sort of paths according to
#' \code{sortValueFUN} function; if not logical, then ordered paths are given as a character vector.
#' @param sortValueFUN Function, that produces single value for one vector, say \code{mean} or \code{sum}.
#' @param sortValueDec Logical, sort decreasing.
#' @param addSum Logical, plot the overall trend.
#' @param paths Character or list or characters, name of paths to plot; if \code{NULL} plot all paths; see examples.
#' @param xlab Character, x-axis label.
#' @param ylab Character, y-axis label; can be a vector of several labels if there are more traits in \code{x} (recycled!).
#' @param xlim Numeric, a vector of two values with x-axis limits; use a list of vectors for more traits.
#' @param ylim Numeric, a vector of two values with y-axis limits; use a list of vectors for more traits.
#' @param color Character, color names; by default a set of 54 colors is predefined from the \pkg{RColorBrewer} package;
#' in addition a black colour is attached at the begining for the overall trend; if there are more paths than
#' colors then recycling occours.
#' @param lineSize Numeric, line width.
#' @param lineType Numeric, line type (recycled); can be used only if lineTypeList=NULL.
#' @param lineTypeList List, named list of numeric values that help to point out a set of paths
#' (distinguished with line type) within upper level of paths (distinguished by,
#' color), e.g., lineTypeList=list("-1"=1, "-2"=2, def=1) will lead to use
#' of line type 1 for paths having "-1" at the end of path name and line type 2,
#' for paths having "-2" at the end of path name, while line type 1 (default) will,
#' be used for other paths; specification of this argument also causes recycling
#' of colors for the upper level of paths; if NULL all lines have a standard line type,
#' otherwise \code{lineType} does not have any effect.
#' @param useDirectLabels Logical, use directlabels package for legend.
#' @param method List, method for direct.label.
#' @param labelPath  Character, legend title; used only if \code{useDirectLabels=FALSE}.
#' @param ...  Arguments passed to other functions (not used at the moment).
#'
#' @example inst/examples/examples_plotSummaryAlphaPart.R
#'
#' @return A list of ggplot objects that can be further modified or displayed.
#' For each trait in \code{x} there is one plot visualising summarized values.
#'
#' @useDynLib AlphaPart, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @export
#' @importFrom directlabels last.qp
#' @importFrom directlabels direct.label
#' @importFrom reshape melt
#' @import ggplot2


plot.summaryAlphaPart <- function (
  x,
  sortValue=TRUE,
  sortValueFUN=sum,
  sortValueDec=TRUE,
  addSum=TRUE,
  paths=NULL,
  xlab=NULL,
  ylab=NULL,
  xlim=NULL,
  ylim=NULL,
  color,
  lineSize=1,
  lineType=1,
  lineTypeList=NULL,
  useDirectLabels=TRUE,
  method=list(last.qp, hjust=0),
  labelPath=NULL,
  ...
) {

  ## --- Setup ---

  if (!("summaryAlphaPart" %in% class(x))) stop("'x' must be of a summaryAlphaPart class")

  by    <- x$info$by
  path  <- x$info$path
  lT    <- x$info$lT
  nT    <- x$info$nT
  nP    <- x$info$nP
  ret   <- vector(mode="list", length=nT)
  names(ret) <- x$info$lT

  ## Axis labels
  if (!is.null(xlab) && length(xlab) > 1) stop("you can provide only one value for 'xlab'")
  if (!is.null(ylab) && length(ylab) < nT) ylab <- rep(ylab, length=nT)

  ## Colors
  if (!missing(color)) {
    if (length(color) < nP) color <- rep(color, length=nP)
      color <- c("black", color)
  } else {
    if (FALSE) { ## Code to generate a bunch of qualitative colors
      requireNamespace("RColorBrewer")
      #library(package="RColorBrewer")
      pals <- c("Set1", "Dark2", "Accent", "Paired", "Set2", "Set3")
      palsN <- brewer.pal.info[pals, "maxcolors"]
      color <- vector(length=sum(palsN))
      j <- 1
      for (i in seq(along=pals)) {
        color[j:(j - 1 + palsN[i])] <- do.call("brewer.pal", args=list(n=palsN[i], name=pals[i]))
        j <- j + palsN[i]
      }
      color <- unique(color)
    }
    color <- c("black",
               "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
               "#A65628", "#F781BF", "#999999", "#1B9E77", "#D95F02", "#7570B3",
               "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#7FC97F",
               "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
               "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
               "#FDBF6F", "#CAB2D6", "#6A3D9A", "#B15928", "#66C2A5", "#FC8D62",
               "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3",
               "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
               "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")
    color <- color[color != "#FFFF33"] ## remove yellow color(s)
  }

  ## Line type
  if (is.null(lineTypeList)) {
    if (length(lineType) < nP) {
      lineType <- c(1, rep(x=lineType, times=nP))
    } else {
      lineType <- c(1, lineType)
    }
  }

  ## --- Create plots ---


  ## Ylim
  ylimT <- ylim


  for (i in 1:nT) { ## loop over traits
    colorI <- color
    lineTypeI <- lineType
    ## Prepare data
    tmp0 <- x[[i]]
    ## Make sure that path has a name "N"
    tmpCol <- colnames(tmp0)
    test <- tmpCol %in% "N"
    if (sum(test) > 1) {
      tmpCol[test & cumsum(test) == 2] <- "N."
      colnames(tmp0) <- tmpCol
      warning("changing path name from 'N' to 'N.'")
    }
    tmp <- melt(tmp0[, !(colnames(tmp0) %in% "N")], id=by)
    colnames(tmp) <- c("by", "path", "trait")
    if (is.logical(sortValue)) {
      if (sortValue) {
        nC <- ncol(tmp0)
        pathStat <- sapply(X=tmp0[, (nC - nP + 1):nC], FUN=sortValueFUN, na.rm=TRUE)
        levs <- names(sort(pathStat, decreasing=sortValueDec))
        tmp$path <- factor(tmp$path, levels=c(x$info$labelSum, levs))
        if (!is.null(lineTypeList)) { ## fiddle with upper (color) and lower (line type) level of paths
          levs2X <- names(lineTypeList); levs2X <- levs2X[levs2X != "def"]
          levs1 <- levs2 <- levs
            for (k in levs2X) {
            j <- paste(k, "$", sep="") ## lower label mark needs to be at the end of path name!!!
            levs2[grep(pattern=j, x=levs1)] <- k
            levs1 <- sub(pattern=j, replacement="", x=levs1)
          }
          levs2[!levs2 %in% levs2X] <- "def"
          levs1X <- unique(levs1)
          colorList <- as.list(color[-1][1:length(levs1X)]); names(colorList) <- levs1X
          colorI <- c("black", unlist(colorList[levs1])); names(colorI) <- NULL
          lineTypeI <- c(lineTypeList$def, unlist(lineTypeList[levs2])); names(lineTypeI) <- NULL
        }
      }
    } else {
      tmp$path <- factor(tmp$path, levels=c(x$info$labelSum, sortValue))
    }


    ## Prepare plot
    #trait in "" since it is not defined
    trait <- tmp$trait
    p <- qplot(x=by, y=trait, group=path, data=tmp, color=path, linetype=path, geom="line")

    p <- p + geom_line(size=lineSize)

    p <- p + xlab(label=ifelse(is.null(xlab), by,    xlab))
    p <- p + ylab(label=ifelse(is.null(ylab), lT[i], ylab[i])) #lT[i] is the TRAIT!!!

    if (!is.null(xlim)) {
      if (is.list(xlim)) {
        xlimI <- xlim[[i]]
      } else {
        xlimI <- xlim
      }
      p <- p + scale_x_continuous(limits=xlimI)
    }

    if (!is.null(ylimT)) {
      if (is.list(ylimT)) {
        ylimI <- ylimT[[i]]
      } else {
        ylimI <- ylimT
      }
      p <- p + scale_y_continuous(limits=ylimI)
    }

    if (useDirectLabels) p <- directlabels::direct.label(p=p, method=method)

    ## This needs to follow direct.label
    p <- p + scale_colour_manual(values=colorI,
    name=ifelse(is.null(labelPath), path, labelPath))

    p <- p + scale_linetype_manual(values=lineTypeI,
    name=ifelse(is.null(labelPath), path, labelPath))
    
    ret[[i]] <- p
  
}

  
  ## --- Return ---

  class(ret) <- c("plotSummaryAlphaPart", class(ret))
  ret

}


