## plot_cross.R

#  CCcolors
#
#' Collaborative Cross colors
#'
#' Get the vector of colors for the Collaborative Cross
#'
#' @param palette Which version of the colors to use? (New or original)
#'
#' @return vector of eight colors
#'
#' @keywords color
#' @export
#' @importFrom grDevices rgb
#' @importFrom graphics par points segments
#' @examples
#' CCcolors()
CCcolors <-
    function(palette=c("new", "original", "official")) {
        palette <- match.arg(palette)

        if(palette=="new") {
            return(
                c(AJ   ="#FFDC00",
                  B6   ="#888888",
                  "129"="#F08080",
                  NOD  ="#0064C9",
                  NZO  ="#7FDBFF",
                  CAST ="#2ECC40",
                  PWK  ="#FF4136",
                  WSB  ="#B10DC9")
            )
        }
        else {
            return(
                c("AJ"  =rgb(240,240,  0,maxColorValue=255),
                  "B6"  =rgb(128,128,128,maxColorValue=255),
                  "129" =rgb(240,128,128,maxColorValue=255),
                  "NOD" =rgb( 16, 16,240,maxColorValue=255),
                  "NZO" =rgb(  0,160,240,maxColorValue=255),
                  "CAST"=rgb(  0,160,  0,maxColorValue=255),
                  "PWK" =rgb(240,  0,  0,maxColorValue=255),
                  "WSB" =rgb(144,  0,224,maxColorValue=255))
            )
        }
    }

#  plot_ind
#
#' Plot an individual
#'
#' Add an individual, as a pair of chromosomes, to a plot
#'
#' @param ind An individual object, as output by
#' [create_parent()] or [cross()]
#' @param center (x,y) vector for the center of the individual
#' @param chrlength Length of chromosomes (Can be a vector of length
#' 2, in which case the two chromosomes will be different lengths,
#' aligned at the top. This is for the X chromosome.)
#' @param chrwidth Width of chromosomes
#' @param gap Gap between chromosomes
#' @param col Vector of colors
#' @param border Color for border
#' @param lend Passed to [graphics::rect()]
#' @param ljoin Passed to [graphics::rect()]
#' @param allborders If TRUE, put borders around all segments
#' @param ... Additional arguments passed to rect()
#'
#' @return None.
#'
#' @importFrom graphics rect
#' @keywords hplot
#' @export
#' @seealso [plot_crosslines()]
#' @examples
#' \dontshow{set.seed(67452378)}
#' mom <- create_parent(100, 1:2)
#' dad <- create_parent(100, 3:4)
#' kid <- cross(mom, dad)
#' plot(0,0, type="n", xlim=c(0, 100), ylim=c(0,100),
#'      xaxt="n", yaxt="n", xlab="", ylab="")
#' loc <- list(c(25,75), c(75,75), c(50,25))
#' plot_ind(mom, loc[[1]])
#' plot_ind(dad, loc[[2]])
#' plot_ind(kid, loc[[3]])
#' plot_crosslines(loc[[1]], loc[[2]], loc[[3]])
plot_ind <-
    function(ind, center, chrlength=30, chrwidth=3, gap=3, col=CCcolors(),
             border="black", lend=1, ljoin=1, allborders=FALSE, ...)
{
    max_alleles <- max(sapply(ind, function(a) max(a$alleles)))
    if(max_alleles > length(col))
        stop("Need more colors: length(col)=", length(col), " but max allele = ", max_alleles)

    if(length(chrlength)==0) stop("chrlength has length 0")
    if(length(chrlength)==1) chrlength <- rep(chrlength, 2)
    if(length(chrlength)>2) {
        warning("chrlength should have length 1 or 2; the first two values will be used")
        chrlength <- chrlength[1:2]
    }

    for(i in 1:2) {
        sgn <- i*2-3
        chr <- ind[[i]]
        pos <- c(0, ind[[i]]$locations)
        allele <- ind[[i]]$alleles

        # rescale pos
        top <- center[2]+max(chrlength)/2
        bottom <- top - chrlength[i]
        if(diff(par("usr")[3:4]) < 0) { # small values at top of figure
            tmp <- top
            top <- bottom
            bottom <- top
        }
        start <- pos[1]
        end <- pos[length(pos)]
        pos <- (pos-start)*(bottom-top)/(end-start) + top

        left <- center[1] + sgn*(gap/2+chrwidth)
        right <- center[1] + sgn*gap/2

        internalborder <- ifelse(allborders, border, NA)
        for(j in 2:length(pos))
            rect(rep(left, length(pos)-1),  pos[-length(pos)],
                 rep(right, length(pos)-1), pos[-1],
                 col=col[allele],
                 border=internalborder, lend=lend, ljoin=ljoin, ...)

        if(!is.na(border) && !is.null(border)) # draw border
            rect(center[1] + sgn*(gap/2+chrwidth), top,
                 center[1] + sgn*gap/2,         bottom,
                 col=NA, border=border, lend=lend, ljoin=ljoin, ...)
    }
    invisible(NULL)
}

#  plot_crosslines
#
#' Plot cross lines
#'
#' Add lines for a cross
#'
#' @param momloc An (x,y) vector with center location for mother
#' @param dadloc An (x,y) vector with center location for mother
#' @param kidsloc Either an (x,y) vector with center location for a
#' kid, or a list of such for multiple kids
#' @param gap Gap arrows and points/rectangles
#' @param cex Character expansion for x point
#' @param chrlength Length of chromosomes
#' @param lwd Line width for points, segments, and arrows
#' @param arrow_length The `length` parameter in the call to
#' [graphics::arrows()]
#' @param col Color of lines and points
#' @param ... Additional arguments passed to arrows() and segments()
#'
#' @return None.
#'
#' @importFrom graphics points arrows
#' @keywords hplot
#' @export
#' @seealso [plot_ind()]
#' @examples
#' \dontshow{set.seed(67452378)}
#' mom <- create_parent(100, 1:2)
#' dad <- create_parent(100, 3:4)
#' kids <- lapply(1:4, function(junk) cross(mom, dad))
#' plot(0,0, type="n", xlim=c(0, 100), ylim=c(0,100),
#'      xaxt="n", yaxt="n", xlab="", ylab="")
#' loc <- list(c(25,75), c(75,75), c(12.5,25), c(37.5,25), c(62.5, 25), c(87.5,25))
#' plot_ind(mom, loc[[1]])
#' plot_ind(dad, loc[[2]])
#' for(i in 1:4) plot_ind(kids[[i]], loc[[i+2]])
#' plot_crosslines(loc[[1]], loc[[2]], loc[3:6])
plot_crosslines <-
    function(momloc, dadloc, kidsloc, gap=3, chrlength=30, cex=1.5,
             lwd=2, arrow_length=0.1, col="white", ...)
{
    stopifnot(length(momloc)==2, length(dadloc)==2)

    point <- colMeans(rbind(momloc, dadloc))
    points(point[1], point[2], pch=4, cex=cex, lwd=2, col=col, ...)

    if(!is.list(kidsloc)) { # 1 kid
        stopifnot(length(kidsloc)==2)
        sgn <- sign(kidsloc[2] - point[2])
        arrows(point[1], point[2]+sgn*gap*2, kidsloc[1], kidsloc[2]-sgn*(chrlength/2+gap),
               lwd=lwd, length=arrow_length, col=col, ...)
    } else { # multiple kids
        if(any(vapply(kidsloc, length, 2) != 2))
            stop("kidsloc must be a list of vectors of length 2")

        kidx <- vapply(kidsloc, "[", 0, 1)
        kidy <- vapply(kidsloc, "[", 0, 2)

        ave <- c(mean(kidx), mean(kidy))
        midpt <- c((point[1]+ave[1])/2, (point[2]+ave[2])/2)
        sgn <- sign(ave[2] - point[2])

        segments(point[1], point[2]+sgn*gap, point[1], midpt[2],
                 lwd=lwd, col=col, ...)

        segments(min(kidx), midpt[2], max(kidx), midpt[2],
                 lwd=lwd, col=col, ...)

        arrows(kidx, rep(midpt[2], length(kidx)),
               kidx, kidy-sgn*(chrlength/2+gap),
               lwd=lwd, length=arrow_length, col=col, ...)
    }

    invisible(NULL)
}
