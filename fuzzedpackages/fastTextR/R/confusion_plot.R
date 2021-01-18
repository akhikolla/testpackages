get_position <- function(x) {
    dat <- x[is.finite(x)]
    ind <- which(is.finite(x), arr.ind = TRUE)
    pos <- ind
    pos[,1] <- ind[,2]
    pos[,2] <- -ind[,1] + 1 + nrow(x)
    pos
}

set_default <- function(x, default) {
    if ( class(x) == class(default) ) {
        return(x)
    }
    return(default)
}

plot_confusion_matrix <- function(x, numbers = list(), label = list(), 
                                  rectangle.color = "grey") {
    row_names <- rownames(x)
    col_names <- colnames(x)

    number.digits <- set_default(numbers$digits, 2)
    number.cex <- set_default(numbers$cex, 0.8)
    number.font <- set_default(numbers$font, 2)
    number.color <- set_default(numbers$color, "red")

    label.cex <- set_default(label$cex, 1)
    label.color <- set_default(label$color, "black")
    label.offset <- 0.4
    label.rotation <- 90

    grid.color <- "grey"

    rectangle.color <- "#66CC00"
    rectangle.border.color <- "white"

    pos <- get_position(x)
    n1 <- min(pos[,2])
    n2 <- max(pos[,2])
    nn <- n2 - n1

    m1 <- min(pos[,1])
    m2 <- max(pos[,1])
    mm <- max(1, m2 - m1)

    plot.new()
    xlab_width <- ylab_width <- 0

    for (i in seq_len(50)) {
        xlim <- c(m1 - 0.5 - xlab_width, m2 + 0.5)
        ylim <- c(n1 - 0.5, n2 + 0.5 + ylab_width)
        plot.window(xlim + c(-0.2, 0.2), ylim + c(-0.2, 0.2), 
                    asp = 1, xaxs = "i", yaxs = "i")
        x.tmp <- max(strwidth(row_names, cex = label.cex))
        y.tmp <- max(strwidth(col_names, cex = label.cex))
        if (min(x.tmp - xlab_width, y.tmp - ylab_width) < 0.0001) {
            break
        }
        xlab_width <- x.tmp
        ylab_width <- y.tmp
    }

    plot.window(xlim = xlim , ylim = ylim,
                asp = 1, xlab = "", ylab = "", xaxs = "i", yaxs = "i")

    ## add rectangles
    rectangle <- matrix(0, nrow(pos), 2)
    colnames(rectangle) <- c("width", "height")
    rectangle[, "width"] <- as.vector(x / rowSums(x))
    rectangle[, "height"] <- as.vector(x / colSums(x))

    symbols(pos, add = TRUE, inches = FALSE,
            rectangles = rectangle, bg = rectangle.color, fg = rectangle.border.color)

    ## add grid
    symbols(pos, add = TRUE, inches = FALSE,  bg = NA,
            squares = rep(1, length(x)), fg = grid.color )

    if ( length(numbers) ) {
        ## add numbers
        i <- which(round(as.vector(x), number.digits) > 0)
        text(pos[i, 1], 0.2 + pos[i, 2], font = number.font, col = number.color,
             labels = round(as.vector(x), number.digits)[i], cex = number.cex)
    }

    ## add row and col labels
    label_offset <- strwidth("W", cex = label.cex) * label.offset
    pos.xlabel <- cbind(m1:m2, n2 + 0.5 + label_offset )
    pos.ylabel <- cbind(m1 - 0.5, n2:n1)

    text(pos.xlabel[,1], pos.xlabel[,2], col_names, srt = label.rotation,
         adj = ifelse(label.rotation == 0, c(0.5,0), c(0,0)),
         col = label.color, cex = label.cex, offset = label.offset)##, ...)
    text(pos.ylabel[,1], pos.ylabel[,2], row_names,
         col = label.color, cex = label.cex, pos = 2, offset = label.offset)##, ...)

    invisible(x)
}
