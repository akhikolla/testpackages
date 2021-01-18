plot.etm <- function(x, tr.choice, xlab = "Time", ylab = "Transition Probability",
                     col = 1, lty, xlim, ylim, conf.int = FALSE, level = 0.95,
                     ci.fun = "linear", ci.col = col, ci.lty = 3,
                     legend = TRUE, legend.pos, curvlab, legend.bty = "n", ...) {

    if (!inherits(x, "etm"))
        stop("'x' must be a 'etm' object")

    is_stratified <- !is.null(x$strata)

    ufrom <- unique(x$trans$from)
    uto <- unique(x$trans$to)
    absorb <- setdiff(uto, ufrom)
    nam1 <- dimnames(x$est)[[1]]
    nam2 <- dimnames(x$est)[[2]]
    pos <- c(paste(nam1[!(nam1 %in% as.character(absorb))],
                   nam2[!(nam2 %in% as.character(absorb))]),
             paste(x$trans$from, x$trans$to))
    if (missing(tr.choice)) {
        if (is_stratified) {
            tr.choice <- pos[1]
        } else {
            tr.choice <- pos
        }
    }

    ref <- sapply(1:length(x$state.names), function(i) {
        paste(x$state.names, x$state.names[i])
    })
    ref <- matrix(ref)
    if (sum(tr.choice %in% ref == FALSE) > 0)
        stop("Argument 'tr.choice' and possible transitions must match")

    if (is_stratified) {

        lstrat <- length(x$strata)
        temp <- lapply(seq_len(lstrat), function(i) {
            tmp <- ci.transfo(x[[i]], tr.choice, level, ci.fun)
            tmp2 <- lapply(tmp, cbind, strata = x$strata[[i]])
            tmp2
        })
        temp <- do.call(c, temp)
    } else {
        temp <- ci.transfo(x, tr.choice, level, ci.fun)
    }

    lt <- length(temp)

    if (missing(lty)) {
        lty <- seq_len(lt)
    }
    else if (length(lty) < lt) {
        lty <- lty * rep(1, lt)
    }
    if (length(col) < lt)
        col <- col * rep(1, lt)

    if (missing(xlim)) {
        xlim <- c(0, max(sapply(temp, function(x) max(x$time))))
    }
    if (missing(ylim)) {
        ylim <- c(0, 1)
    }

    graphics::plot(xlim, ylim, xlab = xlab, ylab = ylab,
                   xlim = xlim, ylim = ylim, type = "n", ...)

    for (i in seq_len(lt)) {
        graphics::lines(temp[[i]]$time, temp[[i]]$P, type = "s",
                        col = col[i], lty = lty[i], ...)
    }

    if (conf.int && !is.null(x$cov)) {
        if (length(ci.col) < lt)
            ci.col <- ci.col * rep(1, lt)
        if (length(ci.lty) < lt)
            ci.lty <- ci.lty * rep(1, lt)
        for (i in seq_len(lt)) {
            graphics::lines(temp[[i]]$time, temp[[i]]$lower, type = "s",
                            col = ci.col[i], lty = ci.lty[i], ...)
            graphics::lines(temp[[i]]$time, temp[[i]]$upper, type = "s",
                            col = ci.col[i], lty = ci.lty[i], ...)
        }
    }

    ## Extend  to deal with the strata
    if (legend) {
        if (missing(legend.pos))
            legend.pos <- "topleft"
        if (missing(curvlab)) {
            if (is_stratified) {
                curvlab <- as.vector(sapply(tr.choice, paste, x$strata))
            } else {
                curvlab <- tr.choice
            }
        }
        if (is.list(legend.pos)) legend.pos <- unlist(legend.pos)
        if (length(legend.pos) == 1) {
            xx <- legend.pos
            yy <- NULL
        }
        if (length(legend.pos) == 2) {
            xx <- legend.pos[1]
            yy <- legend.pos[2]
        }
        args <- list(...)
        ii <- pmatch(names(args),
                     names(formals("legend")[-charmatch("bty",names(formals("legend")))]))
        do.call("legend", c(list(xx, yy, curvlab, col=col, lty=lty, bty = legend.bty),
                            args[!is.na(ii)]))
    }

    invisible()
}
