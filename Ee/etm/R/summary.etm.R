find_times <- function(d, timepoints) {

    ind <- findInterval(timepoints, d$time)
    ind0 <- sum(ind == 0)

    dd <- d[ind, ]

    if (ind0 > 0) {
        tmp <- d[1, , drop = FALSE]
        tmp$P <- round(tmp$P)
        tmp$var <- 0
        tmp$lower <- tmp$upper <- tmp$P
        tmp$n.event <- 0

        for (i in seq_len(ind0)) dd <- rbind(tmp, dd)
    }

    dd$time <- timepoints
    dd$n.event <- cumsum(dd$n.event)

    dd
}



summary.etm <- function(object, tr.choice, ci.fun = "linear", level = 0.95, times, ...) {

    if (!inherits(object, "etm"))
        stop("'object' must be of class 'etm'")

    if (level <= 0 | level > 1) {
        stop ("'level' must be between 0 and 1")
    }

    ref <- c("linear", "log", "cloglog", "log-log")
    if (sum(ci.fun %in% ref == FALSE) != 0) {
        stop("'ci.fun' is not correct. See help page")
    }

    ## Number of strata. Will be computed in this if condition
    ns <- 1
    if (!is.null(object$strata_variable)) {
        ns <- length(object$strata)
        time <- unique(sapply(1:ns, function(i) {
            object[[i]]$time
        }))
    } else {
        time <- object$time
    }

    ## If no event time between s and t, don't need a summary
    if (is.null(time)) stop("no event time")

    ## Derive the transition names we need
    if (missing(tr.choice)) {
        if (!is.null(object$strata_variable)) {

            indi <- lapply(1:ns, function(i) {
                !apply(object[[i]]$est != 0, c(1, 2), function(temp){all(temp == FALSE)})
            })
            indi <- do.call("+", indi) > 0

        } else {

            ind <- object$est != 0
            indi <- !apply(ind, c(1, 2), function(temp){all(temp == FALSE)})
        }

        tmp <- which(indi, arr.ind = TRUE)
        tmp <- tmp[order(tmp[, 1]), ]
        namen <- list(rownames(indi), colnames(indi))
        trs <- lapply(seq_len(NROW(tmp)), function(i) {
            paste(namen[[1]][tmp[i, 1]], namen[[2]][tmp[i, 2]], sep = " ")
        })
        trs <- cbind(trs)
        absorb <- setdiff(as.character(object$tran$to), as.character(object$trans$from))
        for (i in seq_along(absorb))
            trs <- trs[-grep(paste("^", absorb[i], sep =""), trs, perl = TRUE)]

    } else {

        ref <- sapply(1:length(object$state.names), function(i) {
            paste(object$state.names, object$state.names[i])
        })
        ref <- matrix(ref)
        if (sum(tr.choice %in% ref == FALSE) > 0)
            stop("Argument 'tr.choice' and possible transitions must match")
        trs <- tr.choice
    }

    not_missing <- !missing(times)
    if (ns > 1) {

        res <- lapply(seq_len(ns), function(i) {
            tmp <- ci.transfo(object[[i]], trs, level, ci.fun)
            if (not_missing) tmp <- lapply(tmp, find_times, timepoints = times)
            class(tmp) <- "summary.etm"
            tmp
        })
        names(res) <- object$strata
        class(res) <- "summary.etm"

    } else {

        res <- ci.transfo(object, trs, level, ci.fun)
        if (not_missing) res <- lapply(res, find_times, timepoints = times)
        class(res) <- "summary.etm"
    }

    res
}


