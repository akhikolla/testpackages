trprob <- function(x, ...) {
    UseMethod("trprob")
}

trcov <- function(x, ...) {
    UseMethod("trcov")
}

trprob.etm <- function(x, tr.choice, timepoints, ...) {

    if (!inherits(x, "etm"))
        stop("'x' must be a 'etm' object")
    if (!is.character(tr.choice))
        stop("'tr.choice' must be a character vector")
    if (length(tr.choice) != 1)
        stop("The function only extracts 1 transition probability")

    pos <- sapply(1:length(x$state.names), function(i) {
        paste(x$state.names, x$state.names[i])
    })
    pos <- matrix(pos)
    if (!(tr.choice %in% pos))
        stop("'tr.choice' not in the possible transitions")

    trans.sep <- strsplit(tr.choice, " ")
    if (length(trans.sep[[1]]) != 2) {
        tt <- charmatch(trans.sep[[1]], x$state.names, nomatch = 0)
        trans.sep[[1]] <- x$state.names[tt]
    }
    trans.sep <- unlist(trans.sep)

    miss_timepoints <- missing(timepoints)

    ## Number of strata. Will be computed in this if condition
    if (!is.null(x$strata_variable)) {
        ns <- length(x$strata)
    } else {
        if (miss_timepoints) {
            tmp <- x$est[trans.sep[1], trans.sep[2], ]
        } else {
            ind <- findInterval(timepoints, x$time)
            tmp <- numeric(length(timepoints))
            place <- which(ind != 0)
            noplace <- which(ind == 0)
            tmp[place] <- x$est[trans.sep[1], trans.sep[2], ind]
            if (trans.sep[1] == trans.sep[2])
                tmp[noplace] <- 1
        }
        return(tmp)
    }

    ## Compute for the case with strata
    res <- lapply(1:ns, function(i) {
        if (miss_timepoints) {
            tmp <- x[[i]]$est[trans.sep[1], trans.sep[2], ]
        } else {
            ind <- findInterval(timepoints, x[[i]]$time)
            tmp <- numeric(length(timepoints))
            place <- which(ind != 0)
            noplace <- which(ind == 0)
            tmp[place] <- x[[i]]$est[trans.sep[1], trans.sep[2], ind]
            if (trans.sep[1] == trans.sep[2])
                tmp[noplace] <- 1
        }
        tmp
    })

    names(res) <- names(x)[seq_len(ns)]

    res

}

trcov.etm <- function(x, tr.choice, timepoints, ...) {
    if (!inherits(x, "etm"))
        stop("'x' must be a 'etm' object")
    if (!is.character(tr.choice))
        stop("'tr.choice' must be a character vector")
    if (!(length(tr.choice) %in% c(1, 2)))
        stop("'tr.choice' must be of length 1 or 2")
    pos <- sapply(1:length(x$state.names), function(i) {
        paste(x$state.names, x$state.names[i])
    })

    pos <- matrix(pos)
    if (!all((tr.choice %in% pos)))
        stop("'tr.choice' not in the possible transitions")
    if (length(tr.choice) == 1) {
        tr.choice <- rep(tr.choice, 2)
    }

    miss_timepoints <- missing(timepoints)

    if (!is.null(x$strata_variable)) {
        ns <- length(x$strata)
        if (is.null(x[[1]]$cov)) stop("The covariance matrix was not computed")
    } else {
        if (is.null(x$cov)) stop("The covariance matrix was not computed")
        if (miss_timepoints) {
        tmp <- x$cov[tr.choice[1], tr.choice[2], ]
        }
        else {
            ind <- findInterval(timepoints, x$time)
            tmp <- numeric(length(timepoints))
            place <- which(ind != 0)
            tmp[place] <- x$cov[tr.choice[1], tr.choice[2], ind]
        }
        return(tmp)
    }

    res <- lapply(1:ns, function(i) {
        if (miss_timepoints) {
            tmp <- x[[i]]$cov[tr.choice[1], tr.choice[2], ]
        }
        else {
            ind <- findInterval(timepoints, x[[i]]$time)
            tmp <- numeric(length(timepoints))
            place <- which(ind != 0)
            tmp[place] <- x[[i]]$cov[tr.choice[1], tr.choice[2], ind]
        }
        tmp
    })

    names(res) <- names(x)[seq_len(ns)]

    res

}


##############################
### For the stratified etm ###
##############################

"[.etm" <- function(x, ..., drop = FALSE) {

    if (missing(..1)) i <- NULL else i <- ..1

    ## No subscript, do nothing
    if (is.null(i)) return(x)

    ## No strata, do nothing
    if (is.null(x$strata_variable)) return(x)

    if (is.character(i)) {
        ind <- match(gsub(" ", "", i, fixed = TRUE),
                     gsub(" ", "", x$strata, fixed = TRUE))
        if (any(is.na(ind))) stop(paste("subscript(s)",
                                        paste(i[is.na(ind)], sep = " "),
                                        "not matched"))
    } else {
        ind <- i
        if (max(ind) > length(x$strata))
            stop(paste0("There is only ", length(x$strata), " strata"))
    }

    if (length(ind) == 1) {

        res <- x[[ind]]
        res$trans <- x$trans
        res$tra <- x$tra
        res$state.names <- x$state.names
        res$data <- x$data
        res$strat_variable <- NULL
        res$strata <- NULL

    } else {

        res <- unclass(x)[ind]
        res$trans <- x$trans
        res$tra <- x$tra
        res$state.names <- x$state.names
        res$data <- x$data
        res$strat_variable <- NULL
        res$strata <- NULL

    }

    class(res) <- "etm"
    res

}
