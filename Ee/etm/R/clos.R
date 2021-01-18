### Expected change of LoS
### Arthur Allignol <arthur.allignol@uni-ulm.de>
#####################################################################


clos <- function(x, aw, ratio, ...) {
    UseMethod("clos")
}

clos.etm <- function(x, aw = FALSE, ratio = FALSE, ...) {
    if (!inherits(x, "etm")) {
        stop("'x' must be an 'etm' object")
    }
    if (is.null(x$delta.na)) {
        stop("Needs the increment of the Nelson-Aalen estimator")
    }

    if (!is.null(x$strata)) stop("'clos' is not yet implemented for etm objects with strata. Please use e.g., 'clos(etm_object_name[1]'")

    ## test if we have an illness-death model
    dims <- dim(x$est)
    comp.risk <- FALSE
    if (dims[1] == 4) comp.risk <- TRUE
    ## I <- diag(1, dims[1])
    ## tr.mat <- array(apply(x$delta.na, 3, "+", I), dim = dims)
    if (comp.risk) {
        res <- clos.cp(x, aw, ratio)
        ## stop("not yet")
    }
    else res <- clos.nocp(x, aw, ratio)

    class(res) <- "clos.etm"
    res
}


clos.msfit <- function(x, aw = FALSE, ratio = FALSE, cox_model, ...) {

    if (!inherits(x, "msfit")) {
        stop("'x' must be an 'msfit' object")
    }

    if (missing(cox_model)) {
        stop("cox model fit is missing")
    }

    ## CumHaz <- dplyr::group_by(x$Haz, trans)
    ## CumHaz <- dplyr::mutate(CumHaz, dhaz = diff(c(0, Haz)))
    CumHaz <- data.table(x$Haz)
    CumHaz[, dhaz := diff(c(0, Haz)), by = trans]

    trans <- x$trans[!is.na(x$trans)]
    ltrans <- dim(x$trans)
    ltimes <- unique(CumHaz[, length(time), by = trans][, V1])

    times <- sort(unique(CumHaz$time))


    dna <- nev <- array(0, dim = c(ltrans[1], ltrans[1], ltimes))
    nrisk <- matrix(0, nrow = ltimes, ncol = ltrans[1])

    ## Take care of the cox model and do the transformations in the
    ## same loop
    temp_surv <- survival::survfit(cox_model)
    dat_surv <- data.frame(time = temp_surv$time,
                           n.risk = temp_surv$n.risk,
                           n.event = temp_surv$n.event,
                           trans = rep(trans, temp_surv$strata), stringsAsFactors=TRUE)

    for (i in trans) {
        aa <- which(x$trans == i, arr.ind = TRUE)
        dna[aa[1], aa[2], ] <- CumHaz$dhaz[CumHaz$trans == i]

        ## fill nev and nrisk
        dat_temp <- dat_surv[dat_surv$trans == i, ]
        ind <- findInterval(times, dat_temp$time)
        place <- which(ind != 0)
        tmp <- integer(ltimes)
        tmp <- cumsum(dat_temp$n.event)[ind]

        nev[aa[1], aa[2], place] <- c(tmp[1], diff(tmp))
        nrisk[place, aa[1]] <- dat_temp$n.risk[ind]

        tt <- nrisk[2:ltimes, aa[1]] - nev[aa[1], aa[2], 1:(ltimes-1)]

        if (!all((tt == 0) == FALSE)) {
            uu <- max(times[tt == 0])

            if (uu < max(times)) {
                vv <- which(times == uu)
                nrisk[(vv + 1):nrow(nrisk), aa[1]] <- 0
            }
        }
    }

    ## Ugly risk set fix for illness-death models where nobody starts
    ## in state 1 at time 0
    ind <- which(nrisk[, 2] != 0)[1]
    dni <- which(nev[1, 2, ] != 0)[1]

    if ((ind - 1) > dni) {
        nrisk[(dni + 1):(ind - 1), 2] <- cumsum(nev[1, 2, dni:(ind - 2)])
    }

    ii <- seq_len(ltrans[1])
    for (i in seq_along(times)) {
        dna[cbind(ii, ii, i)] <- -(.rowSums(dna[, , i], ltrans[1], ltrans[1], FALSE))
    }

    ## Need to compute AJ (again)
    which.compute <- rep(1, ltimes)
    first <- 1; last <- ltimes
    est <- prodint(dna, times, first, last, which.compute)

    comp.risk <- FALSE
    if (ltrans[1] == 4) comp.risk <- TRUE
    ## I <- diag(1, ltrans[1])
    ## tr.mat <- array(apply(dna, 3, "+", I), dim = c(ltrans, ltimes))

    d_tmp <- data.frame(exit = cox_model$y[, 2])

    zzz <- list(est = est$est,
                delta.na = dna,
                time = est$time,
                n.event = nev,
                n.risk = nrisk[, 1:2],            ## dirty. but these should be illness-death models
                data = d_tmp
                )


    if (comp.risk) {
        stop("'clos.msfit' is not yet implemented with competing risks")
    }
    else res <- clos.nocp(zzz, aw, ratio)
    class(res) <- "clos.etm"

    res
}

