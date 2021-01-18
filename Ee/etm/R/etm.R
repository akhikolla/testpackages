etm <- function(data, ...) {
    UseMethod("etm")
}


etm.data.frame <- function(data, state.names, tra, cens.name, s, t = "last",
                           covariance = TRUE, delta.na = TRUE, modif = FALSE,
                           c = 1, alpha = NULL, strata, ...) {

    if (missing(data))
        stop("Argument 'data' is missing with no default")
    x <- data.table(data)
    if (missing(tra))
        stop("Argument 'tra' is missing with no default")
    if (missing(state.names))
        stop("Argument 'state.names' is missing with no default")
    if (missing(cens.name))
        stop("Argument 'cens.name' is missing with no default")
    if (missing(s))
        stop("Argument 's' is missing with no default")
    if (!is.data.frame(x))
        stop("Argument 'x' must be a data.frame")
    if (!(xor(sum(c("id", "from", "to", "time") %in% names(x)) != 4,
              sum(c("id", "from", "to", "entry", "exit") %in% names(x)) != 5)))
        stop("'data' must contain the right variables")
    if (nrow(tra) != ncol(tra))
        stop("Argument 'tra' must be a quadratic  matrix.")
    if (sum(diag(tra)) > 0)
        stop("transitions into the same state are not allowed")
    if (nrow(tra) != length(state.names)) {
        stop("The row number of 'tra' must be equal to the number of states.")
    }
    if (!is.logical(tra)) {
        stop("'tra' must be a matrix of logical values, which describes the possible transitions.")
    }
    if (length(state.names) != length(unique(state.names))) {
        stop("The state names must be unique.")
    }
    if (!(is.null(cens.name))) {
        if (cens.name %in% state.names) {
            stop("The name of the censoring variable just is a name of the model states.")
        }
    }
    ## The stratification variable
    if (missing(strata)) {
        strat_var <- "X"
        x$X <- "1"
        is_stratified <- FALSE
    } else {
        if (!all(strata %in% names(x)))
            stop("Stratification variables not in data")
        strat_var <- strata
        x[, (strat_var) := lapply(.SD, as.character), .SDcols = strat_var]
        is_stratified <- TRUE
    }

    ## if modif TRUE, check that the model is competing risks. else
    ## set to false and issue a warning
    if (modif == TRUE && covariance == TRUE) {
        ## check for competing risks
        tr.cp <- tra_comp(length(state.names) - 1)
        if (any(dim(tra) != dim(tr.cp)) | (all(dim(tra) == dim(tr.cp)) && !all(tra == tr.cp))) {
            covariance <- FALSE
            warning("The variance of the estimator with the Lay and Ying transformation is only computed for competing risks data")
        }
    }

    ## keep only the variables we need
    reg <- names(x)
    names_msm <- intersect(c("id", "entry", "exit", "time", "from", "to", strat_var), reg)
    x <- x[, names_msm, with = FALSE]

    ## Work through the stratification variables
    combi <- unique(x[, strat_var, with = FALSE])
    if (length(strat_var) == 1) {
        conditions <- lapply(seq_len(nrow(combi)), function(i) {
                                 parse(text = paste0(strat_var, " ==  '", combi[i], "'"))
                             })
    } else {
        conditions <- lapply(seq_len(nrow(combi)), function(i) {
                                 parse(text = paste(sapply(strat_var,
                                       function(j) {
                                           paste0(j, "== '", combi[i, j, with = FALSE], "'")
                                       }),
                                       collapse = " & "))
                             })
    }

    ## transitions
    colnames(tra) <- rownames(tra) <- state.names
    t.from <- lapply(1:dim(tra)[2], function(i) {
        rep(rownames(tra)[i], sum(tra[i, ]))
    })
    t.from <- unlist(t.from)
    t.to <- lapply(1:dim(tra)[2], function(i) {
        colnames(tra)[tra[i, ]==TRUE]
    })
    t.to <- unlist(t.to)
    trans <- data.frame(from=t.from, to=t.to, stringsAsFactors=TRUE)
    namen <- paste(trans[, 1], trans[, 2])

    ## test on transitions
    test <- x[, unique(paste(from, to))]
    if (!(is.null(cens.name))) {
        ref <- c(paste(trans$from, trans$to), paste(unique(trans$from), cens.name))
    } else {
        ref <- paste(trans$from, trans$to)
    }
    ref.wo.cens <- paste(trans$from, trans$to)
    if (!(all(test %in% ref)==TRUE))
        stop("There is undefined transitions in the data set")
    if (x[, sum(as.character(from) == as.character(to))] > 0)
        stop("Transitions into the same state are not allowed")
    if (!(all(ref.wo.cens %in% test) == TRUE))
        warning("You may have specified more possible transitions than actually present in the data")

    n <- length(unique(x$id))

### data.table transformation
    x[, id := if (is.character(id)) as.factor(id) else id]
    if (!(is.null(cens.name))) {
        x[, from := factor(from, levels = c(cens.name, state.names), ordered = TRUE)]
        levels(x$from) <- 0:length(state.names)
        x[, to := factor(to, levels = c(cens.name, state.names), ordered = TRUE)]
        levels(x$to) <- 0:length(state.names)
    } else{
        x[, from := factor(from, levels = state.names, ordered = TRUE)]
        levels(x$from) <- 1:length(state.names)
        x[, to := factor(x$to, levels = state.names, ordered = TRUE)]
        levels(x$to) <- 1:length(state.names)
    }

    ## if not, put like counting process data
    if ("time" %in% names(x)) {
        setorder(x, id, time)
        x[, idd := as.integer(id)]
        x[, masque := rbind(1, apply(as.matrix(idd), 2, diff))]
        x[, entree := c(0, time[1:(length(time) - 1)]) * (masque == 0)]
        x[, ':='(entry = entree,
                 exit = time,
                 entree = NULL,
                 time = NULL,
                 masque = NULL)]
        if (sum(x$entry < x$exit) != nrow(x))
            stop("Exit time from a state must be > entry time")
    } else {
        if (sum(x$entry < x$exit) != nrow(x))
            stop("Exit time from a state must be > entry time")
    }

    x[, from := as.integer(as.character(from))]
    x[, to := as.integer(as.character(to))]

    if (t=="last") t <- max(x$exit)
    if (!(0 <= s & s < t))
        stop("'s' and 't' must be positive, and s < t")
    if (t < x[, min(exit)] | s >= x[, max(exit)])
        stop("'s' or 't' is an invalid time")

    ## remove the lines in which transition before s
    x <- x[exit > s]
    ## remove the entries after t
    x <- x[entry < t]

    ## The Lai and Ying modification (if any)
    if (modif) {
        if (is.null(alpha)) {
            if (length(c) == 1) {
                c_modif <- c
            } else {
                if (length(c) != length(state.names)) {
                    stop("if specifying a unique c for each transient state, 'c' should be the same length as the 'state.names'")
                } else {
                    c_modif <- c
                }
            }
        } else {
            ## the original lay and ying proposal
            c_modif <- c * n^alpha
        }
    } else {
        c_modif <- 0
    }

    res <- lapply(conditions, function(ll)
    {
        zzz <- .etm(entry = x[eval(ll), entry],
                    exit = x[eval(ll), exit],
                    from = x[eval(ll), from],
                    to = x[eval(ll), to],
                    nstate = dim(tra)[1],
                    s,
                    t,
                    covariance,
                    c_modif)

        nrisk <- zzz$n.risk

        est <- zzz$est
        nev <- zzz$n.event
        var_aj <- zzz$cov

        dimnames(est) <- list(state.names, state.names, zzz$time)

        if (!is.null(zzz$n.risk)) {
            colnames(nrisk) <- state.names
            nrisk <- nrisk[, !(colnames(nrisk) %in%
                               setdiff(unique(trans$to), unique(trans$from))),
                           drop = FALSE]
            dimnames(est) <- list(state.names, state.names, zzz$time)
            dimnames(nev) <- list(state.names, state.names, zzz$time)

            if (covariance) {
                pos <- sapply(1:length(state.names), function(i) {
                    paste(state.names, state.names[i])
                })
                pos <- matrix(pos)
                dimnames(var_aj) <- list(pos, pos, zzz$time)
                var_aj[var_aj < 0] <- 0
            } else {
                var_aj <- NULL
            }
        } else {
            var_aj <- NULL
        }

        res <- list(est = est,
                    cov = var_aj,
                    time = zzz$time,
                    n.risk = nrisk,
                    n.event = nev,
                    delta.na = zzz$dna,
                    s = s,
                    t = t)
        res
    })

    if (is_stratified) {
        names(res) <- do.call('c', lapply(conditions, as.character))
        res$trans <- trans
        res$tra <- tra
        res$state.names <- state.names
        res$data <- x
        res$strata_variable <- strata
        res$strata <- do.call('c', lapply(conditions, as.character))
        class(res) <- "etm"
    } else {
        res <- res[[1]]
        res$trans <- trans
        res$tra <- tra
        res$state.names <- state.names
        res$data <- x
        class(res) <- "etm"
    }

    res
}



### with Hist
## etm.formula <- function(formula, data, subset, na.action) {

##     call <- match.call()
##     if (missing(data))
##         data <- environment(formula)
##     if (!missing(subset))
##         data <- subset(data, subset = subset)

##     mod <- EventHistory.frame(formula,
##                               data,
##                               specials = c("strata", "factor"),
##                               stripSpecials = "strata",
##                               stripAlias = list(strata = c("Strata", "factor")),
##                               check.formula = TRUE)

##     mod
## }

