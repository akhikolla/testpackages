#####################################################################
### Pseudo values for excess LoS
### Arthur Allignol <arthur.allignol@uni-ulm.de>
#####################################################################


closPseudo <- function(data, state.names, tra, cens.name, s = 0,
                       formula, na.action,
                       aw = FALSE, ratio = FALSE, ncores = 1,
                       trick_ties= FALSE) {
    
    stopifnot("data.frame" %in% class(data))
    data <- data.table(data)
    
    ## take care of the formula argument
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    temp <- c("", "formula", "data", "id", "subset", "na.action")
    m <- m[match(temp, names(m), nomatch = 0)]
    Terms <- if (missing(data)) terms(formula)
             else terms(formula, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())

    m <- data.table(cbind(id = data$id, m))
    
    n <- length(unique(data[, id]))
    
### get a minimal data set for computing the pseudo values
    reg <- names(data)
    names_msm <- intersect(c("id", "entry", "exit", "time", "from", "to"), reg)
    dat_clos <- data[, names_msm, with = FALSE]
    
    ## theta. From there we will see what kind of model it is
    ## is no alternative weights, NULL
    ## No competing risks: not in the list
    theta <- unlist(etm::clos(etm::etm(dat_clos, state.names = state.names, tra = tra,
                                       cens.name = cens.name, s = 0, covariance = FALSE),
                              aw = aw, ratio = ratio)[c("e.phi", "e.phi.weights.1",
                                                        "e.phi.weights.other",
                                                        "e.phi2", "e.phi3")])
    
    competing <- "e.phi2" %in% names(theta)
    
    ## Compute pseudo values, and store results depending of competing
    ## and aw
    namen <- c("ps.e.phi", "ps.e.phi.weights.1", "ps.e.phi.weights.other",
               "ps.e.phi2", "ps.e.phi3")
    
    if (trick_ties) {
        
        ## we want to find all patients that have the same "dynamic"
        ## and get id's of some of them to compute PS
        make_cat <- function(entry, exit, from, to) {
            if (length(from) == 1) {
                cat <- paste(entry, exit, from, to, sep = "_")
            } else {
                cat <- paste(entry[1], exit[1], from[1], to[1],
                             entry[2], exit[2], from[2], to[2],
                             sep = "_")
            }
            
            list(categs = cat)
        }
        
        cat_dyn <- dat_clos[, make_cat(entry, exit, from, to), by = "id"]
        cat_dyn_red <- unique(data.table(cat_dyn, key = "categs"))
        ids <- cat_dyn_red[, id]
        
    } else {
        ids <- unique(data$id)
    }
    
    psMatrix <- parallel::mclapply(seq_along(ids), function(i) {
        temp <- clos(etm(dat_clos[!(id %in% ids[i])],
                         state.names = state.names, tra = tra,
                         cens.name = cens.name, s = 0, covariance = FALSE),
                     aw = aw, ratio = ratio)
        
        data.table(cbind(temp$e.phi, temp$e.phi.weights.1, temp$e.phi.weights.other,
                         temp$e.phi2, temp$e.phi3))
    }, mc.cores = ncores)

    psMatrix <- rbindlist(psMatrix)
    
    psMatrix <- lapply(seq_along(psMatrix), function(i) {
        n * theta[i] - (n - 1) * psMatrix[, i, with = FALSE]
    })
    psMatrix <- do.call(cbind, psMatrix)
    setnames(psMatrix, namen[c(TRUE, aw, aw, competing, competing)])
    
    ## if trick, we need to merge intelligently
    if (trick_ties) {
        
        bouh <- cbind(cat_dyn_red, psMatrix)
        setkeyv(cat_dyn, "categs")
        psMatrix <- merge(bouh, cat_dyn, by = "categs", all.y = TRUE)
        psMatrix <- psMatrix[, c("id.y", namen[c(TRUE, aw, aw, competing, competing)]),
                             with = FALSE]
        setnames(psMatrix, c("id", colnames(psMatrix)[-1]))
        setkeyv(psMatrix, "id")
    } else {
        psMatrix <- data.frame(ids, psMatrix, stringsAsFactors=TRUE)
        names(psMatrix) <- c("id", names(psMatrix)[-1])
    }
    
    cov <- unique(data.table(m, key = "id"))
    
    theta <- matrix(theta, nrow = 1)
    colnames(theta) <- c("e.phi", "e.phi.weights.1",
                         "e.phi.weights.other", "e.phi2",
                         "e.phi3")[c(TRUE, aw, aw, competing, competing)]
    
    pseudoData <- merge(psMatrix, cov)
    
    zzz <- list(pseudoData = pseudoData,
                theta = theta, aw = aw, call = call)
    class(zzz) <- "closPseudo"
    
    zzz
}
