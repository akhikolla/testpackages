clos_pure <- function(x, aw = FALSE) {
    if (!inherits(x, "etm")) {
        stop("'x' must be an 'etm' object")
    }
    if (is.null(x$delta.na)) {
        stop("Needs the increment of the Nelson-Aalen estimator")
    }
    absorb <- setdiff(levels(x$trans$to), levels(x$trans$from))
    transient <- unique(x$state.names[!(x$state.names %in% absorb)])
    if (!(length(transient) == 2 && length(absorb) %in% c(1, 2)))
        stop("The multistate model must have 2 transient states \n and 1 or 2 absorbing states")
    dims <- dim(x$est)
    comp.risk <- FALSE
    if (dims[1] == 4) comp.risk <- TRUE
    I <- diag(1, dims[1])
    tr.mat <- array(apply(x$delta.na, 3, "+", I), dim = dims)
    #tr.mat <- x$delta.na
    if (comp.risk) {
        res <- clos.cp(x, tr.mat, aw)
    }
    else res <- clos.nocp(x, tr.mat, aw)
    class(res) <- "clos.etm"
    res
}

clos.cp_pure <- function(x, tr.mat, aw) {
    dims <- dim(x$est)
    times <- if (sum(x$n.event[, , dims[3]]) != 0) x$time else x$time[-length(x$time)]
    los <- matrix(rep(times, 3), ncol = 3, byrow = FALSE)
    phi2 <- phi3 <- matrix(0, nrow(los), 3)
    phi2[, 1] <- times; phi3[, 1] <- times
    lt <- length(times)
    tau <- max(times)
    los[length(times)-1, 2:3] <- rep(tau, 2)
    aj <- array(NA, c(4, 4, 1))
    aj[, , 1] <- diag(1, 4, 4)
    
    funn <- function(x, y) {
        x %*% aj[ , , y]
    }
    
    for (i in (lt - 2):1) {
        diffs <- diff(times[(i + 1):length(times)])
        mat <- tr.mat[ , , length(x$time[x$time <= times[i + 1]])]
        aj <- array(apply(X = diag(1:dim(aj)[3]), 1, funn, x = mat), c(4, 4, dim(aj)[3]))
        los[i, 3] <- times[i + 1] + matrix(diffs, nrow=1) %*%
            matrix(aj[2, 2, ], ncol = 1)
        los[i, 2] <- times[i + 1] + matrix(diffs, nrow=1) %*%
            matrix((aj[1, 1, ] + aj[1, 2, ]),ncol=1)
        phi2[i, 2] <- aj[2, 3, dim(aj)[3]] * los[i, 3]
        aj <- array(c(diag(1, 4, 4), aj), c(4, 4, (dim(aj)[3] + 1)))
    }
    
##     phi2[, 3] <- out$phi2case; phi2[, 2] <- out$phi2control
##     phi3[, 3] <- out$phi3case; phi3[, 2] <- out$phi3control
    indi <- apply(x$n.event, 3, function(x) {sum(x[1, ]) != 0})
    wait.times <- x$time[indi]
    wait.prob <- x$est["0", "0", ][indi]
    my.weights <- diff(c(0, 1 - wait.prob))
    
    pp <- x$n.risk[-1, ]
    ev.last <- apply(x$n.event[, , dims[3]], 1, sum)[1:2]
    pp <- rbind(pp, pp[nrow(pp), ] - ev.last)
    filtre <- pp[, 1] <= 0 | pp[, 2] <= 0

    tmp <- list(los)#, phi2, phi3)
    estimates <- lapply(tmp, function(z) {
        ldiff <- z[, 3] - z[, 2]
        ldiff[filtre] <- 0
        estimate <- matrix(ldiff[is.element(z[, 1], wait.times)], nrow = 1) %*%
            matrix(my.weights, ncol=1)
        estimate
    })
    
    e.phi.w1 <- e.phi.w23 <- my.weights1 <- my.weights23 <- NULL
    if (aw) {
        cif1 <- cumsum(c(1, x$est["0", "0", 1:(dims[3] - 1)]) * tr.mat[1, 2, ])
        my.weights1 <- diff(c(0, cif1[indi])) / cif1[length(cif1)]
        cif23 <- cumsum(c(1, x$est["0", "0", 1:(dims[3] - 1)]) *
                        (tr.mat[1, 3, ] + tr.mat[1, 4, ]))
        my.weights23 <- diff(c(0, cif23[indi])) / cif23[length(cif23)]
        weights.aw <- list(my.weights1, my.weights23)
        estimates.aw <- lapply(weights.aw, function(z) {
            ldiff <- los[, 3] - los[, 2]
            ldiff[filtre] <- 0
            estimate <- matrix(ldiff[is.element(los[, 1], wait.times)], nrow = 1) %*%
                matrix(z, ncol = 1)
            estimate
        })
        e.phi.w1 <- estimates.aw[[1]]
        e.phi.w23 <- estimates.aw[[2]]
    }
    res <- list(e.phi = estimates[[1]], phi.case = los[, 3],
                phi.control = los[, 2], ## e.phi2 = estimates[[2]],
##                phi2.case = phi2[, 3]
##                phi2.control = phi2[, 2]
##                 e.phi3 = estimates[[3]], phi3.case = phi3[, 3],
##                 phi3.control = phi3[, 2],
                weights = my.weights,
                w.time = wait.times, time = x$time, e.phi.weights.1 = e.phi.w1,
                e.phi.weights.other = e.phi.w23, weights.1 = my.weights1,
                weights.other = my.weights23)
    res
}

    

clos.nocp_pure <- function(x, tr.mat, aw) {
    dims <- dim(x$est)
    times <- if (sum(x$n.event[, , dims[3]]) != 0) x$time else x$time[-length(x$time)]
    los <- matrix(rep(times, 3), ncol = 3, byrow = FALSE)
    lt <- length(times)
    tau <- max(times)
    los[length(times)-1, 2:3] <- rep(tau, 2)
    aj <- array(NA, c(3, 3, 1))
    aj[, , 1] <- diag(1, 3, 3)
    
    funn <- function(x, y) {
        x %*% aj[ , , y]
    }
    
    for (i in (lt - 2):1) {
        print(i)
        diffs <- diff(times[(i + 1):length(times)])
        mat <- tr.mat[ , , length(x$time[x$time <= times[i + 1]])]
        aj <- array(apply(X = diag(1:dim(aj)[3]), 1, funn,
                          x = mat), c(3, 3, dim(aj)[3]))
        los[, 3][i] <- times[i + 1] + matrix(diffs, nrow=1) %*%
            matrix(aj[2, 2, ], ncol = 1)
        los[, 2][i] <- times[i + 1] + matrix(diffs, nrow=1) %*%
            matrix((aj[1, 1, ] + aj[1,2,]),ncol=1)
        aj <- array(c(diag(1, 3, 3), aj), c(3, 3, (dim(aj)[3] + 1)))
    }

    pp <- x$n.risk[-1, ]
    ev.last <- apply(x$n.event[, , dims[3]], 1, sum)[1:2]
    pp <- rbind(pp, pp[nrow(pp), ] - ev.last)
    filtre <- pp[, 1] <= 0 | pp[, 2] <= 0
    
    indi <- apply(x$n.event, 3, function(x) {sum(x[1, ]) != 0})
    wait.times <- x$time[indi]
    wait.prob <- x$est["0", "0", ][indi]
    los.diff <- los[, 3] - los[, 2]
    los.diff[filtre] <- 0
    my.weights <- diff(c(0, 1 - wait.prob))
    estimate <- matrix(los.diff[is.element(los[, 1], wait.times)], nrow = 1) %*%
        matrix(my.weights, ncol=1)
    e.phi.w1 <- e.phi.w2 <- my.weights1 <- my.weights2 <- NULL
    if (aw) {
        cif1 <- cumsum(c(1, x$est["0", "0", 1:(dims[3] - 1)]) * tr.mat[1, 2, ])
        my.weights1 <- diff(c(0, cif1[indi])) / cif1[length(cif1)]
        cif2 <- cumsum(c(1, x$est["0", "0", 1:(dims[3] - 1)]) * tr.mat[1, 3, ])
        my.weights2 <- diff(c(0, cif2[indi])) / cif2[length(cif2)]
        weights.aw <- list(my.weights1, my.weights2)
        estimates.aw <- lapply(weights.aw, function(z) {
            ldiff <- los[, 3] - los[, 2]
            ldiff[filtre] <- 0
            estimate <- matrix(ldiff[is.element(los[, 1], wait.times)], nrow = 1) %*%
                matrix(z, ncol = 1)
            estimate
        })
        e.phi.w1 <- estimates.aw[[1]]
        e.phi.w2 <- estimates.aw[[2]]
    }
    res <- list(e.phi = estimate[[1]], phi.case = los[, 3],
                phi.control = los[, 2], weights = my.weights,
                w.time = wait.times, time = x$time, e.phi.weights.1 = e.phi.w1,
                e.phi.weights.other = e.phi.w2, weights.1 = my.weights1,
                weights.other = my.weights2, aj = aj)
    res
}


boot.clos <- function(data, state.names, tra, cens.name, s = 0, nboot) {
     res <- double(nboot)
     for (i in seq_len(nboot)) {
          index <- sample(unique(data$id), replace = TRUE)
          inds <- new.id <- NULL
          for (j in seq_along(index)){
               ind <- which(data$id == index[j])
               new.id <- c(new.id, rep(j, length(ind)))
               inds <- c(inds, ind)
          }
          dboot <- cbind(data[inds, ], new.id)
          dboot[, which(names(dboot) == "id")]
          dboot$id <- dboot$new.id
          tr.prob <- etm(dboot, state.names, tra, cens.name, s, cova = FALSE)
          res[i] <- clos(tr.prob)$e.phi
     }
     res
}

