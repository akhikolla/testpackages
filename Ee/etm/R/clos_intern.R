### To be used for single endpoint
clos.nocp <- function(x, aw, ratio) {

    dims <- dim(x$est)
    los <- matrix(rep(x$time, 3), ncol = 3, byrow = FALSE)
    tau <- max(x$data$exit)
    times <- x$time
    surv <- x$est[1, 1, ]

    ## Call to C++ function
    out <- .Call("los_nocp",
                 times,
                 x$delta.na,
                 tau)

    los[, 2] <- out$los0
    los[, 3] <- out$los1
    indi <- apply(x$n.event, 3, function(x) {sum(x[1, ]) != 0})
    wait.times <- x$time[indi]
    wait.prob <- x$est[1, 1, ][indi]

    pp <- x$n.risk[-1, ]
    ev.last <- apply(x$n.event[, , dims[3]], 1, sum)[1:2]
    pp <- rbind(pp, pp[nrow(pp), ] - ev.last)
    filtre <- pp[, 1] <= 0 | pp[, 2] <= 0

    if (ratio) {
        los.diff <- los[, 3] / los[, 2]
    } else {
        los.diff <- los[, 3] - los[, 2]
    }
    los.diff[filtre] <- 0
    my.weights <- diff(c(0, 1 - wait.prob))
    estimate <- matrix(los.diff[is.element(los[, 1], wait.times)], nrow = 1) %*%
        matrix(my.weights, ncol=1)
    
    e.phi.w1 <- e.phi.w2 <- my.weights1 <- my.weights2 <- NULL
    if (aw) {
        I <- diag(1, dims[1])
        tr.mat <- array(apply(x$delta.na, 3, "+", I), dim = dims)
        cif1 <- cumsum(c(1, x$est[1, 1, 1:(dims[3] - 1)]) * tr.mat[1, 2, ])
        my.weights1 <- diff(c(0, cif1[indi])) / cif1[length(cif1)]
        cif2 <- cumsum(c(1, x$est[1, 1, 1:(dims[3] - 1)]) * tr.mat[1, 3, ])
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
                weights.other = my.weights2)
    res
}


#######################################
## The competing risks version
#######################################
clos.cp <- function(x, aw, ratio) {

    dims <- dim(x$est)
    los <- matrix(rep(x$time, 3), ncol = 3, byrow = FALSE)
    phi2 <- matrix(data=c(x$time, rep(0, dims[3]), rep(0, dims[3])),
                   ncol=3, byrow=FALSE)
    phi3 <- matrix(data=c(x$time, rep(0, dims[3]), rep(0, dims[3])),
                   ncol=3, byrow=FALSE)
    ind.cens <- apply(x$n.event, 3, function(r) all(r == 0))
    times <- x$time
    tau <- max(x$time)
    
    out <- .Call("los_cp",
                 times,
                 x$delta.na,
                 tau)

    los[, 2] <- out$los0
    los[, 3] <- out$los1
    phi2[, 3] <- out$phi2case; phi2[, 2] <- out$phi2control
    phi3[, 3] <- out$phi3case; phi3[, 2] <- out$phi3control
    indi <- apply(x$n.event, 3, function(x) {sum(x[1, ]) != 0})
    wait.times <- x$time[indi]
    wait.prob <- x$est[1, 1, ][indi]
    my.weights <- diff(c(0, 1 - wait.prob))
    
    pp <- x$n.risk[-1, ]
    ev.last <- apply(x$n.event[, , dims[3]], 1, sum)[1:2]
    pp <- rbind(pp, pp[nrow(pp), ] - ev.last)
    filtre <- pp[, 1] <= 0 | pp[, 2] <= 0
    
    tmp <- list(los, phi2, phi3)
    estimates <- lapply(tmp, function(z) {
        if (ratio) {
            ldiff <- z[, 3] / z[, 2]
        } else {
            ldiff <- z[, 3] - z[, 2]
        }
        ldiff[filtre] <- 0
        estimate <- matrix(ldiff[is.element(z[, 1], wait.times)], nrow = 1) %*%
            matrix(my.weights, ncol=1)
        estimate
    })
    
    e.phi.w1 <- e.phi.w23 <- my.weights1 <- my.weights23 <- NULL
    if (aw) {
        I <- diag(1, dims[1])
        tr.mat <- array(apply(x$delta.na, 3, "+", I), dim = dims)
        cif1 <- cumsum(c(1, x$est[1, 1, 1:(dims[3] - 1)]) * tr.mat[1, 2, ])
        my.weights1 <- diff(c(0, cif1[indi])) / cif1[length(cif1)]
        cif23 <- cumsum(c(1, x$est[1, 1, 1:(dims[3] - 1)]) *
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
                phi.control = los[, 2], e.phi2 = estimates[[2]],
                phi2.case = phi2[, 3], phi2.control = phi2[, 2],
                e.phi3 = estimates[[3]], phi3.case = phi3[, 3],
                phi3.control = phi3[, 2], weights = my.weights,
                w.time = wait.times, time = x$time, e.phi.weights.1 = e.phi.w1,
                e.phi.weights.other = e.phi.w23, weights.1 = my.weights1,
                weights.other = my.weights23)
    res
}

    
