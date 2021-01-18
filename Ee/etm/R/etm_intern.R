#####################################################################
### The etm machinery
### Arthur Allignol <arthur.allignol@uni-ulm.de
#####################################################################


.etm <- function(entry,
                 exit,
                 from,
                 to,
                 nstate,
                 s,
                 t,
                 covariance,
                 const_modif) {

    times <- sort(unique(exit[to != 0]))
    times <- times[times > s & times <= t]

    if (length(times) == 0) {
        zzz <- list()
        zzz$est <- array(diag(1, nrow = nstate, nstate), dim = list(nstate, nstate, 1))
        zzz$dna <- array(matrix(0, nrow = nstate, ncol = nstate), dim = list(nstate, nstate, 1))
        zzz$n.risk <- zzz$time <- zzz$n.event <- NULL
        return(zzz)
    } else {

        c_modif <- matrix(const_modif, length(times), nstate, byrow = TRUE)

        zzz <- .Call("gen_msm", times, entry, exit, from, to, nstate, c_modif)

    }

    if (covariance) {

        cov_etm <- .Call("cov_aj",
                         zzz$time,
                         zzz$est,
                         zzz$n.risk,
                         zzz$n.event,
                         zzz$dna)

        zzz$cov <- cov_etm

    }

    zzz
}

