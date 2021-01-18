##########################################################################
#' @title Fitting Log-Linear Models by Message Passing
#' 
#' @description Fit log-linear models to multidimensional contingency
#'     tables by Iterative Proportional Fitting.
#'
#' @name loglin-effloglin
#' 
##########################################################################
#' 
#' @details The function differs from \code{loglin} in that 1) data
#'     can be given in the form of a list of sufficient marginals and
#'     2) the model is fitted only on the cliques of the triangulated
#'     interaction graph of the model. This means that the full table
#'     is not fitted, which means that \code{effloglin} is efficient
#'     (in terms of storage requirements). However \code{effloglin} is
#'     implemented entirely in R and is therefore slower than
#'     \code{loglin}. Argument names are chosen so as to match those
#'     of loglin()
#' @param table A contingency table
#' @param margin A generating class for a hierarchical log--linear model
#' @param fit If TRUE, the fitted values are returned.
#' @param eps Convergence limit; see 'details' below.
#' @param iter Maximum number of iterations allowed
#' @param print If TRUE, iteration details are printed.
#' @return A list.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{loglin}}
#'
#' @references Radim Jirousek and Stanislav Preucil (1995). On the
#'     effective implementation of the iterative proportional fitting
#'     procedure. Computational Statistics & Data Analysis Volume 19,
#'     Issue 2, February 1995, Pages 177-189
#' 
#' @keywords models
#' @examples
#' 
#' data(reinis)
#' glist <-list(c("smoke", "mental"), c("mental", "phys"),
#'              c("phys", "systol"), c("systol", "smoke"))
#' 
#' stab <- lapply(glist, function(gg) tabMarg(reinis, gg))
#' fv3 <- effloglin(stab, glist, print=FALSE)
#' 
#' @export effloglin
effloglin <- function(table, margin, fit=FALSE, eps=0.01, iter=20, print=TRUE){

    amat <- ugList(margin, result="matrix")
    vn   <- colnames(amat)
    tri  <- triangulateMAT(amat)
    rip  <- ripMAT(tri)

    cliq     <- rip$cliques
    len.cliq <- length(cliq)

    ## get "host clique" for each generator
    ## FIXME: (effloglin) use general "get.host.clique" function.
    ghost   <- rep(NA, length(margin))
    seqcliq <- seq_along(cliq)
    for (kk in 1:length(margin)){
        ##cat("kk:", kk,"\n")
        gg <- margin[[kk]]
        for (ii in seqcliq){
            ##cat ("ii", ii, "\n")
            zz <- match(gg, cliq[[ii]])
            if (!any(is.na(zz))){
                ghost[kk] <- ii
                break
            }
        }
    }

    if (is.array(table)){
        Nobs   <- sum(table)
        ##stlist <- lapply(margin, function(xx) {tabMarg(table, xx)})
        stlist <- lapply(margin, function(xx) tabMarg(table, xx))
    } else {
        Nobs   <- sum(table[[1]])
        stlist <- table
    }

    zzz       <- unlist(lapply(stlist, dimnames), recursive=FALSE)
    vl        <- zzz[unique.default(names(zzz))]
    pot.list  <- lapply(cliq, function(cq)
                        parray(cq, levels=vl[cq], values=1, normalize="all"))


    ##   cat("effloglin\n")
    ##   print(as.data.frame.table(pot.list[[1]]))

    ## ## Potential list over cliques
    ## Clique marginals
    prob.list  <- propagateLS(pot.list, rip, initialize=TRUE)

    itcount  <- 1L
    logL     <- 0
    max.dif  <- vector("numeric", length(margin))
    repeat{
        ##cat(sprintf("---------- iteration: %i -----------\n", itcount))
        for (ss in seq_along(margin)){
            gg      <- margin[[ss]]
            st      <- stlist[[ss]]
            cq      <- cliq[[ghost[ss]]]
            cq.idx  <- ghost[ss]
            cpot    <- prob.list[[cq.idx]]
            ##adjust  <- tableOp(st, tabMarg(cpot, gg)*Nobs, "/")

            ##tm      <- tabMarg(cpot, gg)*Nobs
            tm      <- tabMarg(cpot, gg) * Nobs
            adjust  <- st / tm
            max.dif[ss] <- max(abs(log(adjust)))
            ##max.dif[ss] <- max(abs(st-tm))
            logL    <- logL + sum(st * log(adjust))
            ##pot.list[[cq.idx]] <- tableOp(pot.list[[cq.idx]], adjust, "*")
            ##pot.list[[cq.idx]] <- tableOp2(pot.list[[cq.idx]], adjust, `*`)
            pot.list[[cq.idx]] <- tabProd(pot.list[[cq.idx]], adjust)
            prob.list          <- propagateLS(pot.list, rip, initialize=TRUE)
        }

        if (print)
            cat("max deviation (obs-fitted):", max(max.dif), "\n")
        if (max(max.dif) < eps || itcount >= iter)
            break()
        itcount <- itcount + 1L
    }

    vl    <- unlist(lapply(stlist, dimnames), recursive=FALSE)[vn]
    nlev  <- unlistPrim(lapply(vl, length))
    gn    <- lapply(margin, match, vn)
    nparm <- .dim_loglin(gn, nlev)
    df    <- prod(nlev) - 1 - nparm
    
    ans <- list(potlist=pot.list, margin=margin, vn=vn, rip=rip, ghost=ghost,
                stlist=stlist, logL=logL, nparm=nparm, df=df)
    
### Create full joint:
    if (fit){
        pjoint <- prob.list[[1]]
        if (length(prob.list)>1){
            for (ii in 2:length(prob.list)){
                pjoint <- tableOp(pjoint, tableOp(prob.list[[ii]],
                                                  tabMarg(prob.list[[ii]], rip$sep[[ii]]),
                                                  ##tabMarg(prob.list[[ii]], rip$sep[[ii]]),
                                                  "/"),"*")
            }
        }
        ##pjoint <- tablePerm(pjoint, vn)*Nobs
        pjoint <- tabPerm(pjoint, vn) * Nobs
        ans <- c(ans, list(fit=pjoint))
    }
    ## class(ans) <- "effloglin"
    return(ans)
}





