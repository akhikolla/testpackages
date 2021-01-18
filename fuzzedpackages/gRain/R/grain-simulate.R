####
#### Re-implementation of simulate() function - quite fast...
#### Bristol, March 2008
####

## ##################################################################
##
#' @title Simulate from an independence network
#' @description Simulate data from an independence network.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @name grain-simulate
##
## ##################################################################
#'
#' @param object An independence network.
#' @param nsim Number of cases to simulate.
#' @param seed An optional integer controlling the random number
#'     generatation.
#' @param \dots Not used.
#' @return A data frame

#' @references Søren Højsgaard (2012). Graphical Independence Networks
#'     with the gRain Package for R. Journal of Statistical Software,
#'     46(10), 1-26.  \url{http://www.jstatsoft.org/v46/i10/}.
#' @keywords models
#' @examples
#' 
#' tf <- system.file("huginex", "chest_clinic.net", package = "gRain")
#' 
#' chest <- loadHuginNet(tf, details=1)
#' simulate(chest,n=10) 
#' 
#' chest2 <- setFinding(chest, c("VisitToAsia", "Dyspnoea"),
#'                             c("yes", "yes"))
#' simulate(chest2, n=10)
#' @export simulate.grain

#' @export 
simulate.grain <- function(object, nsim=1, seed=NULL, ...){

    if (!isCompiled(object)){
        object <- compile(object, propagate=TRUE)
    } else {
        if (!isPropagated(object)){
            object <- propagate(object)
        }
    }
    
    plist  <- getgrain(object, "pot_equi")
    rp     <- rip(object)
    cqlist <- rp$cliques
    seplist <- rp$separators

    ## Init
    ans           <- matrix(0, nrow=nsim, ncol=length(nodeNames(object)))
    colnames(ans) <- nodeNames(object)
        
    ## Iterate
    for (ii in seq_along(cqlist)){
        ctab <- plist[[ii]]
        sep  <- seplist[[ii]]    ## What we condition on
        if (length(sep) == 0)
            res <- simulateArray(ctab, nsim=nsim)
        else {                            
            mtab <- tableMargin(ctab, sep)     ## FIXME: Old table-function
            ctab <- tableOp2(ctab, mtab, `/`)  ## FIXME: Old table-function

            vn   <- varNames(ctab)
            rest <- setdiff(vn, sep) ## Variables to be simulated
            ## cat(" ii: ", ii, " vn: ", vn, " rest:", rest, " sep:", sep, "\n")
           
            sepidx <- match(sep, vn) ## NOTE: Must do after table-division.
            res    <- matrix(0, nrow=nsim, ncol=length(rest))
            colnames(res) <- rest
           
            given    <- ans[, sep, drop=FALSE]
            ##cat("given:\n"); print(given)
            vals  <- unique(given)
            sc    <- cumprod(apply(vals, 2, max) )
            sc    <- c(1, sc)[1:length(sc)]
            key   <- ((given - 1) %*% sc) + 1
            ##cat(sprintf("key=%s\n", toString(key)));  browser()
            for(kk in unique(key)){
                n_sim   <- sum(kk == key)
                idx  <- given[match(kk, key), ]
                ## dd <<- list(x=ctab, nsim=n_sim, margin=sepidx, value.margin=idx)
                res[kk == key, ] <- simulateArray(ctab, nsim=n_sim,
                                                  margin=sepidx, value.margin=idx)
            }
        }         
        ans[, colnames(res)] <- res
    } ## for
    
    ns  <- nodeStates(object)
    vn  <- colnames(ans)
    out <- vector("list", ncol(ans))
    names(out) <- vn
    for (jj in 1:ncol(ans)){
        out[[jj]] <- factor(ans[,jj], levels=seq(ns[[jj]]))
        levels(out[[jj]]) <- ns[[jj]]
    }
    out <- as.data.frame(out)
    names(out) <- vn
    out
}
