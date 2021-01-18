################################################################################
#' @title Mean, covariance and counts for grouped data (statistics for
#'     conditional Gaussian distribution).
#' 
#' @description \code{CGstats} provides what corresponds to calling
#'     \code{cow.wt} on different strata of data where the strata are defined by
#'     the combinations of factors in data.
#'
#' @name cg-stats
################################################################################
#'
#' @aliases CGstats CGstats.data.frame 
#' 
#' @param object A dataframe.
#' @param varnames Names of variables to be used.
#' @param homogeneous Logical; if TRUE a common covariance matrix is reported.
#' @param simplify Logical; if TRUE the result will be presented in a simpler
#'     form.
#' @return A list whose form depends on the type of input data and the varnames.
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{cov.wt}}
#' @keywords utilities
#' @examples
#' 
#' data(milkcomp)
#' # milkcomp <- subset(milkcomp, (treat %in% c("a", "b")) & (lactime %in% c("t1", "t2")))
#' # milkcomp <- milkcomp[,-1]
#' # milkcomp$treat 	<- factor(milkcomp$treat)
#' # milkcomp$lactime 	<- factor(milkcomp$lactime)
#' 
#' CGstats(milkcomp)
#' CGstats(milkcomp, c(1, 2))
#' CGstats(milkcomp, c("lactime", "treat"))
#' CGstats(milkcomp, c(3, 4))
#' CGstats(milkcomp, c("fat", "protein"))
#' 
#' CGstats(milkcomp, c(2, 3, 4), simplify=FALSE)
#' CGstats(milkcomp, c(2, 3, 4), homogeneous=FALSE)
#' CGstats(milkcomp, c(2, 3, 4), simplify=FALSE, homogeneous=FALSE)
#' 
#' 
#' 
#' @export CGstats
CGstats <- function(object, varnames=NULL, homogeneous=TRUE, simplify=TRUE){
  UseMethod("CGstats")
}

## FIXME: CGstats: Determination of cont.names, cont.idx is fragile
#' @export
CGstats.data.frame <- function(object, varnames=NULL, homogeneous=TRUE, simplify=TRUE){

    if (is.null(varnames)) varnames <- names(object)
    
    use.idx <- if (is.numeric(varnames)) varnames
               else match(varnames, names(object))
    
    zz <- unlist(lapply(object, is.numeric))
    cont.idx <- intersect(which(zz),  use.idx)
    disc.idx <- intersect(which(!zz), use.idx)
  
    cont.names <- names(zz)[cont.idx]
    disc.names <- names(zz)[disc.idx]
    
    CGstats_internal(object, disc.names, cont.names, homogeneous, simplify)
}

## October 2015: New implementation of CGstats_internal

## #' @rdname CGstats
## #'
## #' @param disc.names Vector of names of discrete variables
## #' @param cont.names Vector of names of continuous variables
CGstats_internal <- function(object, disc.names=NULL, cont.names=NULL, homogeneous=TRUE, simplify=TRUE){
    
    if (length(cont.names) == 0) {        
        xt    <- xtabs(~., data=object[, disc.names, drop=FALSE])
        ans   <- list(n.obs=xt)
        disc.levels <- dim(xt)
    } else {
        if (length(disc.names) == 0) {
            ans <- stats::cov.wt(object[, cont.names, drop=FALSE], method="ML")
            disc.levels <- NULL
        } else {
            ##cat("the mixed case\n")
            ans <- .moments.by( object, disc.names, cont.names )
            disc.levels <- dim(ans$n.obs)
            
            if (homogeneous)
                ans$cov <- .makeHomogeneous(ans$cov, ans$n.obs)
            
            if (simplify){
                ans$center <- t.default(do.call(rbind, ans$center))
                rownames(ans$center) <- cont.names
                if (!homogeneous){
                    ans$cov <- t.default(do.call(rbind, lapply(ans$cov, as.numeric)))
                }
            }
        }
    }
  c(ans, list(cont.names=cont.names, disc.names=disc.names, disc.levels=disc.levels))
  # res  <- c(ans, list(cont.names=cont.names, disc.names=disc.names, disc.levels=disc.levels))
  # ##class(res) <- "CGstats"
  # return(res)
}

## deprecated
print.CGstats <- function(x, ...){
    print.default(x[1:3])
    return(invisible(x))
}

## This is fine
.cov.wt <- function(x, method="ML"){
    if (!is.matrix(x)) x <- as.matrix(x)
    
    n.obs  <- nrow(x)
    center <- colSums(x) / n.obs
    N   <- if (method == "ML") n.obs else n.obs-1
    zz  <- x - rep(center, each=nrow(x))
    cov <- crossprod(zz) / N
    list(cov=cov, center=center, n.obs=n.obs)
}

## Should perhaps detect disc.names (as the factors in object) if not given
## Should check if object is dataframe.
.make.split.info <- function(object, disc.names){
    .dfcols2namevec <- function(x, sep='|'){
        ##x<-lapply(x, as.character)
        do.call(function(...)paste(..., sep=sep), x)
    }

    xt  <- xtabs(~., object[,disc.names,drop=FALSE])
    zz  <- as.data.frame.table(xt, useNA="always")[, seq_along(disc.names), drop=FALSE]
    uniqval <- .dfcols2namevec(zz)
    facstr  <- .dfcols2namevec(object[, disc.names, drop=FALSE])

    list(facstr=facstr, uniqval=uniqval, xt=xt)    
}


## Should check if object is dataframe.
.split.data.by <- function(object, disc.names){
    spl <- .make.split.info(object, disc.names)
    
    out <- vector("list", length(spl$uniqval))
    names(out)  <- spl$uniqval
    len.uniqval <- length(spl$uniqval)
    i=1
    while(i <= len.uniqval ){
        out[[i]] <- object[spl$uniqval[i] == spl$facstr, , drop=FALSE]
        i = i + 1
    }
    out
}

## Check of input
## Could have possibility of simply invoking split.data.by
.moments.by <- function(object, disc.names, cont.names){

    spl <- .make.split.info(object, disc.names)    

    len.uniqval <- length(spl$uniqval)
    SS.mean <- SS.cov <- vector("list", len.uniqval)
    names(SS.mean) <- names(SS.cov) <- spl$uniqval
    
    for (ii in seq_along(spl$uniqval)){
        xx <- object[spl$uniqval[ii] == spl$facstr, cont.names, drop=FALSE]
        zz <- .cov.wt(xx, method="ML" )
        SS.mean[[ii]] <- zz$center
        SS.cov[[ii]]  <- zz$cov
    }
    
    list(n.obs=spl$xt, center=SS.mean, cov=SS.cov)
}

## Works only for heterogeneous, mixed variables and simplify=TRUE
.extendCGstats <- function(CGstats) {
  n.i    <- as.numeric(CGstats[['n.obs']])
  N      <- sum(n.i)
  Q      <- nrow(CGstats[['center']])
    
  SSD    <- matrix(rowSums(.colmult(n.i, CGstats[['cov']])), nrow=Q)
  S      <- SSD / N
    
  mmm    <- CGstats[['center']]
  ttt    <- .colmult(n.i, mmm)
  quad   <- ttt %*% t(mmm)
  SS     <- SSD + quad
  
  ans <- c(CGstats, list(N=N, SSD=SSD, SS=SS))
  ##class(ans) <- "CGstats"
  ans
}

## Create homogeneous variance
## input: VV: list of covariance matrices; nn: vector of counts
.makeHomogeneous <- function(VV, nn){
    WW   <- VV[[1]]
    WW[] <- 0
    len.VV <- length(VV)
    i <- 1
    while(i <= len.VV){
        if (nn[i] > 0) WW <- WW + VV[[i]] * nn[i]
        i <- i + 1
    }
    WW <- WW / sum(nn)
    WW
}








