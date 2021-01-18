## FIXME: citest_array: Clean up code; use modern tabXXX functions; add parametric bootstrap

###################################################################################
#'
#' @title Test for conditional independence in a contingency table
#' @description Test for conditional independence in a contingency table
#'     represented as an array.
#' @name citest-array
#' 
###################################################################################
#'
#' @param x An array of counts with named dimnames.
#' @param set A specification of the test to be made. The tests are of the form u
#'     and v are independent condionally on S where u and v are variables and S
#'     is a set of variables. See 'details' for details about specification of
#'     \code{set}.
#' @param statistic Possible choices of the test statistic are \code{"dev"} for
#'     deviance and \code{"X2"} for Pearsons X2 statistic.
#' @param method Method of evaluating the test statistic. Possible choices are
#'     \code{"chisq"}, \code{"mc"} (for Monte Carlo) and \code{"smc"} for
#'     sequential Monte Carlo.
#' 
#' @param adjust.df Logical. Should degrees of freedom be adjusted for
#'     sparsity?
#' @param slice.info Logical. Should slice info be stored in the
#'     output?
#' @param L Number of extreme cases as stop criterion if method is
#'     \code{"smc"} (sequential Monte Carlo test); ignored otherwise.
#' @param B Number (maximum) of simulations to make if method is
#'     \code{"mc"} or \code{"smc"} (Monte Carlo test or sequential
#'     Monte Carlo test); ignored otherwise.
#' @param ...  Additional arguments.

#' @details \code{set} can be 1) a vector or 2) a right-hand sided formula in
#'     which variables are separated by '+'. In either case, it is tested if the
#'     first two variables in the \code{set} are conditionally independent given
#'     the remaining variables in \code{set}.  (Notice an abuse of the '+'
#'     operator in the right-hand sided formula: The order of the variables does
#'     matter.)
#' 
#' If \code{set} is \code{NULL} then it is tested whether the first two
#' variables are conditionally independent given the remaining variables.
#' 

#' @return An object of class `citest` (which is a list).
#' @author Søren Højsgaard, \email{sorenh@@math.aau.dk}
#' @seealso \code{\link{ciTest}}, \code{\link{ciTest_df}},
#'     \code{\link{ciTest_mvn}}, \code{\link{chisq.test}}
#' @keywords htest
#' @examples
#'
#' data(lizard)
#'
#' ## lizard is has named dimnames
#' names( dimnames( lizard ))
#' ## checked with
#' is.named.array( lizard )
#' 
#' ## Testing for conditional independence:
#' # the following are all equivalent:
#' ciTest(lizard, set=~diam + height + species)
#' # ciTest(lizard, set=c("diam", "height", "species"))
#' # ciTest(lizard, set=1:3)
#' # ciTest(lizard)
#' # (The latter because the names in lizard are as given above.)
#'
#' ## Testing for marginal independence
#' ciTest(lizard, set=~diam + height)
#' ciTest(lizard, set=1:2)
#'
#' ## Getting slice information:
#' ciTest(lizard, set=c("diam", "height", "species"), slice.info=TRUE)$slice
#' 
#' ## Do Monte Carlo test instead of usual likelihood ratio test. Different
#' # options:
#'
#' # 1) Do B*10 simulations divided equally over each slice: 
#' ciTest(lizard, set=c("diam", "height", "species"), method="mc", B=400)
#' # 2) Do at most B*10 simulations divided equally over each slice, but stop
#' # when at most L extreme values are found
#' ciTest(lizard, set=c("diam", "height", "species"), method="smc", B=400)

#' @rdname citest-array
ciTest_table <- function(x, set=NULL, statistic="dev", method="chisq", adjust.df=TRUE, slice.info=TRUE, L=20, B=200, ...){

    statistic <- match.arg(toupper(statistic), c("DEV",   "X2"))
    method    <- match.arg(toupper(method),    c("CHISQ", "MC", "SMC"))
    
    if (is.null(set)){
        set <- names( dimnames(x) )
    } else {
        if ( inherits(set, "integer") || inherits(set, "numeric") ){
            x   <- tabMarg(x, set)
        } else
            if (inherits(set,c("formula","character"))){
                set <- unlist(rhsFormula2list(set))
                vn  <- names(dimnames(x))
                set <- vn[pmatch(set, vn)]
                x   <- tabMarg(x, set)
            }
    }

    switch(method,
           "CHISQ"={
               .CI_X2_prim(x, statistic=statistic, adjust.df=adjust.df,
                           slice.info=slice.info)
           },
           "MC"=,"SMC"={
               .CI_SMC_prim(x, statistic=statistic, method=method,
                            slice.info=slice.info, L=L, B=B)
           }
           )
}


###
### CIP test; asymptotic, based on either deviance or Pearsons X2
###

# 'x' is a named array; let u be the first name, w be the second and R
# denote the 'rest'. The function tests u _|_ w | R.

.CI_X2_prim <- function(x, statistic="DEV", adjust.df=TRUE, slice.info=TRUE){
    
    statistic <- match.arg(toupper(statistic), c("DEV",   "X2"))

    vn    <- names(dimnames(x))
    di    <- dim(x)
    u     <- vn[1]
    w     <- vn[2]
    R     <- vn[-(1:2)]
    dim.u <- di[1]
    dim.w <- di[2]
    dim.R  <- prod(di[-(1:2)])
    
    t.uR <- tabMarg(x, c(u, R))
    t.wR <- tabMarg(x, c(w, R))

    ##str(list(t.uR=t.uR, t.wR=t.wR, R=R, vn=vn))
    fit.table <- fit2way_(t.uR, t.wR, R, vn)
    
    ## Evaluate test statistic
    ## FIXME There are functions for that in other functions
    if (statistic == "DEV"){           ## Deviance
        tobs  <- 2 * x * log(x / fit.table)
    } else {                           ## Pearson X2
        tobs <- (x - fit.table)^2 / fit.table   
    }
    tobs[!is.finite(tobs)] <- 0

    tobsGlobal <- sum(tobs)

    ## Calculate df with or without adjustment for sparsity
    if (!adjust.df) dofSlice <- rep.int((dim.u - 1) * (dim.w - 1), dim.R)
    else {
        t.uRmat    <- matrix(t.uR, nrow=dim.R, byrow=TRUE)
        t.wRmat    <- matrix(t.wR, nrow=dim.R, byrow=TRUE)
        
        z <- (t.uRmat > 0) * 1
        dim.u.adj <- if (!is.null(dim(z))) rowSums(z) else sum(z)
        
        z <- (t.wRmat > 0) * 1
        dim.w.adj <- if (!is.null(dim(z))) rowSums(z) else sum(z)
        
        d1         <- dim.u.adj - 1
        d1[d1 < 0] <- 0
        d2         <- dim.w.adj - 1
        d2[d2 < 0] <- 0
        dofSlice   <- d1 * d2
    } 
    
    dofGlobal <- sum(dofSlice)
    pGlobal   <- 1 - pchisq(tobsGlobal, dofGlobal)

    if (length(R) && slice.info){
        tobsSlice <- rowSums(matrix(tobs, nrow=dim.R, byrow=TRUE))
        pSlice    <- 1 - pchisq(tobsSlice, df=dofSlice)
        sliceInfo <- list(statistic=tobsSlice, p.value=pSlice, df=dofSlice)
        des       <- expand.grid(dimnames(x)[-(1:2)])
        slice     <- cbind(as.data.frame(sliceInfo[1:3]), des)
    } else slice <- NULL
    
    ans <- list(statistic=tobsGlobal, p.value=pGlobal, df=dofGlobal, statname=statistic,
                method="CHISQ", adjust.df=adjust.df, varNames=vn, slice=slice)

    class(ans) <- "citest"
    ans
}


###
### CIP test; exact, based on sequential monte carlo
###

.CI_SMC_prim <- function(x, statistic="DEV", method="SMC", L=50, B=200, slice.info=FALSE){
    
    statistic <- match.arg(toupper(statistic), c("DEV",   "X2"))
    switch(method,
           "MC"={
               zzz <- .CI_MC_prim(x, statistic=statistic, B=10*B, slice.info=slice.info)
           },
           "SMC"={
               zzz <- .CI_MC_prim(x, statistic=statistic, B=B,    slice.info=slice.info)
               tot <- as.numeric(zzz[c("n.extreme","B")])
               if (slice.info){
                   slice <- zzz$slice
               }
               repeat{
                   if (tot[1] > L) break
                   zzz<- .CI_MC_prim(x, statistic=statistic, B=B, slice.info=slice.info)
                   tot <- tot + as.numeric(zzz[c("n.extreme","B")])
                   if (slice.info){
                       slice[,"n.extreme"] <- slice[,"n.extreme"] +
                           zzz$slice[,"n.extreme"]
                   }
               }
               
               zzz[c("p.value", "n.extreme", "B")] <- c(tot[1] / tot[2], tot)
               zzz["method"] <- "SMC"
               
               if (slice.info){
                   slice[,"p.value"] <- slice[,"n.extreme"] / tot[2]
                   zzz[["slice"]]    <- slice
               }
           })
    zzz
}


###
### CIP test; exact, based on monte carlo
###
.CI_MC_prim <- function(x, statistic="DEV", B=100, slice.info=FALSE){

    statistic <- match.arg(toupper(statistic), c("DEV",   "X2"))

    # Calculates deviance for independence model in r x c table
    .devFun2 <- function(obs, fit){ 
        ii  <- obs * fit > 0
        2 * sum(obs[ii] * log(obs[ii] / fit[ii]))
    }

    # Calculates deviance for independence model in r x c table
    .X2Fun2 <- function(obs, fit){ 
        ii  <- obs * fit > 0
        a   <- (obs - fit)^2 / fit
        sum(a[ii])
    }

    .statFun <- if (statistic=="DEV") .devFun2 else .X2Fun2
    
    dn     <-  dim(x)
    v.idx  <-  seq_len(length(dn))
    v1R    <-  v.idx[-2]
    v2R    <-  v.idx[-1]
    dim12  <-  dim(x)[1:2]
    dim.R  <-  prod(dn[-(1:2)]) ## Careful when R is empty
    
    ## Marginal tables for (v1,R) and (v2,R) as matrices. Each row
    ## is a configuration of R
    t1R  <-  tabMarg(x, v1R)
    t2R  <-  tabMarg(x, v2R)
    t1R  <-  matrix(t1R, nrow=dim.R, byrow=TRUE)
    t2R  <-  matrix(t2R, nrow=dim.R, byrow=TRUE)
    xmat <-  matrix(x,   nrow=dim.R, byrow=TRUE)
    
    ## Find observed statistics for each slice
    tobs.slice <- vector("numeric", dim.R)
    for (ii in seq_len(dim.R)){
        r.sum    <- t1R[ii, ]
        c.sum    <- t2R[ii, ]
        expected <- outerPrim(r.sum, c.sum) / sum(r.sum)
        mm       <- xmat[ii, ]
        dim(mm)  <- dim12
        tobs.slice[ii] <- .statFun(mm, expected)
    }
    
    ## Find reference distribution for each slice
    tref.slice      <- matrix(NA, nrow=dim.R, ncol=B)
    n.extreme.slice <- vector("numeric", dim.R)
    for (ii in seq_len(nrow(t1R))){
        r.sum    <- t1R[ii,]
        c.sum    <- t2R[ii,]
        expected <- outerPrim(r.sum, c.sum) / sum(r.sum)
        zzz      <- r2dtable(B, r.sum, c.sum)
        for (kk in seq_len(B))
            tref.slice[ii, kk] <- .statFun(zzz[[kk]],expected)
        n.extreme.slice[ii] <- sum(tobs.slice[ii] < tref.slice[ii, ])
    }
    
    tref.total  <- colSums(tref.slice)
    tobs.total  <- sum(tobs.slice)
    n.extreme   <- sum(tobs.total < tref.total)
    p.value.slice <- n.extreme.slice / B
    p.value.total <- n.extreme / B
    
    if (slice.info){
        des   <- expand.grid(dimnames(x)[-(1:2)])
        slice <- cbind(data.frame(statistic=tobs.slice, n.extreme=n.extreme.slice,
                                  p.value=p.value.slice, df=NA), des)
    } else {
        slice=NULL
    }
    
    ans <- list(statistic=tobs.total, p.value=p.value.total, df=NA, statname=statistic,
                method="MC", varNames=names(dimnames(x)), n.extreme=n.extreme, B=B, slice=slice)
    class(ans) <- "citest"
    ans
}

