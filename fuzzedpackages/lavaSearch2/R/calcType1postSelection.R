### calcType1postSelection.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 31 2017 (16:42) 
## Version: 
## last-updated: feb  4 2019 (09:40) 
##           By: Brice Ozenne
##     Update #: 130
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Documentation - calcType1postSelection 
##' @title Compute the Type 1 Error After Selection [EXPERIMENTAL]
##' @description Compute the type 1 error after selection [EXPERIMENTAL].
##' @name calcType1postSelection
##' 
##' @param level [numeric 0-1] expected coverage.
##' @param mu [numeric vector] the expectation of the joint distribution of the test statistics
##' @param Sigma [matrix] the variance-covariance of the joint distribution of the test statistics.
##' @param quantile.previous [numeric] significance quantile used at the previous step.
##' @param df [integer > 0] the degree of freedom of the joint Student's t distribution.
##' Only used when \code{distribution="pvmt"}.
##' @param n [integer > 0] number of points for the numerical integration
##' @param distribution [character] distribution of the test statistics.
##' Can be \code{"pmvnorm"} (normal distribution) or \code{"pvmt"} (Student's t distribution)
##' @param correct [logical] if true, correct the level to account for previous testings.
##' @param ... arguments passed to \code{\link{intDensTri}}.  
##'
##' @details The number of tests at the current step (i.e. after selection) is assumed to be
##' one less than the number of tests at the previous step (i.e. before selection).
##'
##' Arguments \code{mu} and \code{Sigma} must contain the moments for the vector of test statistics
##' before and after selection (in that order).
##' 
##' @return [numeric] the type 1 error.
##' @author Brice Ozenne
##' @examples
##' library(mvtnorm)
##' n <- 350
##' 
##' #### only 2 tests
##' Sigma <- rbind(c(1,0,0),c(0,1,1),c(0,1,1))
##' z2 <- qmvnorm(0.95, mean = rep(0,2), sigma = Sigma[1:2,1:2], tail = "both.tails")$quantile
##' 
##' ## no selection since strong effect
##' mu <- c(10,0,0)
##' calcType1postSelection(0.95, quantile.previous = z2, distribution = "gaussian",
##'                         mu = mu, Sigma = Sigma, correct = TRUE)
##'
##' ## strong selection
##' \dontrun{
##' mu <- c(0,0,0)
##' levelC <- calcType1postSelection(0.95, quantile.previous = z2, distribution = "gaussian",
##'                         mu = mu, Sigma = Sigma)
##' print(levelC) # more liberal than without selection
##' calcType1postSelection(levelC, quantile.previous = z2, distribution = "gaussian",
##'                         mu = mu, Sigma = Sigma, correct = FALSE)
##' }
##'
##' #### 3 tests
##' Sigma <- diag(1,5,5)
##' Sigma[4,2] <- 1
##' Sigma[2,4] <- 1
##' Sigma[5,3] <- 1
##' Sigma[3,5] <- 1
##' 
##' z2 <- qmvnorm(0.95, mean = mu[1:3], sigma = Sigma[1:3,1:3], tails = "both.tails")$quantile
##' 
##' ## no selection since strong effect
##' \dontrun{
##' mu <- c(10,0,0,0,0)
##' calcType1postSelection(0.95, quantile.previous = z2, distribution = "gaussian",
##'                         mu = mu, Sigma = Sigma, correct = TRUE)
##' 
##' ## strong selection
##' mu <- c(0,0,0,0,0)
##' levelC <- calcType1postSelection(0.95, quantile.previous = z2,
##'                         mu = mu, Sigma = Sigma, distribution = "gaussian")
##' calcType1postSelection(levelC, quantile.previous = z2, distribution = "gaussian",
##'                         mu = mu, Sigma = Sigma, correct = FALSE)
##' }
##'
##' @concept post-selection inference


## * calcType1postSelection
##' @rdname calcType1postSelection
##' @export
calcType1postSelection <- function(level, mu, Sigma, quantile.previous, distribution, df,
                                    n = 10, correct = TRUE, ...){

    ## ** normalisation of the arguments
    p <- length(mu)
    if(p %% 2 == 0){
        stop("\'mu\' must have uneven length\n",
             "since it contains the means for the test statistics at the previous step (p+1)/2",
             "and those at the current step (p-1)/2")
    }
    if(NCOL(Sigma)!=p || NROW(Sigma)!=p){
        stop("\'Sigma\' must be a pxp matrix where p is the length of \'mu\' \n")
    }
    nTest <- (p-1)/2
    distribution <- match.arg(distribution, choices = c("student","gaussian"))
    
    if(distribution=="student"){
        pdist <- "pmvt"
        qdist <- "qmvt"
        args.dist <- list(delta = rep(0,nTest),
                          sigma = Sigma[1:nTest,1:nTest,drop=FALSE],
                          tail = "both.tails", df = df)
    }else if(distribution=="gaussian"){
        pdist <- "pmvnorm"
        qdist <- "qmvnorm"
        args.dist <- list(delta = rep(0,nTest),
                          sigma = Sigma[1:nTest,1:nTest,drop=FALSE],
                          tail = "both.tails")
    }    

    ## ** wrapper
    warper <- function(x){ # x <- 0.95## level
        # P[A|B] = 1-P[\bar{A}|B]  = 1-P[\bar{A},B]/P[B]
        ## Current quantile
        args.dist$p <- x
        newQuantile <- do.call(qdist, args = args.dist)$quantile
        
        ## Compute type one error
        num <- intDensTri(mu = mu, Sigma = Sigma, df = df, n=n,
                          x.min = quantile.previous, z.max = rep(newQuantile, nTest),
                          distribution = pdist, ...) # P[\bar{A},B]
        # autoplot(num, coord.plot = c("x","z"))
        denum <- intDensTri(mu = mu[1:(nTest+1)], Sigma = Sigma[1:(nTest+1),1:(nTest+1)], df = df, n=n,
                            x.min = quantile.previous, z.max = NULL,
                            distribution = pdist, ...) # P[B]
        out <- 1-num$value/denum$value
        return(out)
    }

    ## ** compute type1 error
    if(correct){
        out <- stats::optim(par = level, fn = function(x){
            obj <- abs((1-level)-warper(x))
            # cat("level:",x," obj:",obj,"\n")
            return(obj)
        },
        method = "Brent", lower = 0.51, upper = 0.999,
        control = list(abstol = 1e-4, reltol = 1e-4)
        )$par

    }else{
        out <- warper(level)
    }
    
    return(out)
}

#----------------------------------------------------------------------
### calcType1postSelection.R ends here
