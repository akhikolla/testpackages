### InTriGauss.R --- 
##----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 14 2017 (11:49) 
## Version: 
## last-updated: maj 23 2018 (09:42) 
##           By: Brice Ozenne
##     Update #: 495
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:



## * Documentation - intDensTri
##' @title Integrate a Gaussian/Student Density over a Triangle
##' @name intDensTri
##' @description Consider a univariate random variable X,
##' two multivariate random variables Y and Z,
##' and t1 and t2 two real numbers.
##' This function can compute either
##' P[|X|>t1,|X]>|Y1|,...,|X]>|Yp|] if zmin is not specified,
##' P[|Z1|<t2,...,|Zq|<t2,|X|>t1,|X]>|Y1|,...,|X]>|Yp|] if zmin is specified.
##' 
##' @param mu [numeric vector] the expectation of the joint distribution.
##' @param Sigma [matrix] the variance-covariance of the joint distribution.
##' @param df [integer > 0] the degree of freedom of the joint Student's t distribution.
##' Only used when \code{distribution="pvmt"}.
##' @param n [integer > 0] number of points for the numerical integration.
##' @param x.min [numeric] the minimum value along the x axis.
##' @param z.max [numeric vector, optional] the maximum value along the z axis.
##' Define the dimension of Z.
##' @param type [character] the type of mesh to be used.
##' Can be \code{\"raw\"}, \code{\"double\"}, or \code{\"fine\"}.
##' @param prune [integer >0] number of standard deviations after which the domain ends along the x axis. 
##' @param proba.min [numeric 0-1] the probability used to find the maximum value along the x axis.
##' Only used if \code{prune} is not specified.
##' @param distribution [character] type of joint distribution.
##' Can be \code{"pmvnorm"} (normal distribution) or \code{"pvmt"} (Student's t distribution)
##' @details
##' Argument \code{type}: \itemize{
##' \item \code{\"raw\"}: mesh with points inside the domain
##' \item \code{\"double\"}: mesh with points outside the domain
##' \item \code{\"fine\"}: mesh with points inside the domain plus additional rectangles trying to fill the missing domain.
##' }
##'
##' Argument \code{Sigma} and \code{mu}:
##' define the mean and variance-covariance of the random variables X, Y, Z
##' (in this order). The length of the argument \code{z.max} is used to define the dimension of Z.
##' The dimension of X is always 1.
##'
##' @return A numeric.
##' 
##' @examples
##' library(mvtnorm)
##' 
##' p <- 2
##' Sigma <- diag(p)
##' mu <- rep(0, p)
##' 
##' ## bivariate normal distribution
##' z2 <- qmvt(0.975, mean = mu, sigma = Sigma, df = 1e3)$quantile
##'
##' # compute integral
##' intDensTri(mu = mu, Sigma = Sigma, n=5, x.min=0, type = "fine")$value-1/2
##' intDensTri(mu = mu, Sigma = Sigma, n=30, x.min=0, type = "raw")$value-1/2
##' intDensTri(mu = mu, Sigma = Sigma, n=50, x.min=0, type = "raw")$value-1/2
##'
##' intDensTri(mu = mu, Sigma = Sigma, df = 5, n=5, x.min=0, distribution = "pmvt")$value-1/2
##' res <- intDensTri(mu = mu, Sigma = Sigma, df = 5, n=10, x.min=0, distribution = "pmvt")
##' res$value-1/2
##' ggplot2::autoplot(res)
##' 
##' ## trivariate normal distribution
##' \dontrun{
##' p <- 3
##' Sigma <- diag(p)
##' mu <- rep(0, p)
##'
##' res2 <- intDensTri(mu = mu, Sigma = Sigma, n=5, x.min = 0, z.max = 10)
##' ggplot2::autoplot(res2)
##' ggplot2::autoplot(res2, coord.plot = c("x","z1"))
##' res2
##' }
##' 
##' #### when the distribution is far from 0
##' \dontrun{
##' eq1 <- intDensTri(mu = c(10,0), Sigma = diag(1,2), 
##'                   x.min = 2, n=10)
##' eq1$value-1
##' ggplot2::autoplot(eq1)
##'
##' eq2 <- intDensTri(mu = c(10,0,0), Sigma = diag(1,3),
##'                   x.min=2, z.max = 10, type = "raw",
##'                   n=10)
##' ggplot2::autoplot(eq2, coord.plot = c("y1","z1"))
##' eq2$value-1
##'
##' ## more variables
##' p <- 5
##' Sigma <- diag(p)
##' mu <- rep(0, p)
##'
##' res2 <- intDensTri(mu = mu, Sigma = Sigma, n=5, x.min = 1, z.max = c(2,2))
##' res2$grid
##' }
##'
##' @concept post-selection inference

## * intDensTri
##' @rdname intDensTri
##' @export 
intDensTri <- function(mu, Sigma, df, n, x.min, z.max = NULL,
                       type = "double", proba.min = 1e-6, prune = NULL, distribution = "pmvnorm"){

     interior <- NULL ## [:for CRAN check] subset
    
    ## ** normalize arguments
    type <- match.arg(type, c("raw","fine","double"))
    
    p <- length(mu)
    d.z <- length(z.max)
    d.y <- p - d.z - 1
    z.max <- unique(z.max)
    
    if(is.null(prune)){
        if(distribution=="pmvnorm"){
            prune <- ceiling(abs(qmvnorm(proba.min, sigma = diag(1,d.y))$quantile))
        }else if(distribution=="pmvt"){
            prune <- ceiling(abs(mvtnorm::qmvt(proba.min, sigma = diag(1,d.y), df = df)$quantile))
        }else {
            stop("distribution must be either \"pmvnorm\" or \"pmvt\" \n")
        }
    }
    coordX.max <- (mu[1]+x.min) + prune * sqrt(diag(Sigma)[1])
    
    ls.args <- list(mu,Sigma)    
    if(distribution=="pmvt"){
        names(ls.args) <- c("delta","sigma")
        ls.args$df <- df
    }else{
        names(ls.args) <- c("mean","sigma")
    }

    ## ** create the grid of points to integrate over
    resGrid <- createGrid(n,
                          xmin = x.min, xmax = coordX.max, d.y = d.y, 
                          d.z = d.z, zmax = z.max, 
                          fine = (type=="fine"), double = FALSE
                          )
    grid <- resGrid$grid
    seqNames.min <- resGrid$seqNames.min
    seqNames.max <- resGrid$seqNames.max

    if(type=="double"){
        grid.double <- createGrid(n,
                                  xmin = x.min, xmax = coordX.max, d.y = d.y, 
                                  d.z = d.z, zmax = z.max, 
                                  fine = FALSE, double = TRUE
                                  )$grid
        grid <- rbind(cbind(grid,interior=TRUE),
                      cbind(grid.double,interior=FALSE))
    }

    ## ** integration
    total.area <- 0
    grid$area <- apply(grid, 1, function(x){
        do.call(distribution, args = c(lower = list(c(x[seqNames.min])),
                                       upper = list(c(x[seqNames.max])),
                                       ls.args))
    })


    if(type=="double"){
        grid.interior <- subset(grid, interior==TRUE,select = c("area", "weight"))
        area.interior <- sum(grid.interior$area * grid.interior$weight)
        grid.exterior <- subset(grid, interior==FALSE,select = c("area", "weight"))
        area.exterior <- sum(grid.exterior$area * grid.exterior$weight)
        total.area <- area.interior + (area.exterior-area.interior)/2
    }else{        
        total.area <- sum(grid$area * grid$weight)
    }

    ## ** export
    out <- list()
    out$value <- total.area    
    out$prune <- prune
    out$type <- type
    out$grid <- grid
    class(out) <- "intDensTri"
    return(out)
}

## * autoplot.intDensTri
#' @title 2D-display of the Domain Used to Compute the Integral
#' @description 2D-display of the domain used to compute the integral.
#'
#' @param object output of the function \code{intDensTri}.
#' @param coord.plot [character vector] the x and y coordinates. Can be \code{"x"}, \code{"y1"} to \code{"yd"}, \code{"z"} if \code{zmin} was specified when calling \code{intDensTri}.
#' @param plot [logical] should the plot be displayed?
#' @param ... [internal] Only used by the generic method.
#'
#' @return A \code{ggplot} object.
#' @seealso \code{\link{intDensTri}}
#' 
#' @method autoplot intDensTri
#' @export
autoplot.intDensTri <- function(object, coord.plot=c("x","y1"), plot = TRUE, ...){

    if(length(coord.plot) != 2){
        stop("coord.plot must have length 2 \n")
    }

    gg.data <- object$grid
    gg.data$index <- as.factor(gg.data$index)
    gg.data$weight <- as.factor(gg.data$weight)
    
    if("x" %in% coord.plot == FALSE){
        x.ref <- min(abs(unique(gg.data$x.min)))
        gg.data <- rbind(gg.data[gg.data$x.min == x.ref,,drop=FALSE],
                         gg.data[gg.data$x.max == -x.ref,,drop=FALSE])
    }

    if(object$type=="fine"){
        gg.grid <- ggplot2::ggplot(gg.data, aes_string(xmin = paste0(coord.plot[1],".min"), ymin = paste0(coord.plot[2],".min"),
                                                      xmax = paste0(coord.plot[1],".max"), ymax = paste0(coord.plot[2],".max"),
                                                      fill = "index", alpha = "weight"))
        gg.grid <- gg.grid + ggplot2::scale_alpha_manual(values = c(0.45,1))
    }else if(object$type=="raw"){
        gg.grid <- ggplot2::ggplot(gg.data, aes_string(xmin = paste0(coord.plot[1],".min"), ymin = paste0(coord.plot[2],".min"),
                                                      xmax = paste0(coord.plot[1],".max"), ymax = paste0(coord.plot[2],".max"),
                                                      fill = "index"))        
    }else if(object$type=="double"){
        gg.grid <- ggplot2::ggplot(gg.data, aes_string(xmin = paste0(coord.plot[1],".min"), ymin = paste0(coord.plot[2],".min"),
                                                      xmax = paste0(coord.plot[1],".max"), ymax = paste0(coord.plot[2],".max"),
                                                      fill = "index", alpha = "interior"))
        gg.grid <- gg.grid + ggplot2::scale_alpha_manual(values = c(0.4,1))
    }
    gg.grid <- gg.grid + ggplot2::geom_rect() + ggplot2::xlab(coord.plot[1]) + ggplot2::ylab(coord.plot[2])
    gg.grid <- gg.grid + ggplot2::geom_abline(slope = 1,color = "black") + ggplot2::geom_abline(slope = -1,color = "black")

    if(plot){
        print(gg.grid)
    }
    return(invisible(gg.grid))
}


## * print.intDensTri
#' @method print intDensTri
#' @export
print.intDensTri <- function(x, ...){
    interior <- NULL
    
    x.min <- min(abs(x$grid$x.min))
    x.max <- max(abs(x$grid$x.max))
    n.rectangle <- switch(x$type,
                          "double" = NROW(x$grid[interior==TRUE]),
                          NROW(x$grid)
                          )
    
    cat("integral=",x$value,"\n",
        "computed with ",n.rectangle," rectangles \n",
        "with x ranging from ",x.min," to ",x.max,"\n",
        sep = "")

    return(NULL)
}
#----------------------------------------------------------------------
### InTriGauss.R ends here
