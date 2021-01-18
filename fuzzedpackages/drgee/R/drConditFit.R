expit <- function(x){ (exp(-x) + 1)^(-1) }

getClusterInfo <- function( ya.dt ) {

    ## Clumsy workaround to get through the CRAN check
    id <- NULL
    y <- NULL
    a <- NULL
    y.sum <- NULL
    a.sum <- NULL
    size <- NULL
    y.disc <- NULL
    a.disc <- NULL
    ya.disc <- NULL

    ## Sum by cluster id
    ya.clust.dt <- ya.dt[, list( y.sum = sum(y), a.sum = sum(a), size = .N), by = id][, c("y.sum", "a.sum", "size", "id"), with=F]

    ## Identify outcome discordant clusters
    ya.clust.dt[, y.disc := ifelse( y.sum < size & y.sum > 0, 1, 0)]
    ## Identify exposure discordant clusters
    ya.clust.dt[, a.disc := ifelse( a.sum < size & a.sum > 0, 1, 0)]
    ## Identify doubly discordant clusters
    ya.clust.dt[, ya.disc := y.disc * a.disc]

    return ( ya.clust.dt )

}

pseudoPairsIdx <- function(ya.dt, twice = TRUE) {
    
    ## Find doubly discordant pairs among pseudopairs
    ## If each pseudopair should be used twice (both orders)
    if( twice ){
        
        ya.pairs.dt <- with( ya.dt,
                            merge(ya.dt,
                                  ya.dt,
                                  allow.cartesian = TRUE,
                                  suffixes = c(".1", ".2")
                                  )[(idx.1 != idx.2)]
                            )
        
    ## If each pseudopair should only be used once (in the order they have in original data)
    } else {
        
        ya.pairs.dt <- with( ya.dt,
                            merge(ya.dt,
                                  ya.dt,
                                  allow.cartesian = TRUE,
                                  suffixes = c(".1", ".2")
                                  )[(idx.1 < idx.2)]
                            )

    }

    return( ya.pairs.dt )
}

pseudoPairs <- function(object, ya.pairs.dt) {

    ## Find the ids for...
    
    ## ... the first member of each pair:
    id.1.pp <- object$id[ya.pairs.dt$idx.1]
    ## Find the outcomes, 
    y.1.pp <- object$y[ya.pairs.dt$idx.1,, drop = FALSE]
    colnames(y.1.pp) <- object$y.names
    ## find the exposure, 
    a.1.pp <- object$a[ya.pairs.dt$idx.1,, drop = FALSE]
    colnames(a.1.pp) <- object$a.names
    ## find the interaction terms, 
    x.1.pp <- object$x[ya.pairs.dt$idx.1,, drop = FALSE]
    ## find the covariates in the outcome model, 
    v.1.pp <- object$v[ya.pairs.dt$idx.1,, drop = FALSE]
    ## find the covariates in the exposure model
    z.1.pp <- object$z[ya.pairs.dt$idx.1,, drop = FALSE]
    
    ## ... the second member of each pair, 
    ## Find the outcomes, 
    y.2.pp <- object$y[ya.pairs.dt$idx.2,, drop = FALSE]
    colnames(y.2.pp) <- object$y.names
    ## find the exposure, 
    a.2.pp <- object$a[ya.pairs.dt$idx.2,, drop = FALSE]
    colnames(a.2.pp) <- object$a.names
    ## find the interaction terms, 
    x.2.pp <- object$x[ya.pairs.dt$idx.2,, drop=FALSE]
    ## find the covariates in the outcome model, 
    v.2.pp <- object$v[ya.pairs.dt$idx.2,, drop=FALSE]
    ## find the covariates in the exposure model, 
    z.2.pp <- object$z[ya.pairs.dt$idx.2,, drop=FALSE]
    
    ## Create outcomes with interactions
    yx.1.pp <- x.1.pp * as.vector(y.1.pp)
    yx.2.pp <- x.2.pp * as.vector(y.2.pp)
    ## Create exposures with interactions
    ax.1.pp <- x.1.pp * as.vector(a.1.pp)
    ax.2.pp <- x.2.pp * as.vector(a.2.pp)
    
    ## Create differences
    yx.diff.pp <- yx.1.pp - yx.2.pp
    ax.diff.pp <- ax.1.pp - ax.2.pp
    v.diff.pp <- v.1.pp - v.2.pp
    z.diff.pp <- z.1.pp - z.2.pp
    
    ## find the sum of interaction terms
    x.sum.pp <-  x.1.pp +  x.2.pp
    
    ## The exposure with interaction in the derived prospective model
    a.1.x.sum.pp <- x.sum.pp * as.vector( a.1.pp )
    ## The outcome with interaction in the derived retrospective model
    y.1.x.sum.pp <- x.sum.pp * as.vector( y.1.pp )
    
    ## find the original id for each pair
    id.pp <- object$id[ya.pairs.dt$idx.1]
    
    return( list(id.pp = id.pp, 
                 y.1.pp = y.1.pp,
                 a.1.pp = a.1.pp,
                 x.1.pp = x.1.pp, 
                 v.1.pp = v.1.pp, 
                 z.1.pp = z.1.pp, 
                 y.2.pp = y.2.pp,
                 a.2.pp = a.2.pp,
                 x.2.pp = x.2.pp, 
                 v.2.pp = v.2.pp, 
                 z.2.pp = z.2.pp,
                 yx.1.pp = yx.1.pp, 
                 yx.2.pp = yx.2.pp,                  
                 ax.1.pp = ax.1.pp, 
                 ax.2.pp = ax.2.pp,
                 v.diff.pp = v.diff.pp, 
                 z.diff.pp = z.diff.pp, 
                 ax.diff.pp = ax.diff.pp, 
                 yx.diff.pp = yx.diff.pp, 
                 x.sum.pp =  x.sum.pp, 
                 a.1.x.sum.pp =  a.1.x.sum.pp, 
                 y.1.x.sum.pp =  y.1.x.sum.pp ) )
}

drConditFit <- function(object, rootFinder = findRoots, intercept = TRUE) {

    ## Find the total number of observations
    n.obs <- length(object$id)
    
    ## Clumsy workaround to get through the CRAN check
    idx <- NULL
    id <- NULL
    y <- NULL
    a <- NULL
    y.1 <- NULL
    a.1 <- NULL
    y.2 <- NULL
    a.2 <- NULL
    
    ya.dt <- data.table(idx = seq_len(n.obs),
                        id = object$id,
                        y = as.integer(object$y),
                        a = as.integer(object$a) )
    names(ya.dt) <- c("idx", "id", "y", "a")    
    setkey(ya.dt, "id")

    ## Get info about clusters 
    clusters <- getClusterInfo(ya.dt)
    clusters.info <- c(colSums( clusters[, c("y.disc", "a.disc", "ya.disc"), with=F] ),
                       n.clust = nrow(clusters) )

    ## Select doubly discordant pseudo pairs
    ## Extract indices for pseudopairs    
    ya.pairs.dt <- pseudoPairsIdx(ya.dt)

    ## Select the pairs that are doubly discordant
    ya.dd.pairs.dt <- ya.pairs.dt[(y.1 != y.2) & (a.1 != a.2)]

    dd <- pseudoPairs(object, ya.dd.pairs.dt)
    pp.info <- c(n.pp = nrow(ya.pairs.dt),
                 n.dd.pp = nrow(ya.dd.pairs.dt) )
        
    ## Find the number of main parameters
    n.beta.params <- ncol(object$ax)
    ## Find the number of outcome nuisance parameters
    n.gamma.params <- ncol(object$v)
    ## Find the number of outcome nuisance parameters
    n.alpha.params <- ncol(object$z)

    ## Create a matrix of zeros
    beta.zeros <- matrix(rep(0, n.beta.params^2 ), nrow = n.beta.params)

    ## If we want to estimate nuisance model with intercepts
    if( intercept ) {

        coef.names <- c(object$ax.names,
                        object$ax.names,
                        c(object$x.names, object$v.names),
                        object$ax.names,
                        c(object$x.names, object$z.names))

        icpt <- rep(1, length(dd$y.1.pp))
        
        n.gamma.params <- n.gamma.params + 1
        n.alpha.params <- n.alpha.params + 1

        o.model.matrix <- cbind(dd$a.1.x.sum.pp, icpt, dd$v.diff.pp)
        o.fit <- geeFit(dd$y.1.pp, o.model.matrix, "logit")
        beta1.hat <- o.fit$coefficients[1:n.beta.params]
        gamma.hat <- o.fit$coefficients[-(1:n.beta.params)]
        exp.gamma.v <- as.vector(exp( cbind(icpt, dd$v.diff.pp) %*% gamma.hat ))
        
        e.model.matrix <- cbind(dd$y.1.x.sum.pp, icpt, dd$z.diff.pp)
        e.fit <- geeFit(dd$a.1.pp, e.model.matrix, "logit")
        beta2.hat <- e.fit$coefficients[1:n.beta.params]
        alpha.hat <- e.fit$coefficients[-(1:n.beta.params)]
        exp.alpha.z <- as.vector(exp( cbind(icpt, dd$z.diff.pp) %*% alpha.hat ))
        
    } else {
        
        coef.names <- c(object$ax.names,
                        object$ax.names,
                        object$v.names,
                        object$ax.names,
                        object$z.names)
        
        o.model.matrix <- cbind(dd$ax.diff.pp, dd$v.diff.pp)
        o.fit <- geeFit(dd$y.1.pp, o.model.matrix, "logit")
        beta1.hat <- o.fit$coefficients[1:n.beta.params]
        gamma.hat <- o.fit$coefficients[-(1:n.beta.params)]
        exp.gamma.v <- as.vector(exp( dd$v.diff.pp %*% gamma.hat -
                                      dd$x.2.pp %*% beta1.hat ))
        
        e.model.matrix <- cbind( dd$yx.diff.pp, dd$z.diff.pp )
        e.fit <- geeFit( dd$a.1.pp, e.model.matrix, "logit" )
        beta2.hat <- e.fit$coefficients[1:n.beta.params]
        alpha.hat <- e.fit$coefficients[-(1:n.beta.params)]
        exp.alpha.z <- as.vector(exp( dd$z.diff.pp %*% alpha.hat -
                                      dd$x.2.pp %*% beta2.hat ))

    }

    ## Estimating functions for the nuisance models
    U.onuis <- o.model.matrix * o.fit$res
    U.enuis <- e.model.matrix * e.fit$res

    ## Derivatives of the estimating functions for the nuisance models
    d.U.onuis.beta <- matrix(rep(0, ( n.beta.params + n.gamma.params ) *
                                    n.beta.params),
                             nrow = n.beta.params + n.gamma.params)
    d.U.onuis.beta1.gamma <- crossprod(o.model.matrix, o.fit$d.res)
    d.U.onuis.beta2.alpha <- matrix(rep(0, ( n.beta.params + n.gamma.params ) *
                                           ( n.beta.params + n.alpha.params)),
                                    nrow = n.beta.params + n.gamma.params)
    d.U.onuis <- cbind(d.U.onuis.beta,
                       d.U.onuis.beta1.gamma,
                       d.U.onuis.beta2.alpha )
    
    d.U.enuis.beta <- matrix(rep(0, ( n.beta.params + n.alpha.params ) *
                                    n.beta.params),
                             nrow = n.beta.params + n.alpha.params)
    d.U.enuis.beta1.gamma <- matrix(rep(0, ( n.beta.params + n.alpha.params ) *
                                           ( n.beta.params + n.gamma.params)),
                                    nrow = n.beta.params + n.alpha.params )
    d.U.enuis.beta2.alpha <- crossprod( e.model.matrix, e.fit$d.res )
    d.U.enuis <- cbind(d.U.enuis.beta,
                       d.U.enuis.beta1.gamma,
                       d.U.enuis.beta2.alpha )

    ## DR estimation of the log odds ratio parameter 
    all.args <- c(list(beta.init = beta1.hat,
                       eq.func = dr.logit.eq.func,
                       d.eq.func = NULL,
                       arg.list = list(y.f = dd$y.1.pp,
                                       ax.f = dd$a.1.x.sum.pp,
                                       a.f = dd$a.1.pp,
                                       x.f = dd$x.sum.pp,
                                       exp.gamma.v = exp.gamma.v,
                                       exp.alpha.z = exp.alpha.z)))    
    ## Call equation solver with beta.init as initial guess
    ## and eq.func as estimation function
    root.object <- try( do.call(rootFinder, all.args) )
    
    if (inherits(root.object, "try-error")) {

        coefficients <- c( rep(NA, n.beta.params),
                          beta2.hat,
                          alpha.hat,
                          beta1.hat,
                          gamma.hat )
        
        names(coefficients) <- coef.names
        U <- NULL
        d.U.sum <- NULL
        optim.object.o <- NULL
        optim.object.e <- NULL
        id <- NULL
        id.vcov <- NULL
        
    } else {
    
        beta.hat <- root.object$roots
        optim.object <- root.object$optim.object
    
        exp.beta.ax <- as.vector(exp( dd$a.1.x.sum.pp %*% beta.hat ) )
        exp.beta.x <- as.vector(exp(  dd$x.sum.pp %*% beta.hat ) )
    
        p <- as.vector(exp.beta.x * exp.alpha.z * (1 + exp.gamma.v) )
        q <- as.vector(p + 1 + exp.beta.x * exp.gamma.v)
        e.star <- as.vector(p / q)
    
        res.a.e.star <- as.vector(dd$a.1.pp - e.star)
        res.y <- as.vector(dd$y.1.pp - 1 + 1/(1 + exp.beta.ax * exp.gamma.v))

        ## DR estimating function
        U.dr <- dd$x.sum.pp * res.a.e.star * res.y

        ## derivatives of DR estimating function
        d.res.a.e.star.beta <- -dd$x.sum.pp * (e.star / q) *
            ( q - p - exp.beta.x * exp.gamma.v)
        d.res.y.beta <- -dd$a.1.x.sum.pp * exp.beta.ax *
            exp.gamma.v / (1 + exp.beta.ax * exp.gamma.v)^2
        d.U.dr.beta <- crossprod(dd$x.sum.pp ,
                                 d.res.a.e.star.beta * res.y +
                                 d.res.y.beta * res.a.e.star )

        d.U.dr.beta1 <- beta.zeros
        d.res.a.e.star.gamma <- -o.model.matrix[, -(1:n.beta.params)] *
            exp.gamma.v *
            (exp.beta.x / q) *
            (exp.alpha.z - e.star * (exp.alpha.z + 1))
        d.res.y.gamma <- -o.model.matrix[, -(1:n.beta.params)] *
            exp.beta.ax *
            exp.gamma.v / (1 + exp.beta.ax * exp.gamma.v)^2
        d.U.dr.gamma <- crossprod(dd$x.sum.pp, 
                                  d.res.a.e.star.gamma * res.y +
                                  d.res.y.gamma * res.a.e.star )    
        d.U.dr.beta2 <- beta.zeros
        d.res.a.e.star.alpha <- -e.model.matrix[, -(1:n.beta.params)] *
            e.star * (1 - e.star)
        d.U.dr.alpha <- crossprod(dd$x.sum.pp,
                                  d.res.a.e.star.alpha * res.y )
               
        d.U.dr <- cbind(d.U.dr.beta,
                        d.U.dr.beta1,
                        d.U.dr.gamma,
                        d.U.dr.beta2,
                        d.U.dr.alpha )

        U <- cbind(U.dr, U.onuis, U.enuis)
        d.U.sum <- rbind( d.U.dr, d.U.onuis, d.U.enuis)
        
        coefficients <- c(beta.hat,
                          beta2.hat,
                          alpha.hat,
                          beta1.hat,
                          gamma.hat)
        
        names(coefficients) <- coef.names
        
    }

    result <- list(coefficients = coefficients,
                   coef.names = coef.names, 
                   U = U,
                   d.U.sum = d.U.sum,
                   optim.object = optim.object,
                   optim.object.o = NULL,
                   optim.object.e = NULL,
                   id = dd$id.pp,
                   id.vcov = dd$id.pp,
                   clusters.info = clusters.info,
                   pp.info = pp.info)

    return(result)
        
}
