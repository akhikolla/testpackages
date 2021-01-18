## conditFit.R

## Utility function to obtain the exact score residuals
## from clogit::survival
getScoreResidualsFromClogit <-
    function(
             fit,   ## The fitted object
             coefs, ## The estimated parameters
             y,     ## The observed outcomes
             x,     ## The observed covariates
             id     ## The cluster identification variable
             )
{
    if ( !missing( fit ) ) {
        coefs <- coef(fit)
    } else if ( missing( coefs ) ) {
        stop("\nEither or a fitted object or an estimate of the parameters needs to be supplied");
    }
    
    if ( length(coefs) != ncol(x) ) {
        stop("\nThe number of coefficients do not match the number of columns in the design matrix\n\n")
    }
    
    ## Extract info about clusters        
    clusters.info <- getDiscordantClustersInfo(y, id)
    disc.clusters <- with(clusters.info, clusters.info[which(disc), ])

    resids <- .Call(conditRes,
                    coefs, 
                    disc.clusters$ysum, 
                    disc.clusters$clust.size, 
                    disc.clusters$min.idx,
                    disc.clusters$inv,
                    as.vector(y),
                    x)

    U <- x * as.vector( resids$res )
    
    dU.sum <- crossprod( x , resids$dres )

    n.clust <- length( unique( id ) )
    
    return ( list(U = U, dU.sum = dU.sum, dU.mean = dU.sum / n.clust) )

}


getDiscordantClustersInfo <-
    function(y, ## The observed outcomes
             id ## The cluster identification variable
             )
{
    
    ## Assuming that the observations are sorted by id
    
    ## Find the total number of observations
    n.obs <- length(id)

    ## Indices in the original matrix
    idx <- seq_len(n.obs)

    ## Create data table with variables
    ## idx - row number in the original data vector
    ## id  - cluster indentification variable
    ## y   - the observec outcome
    y.dt <- data.table(idx = idx,
                       id = id,
                       y = y)

    names(y.dt) <- c("idx", "id", "y")

    setkeyv(y.dt, "id")

    ## Create variables
    ## clust.size - for cluster size
    ## ysum       - for cluster sum of y
    ## min.idx    - for first row number where cluster starts
    clusters <- data.table( with(y.dt,
                                 y.dt[, list(clust.size = .N,
                                             ysum = sum(y),
                                             min.idx = min(idx)), 
                                      by = id]) )

    ## Create variable
    ## disc - for outcome discordance
    with(clusters, clusters[, disc := TRUE])
    with(clusters, clusters[ysum == 0 | ysum == clust.size, disc := FALSE])
    
    ## Create variable
    ## inv - indicating whether more than half
    ##       of observations are cases
    with(clusters, clusters[, inv := 0L])
    with(clusters, clusters[ysum >  0.5 * clust.size, inv := 1L])

    ## The result is a data.table object with elements
    ## id         - cluster indentification variable
    ## clust.size - cluster sizes
    ## ysum       - cluster sums of y
    ## min.idx    - first row number where cluster starts
    ## disc       - for outcome discordance
    ## inv        - indicating whether more than half
    ##              of observations are cases

    return( clusters )
    
}

getResidsFromClogit <-
    function(
             estimate,              ## The estimated regression parameters
             clusters.info,         ## A data table containing info about the clusters
             y,                     ## The observed outcomes
             x,                     ## The observed covariates
             discordant.only = FALSE ## Use the discordant pairs only
             )
{
    ## Extract the needed variables
    if (discordant.only) {
        
        disc.clusters <- with(clusters.info, clusters.info[which(disc), ])

        resids <- .Call("conditRes",
                        estimate,
                        disc.clusters$ysum, 
                        disc.clusters$clust.size, 
                        disc.clusters$min.idx,
                        disc.clusters$inv,
                        as.vector(y),
                        x,
                        PACKAGE = "drgee")
        
    } else {

        resids <- .Call("conditRes",
                        estimate,
                        clusters.info$ysum, 
                        clusters.info$clust.size, 
                        clusters.info$min.idx,
                        clusters.info$inv,
                        as.vector(y),
                        x,
                        PACKAGE = "drgee")
        
    }

    return ( resids )

}

conditFit <-
    function(y,
             x,
             y.names = colnames(y),
             x.names = colnames(x), 
             id) {

        ## x.cent <- .Call("center", x, id, PACKAGE = "drgee")
        x.cent <- x
        colnames(x.cent) <- x.names

        ## Assuming that the observations are sorted by id
        
        ## Find the total number of observations
        n.obs <- nrow(x.cent)

        ## Find the number of parameters
        n.params <- ncol(x.cent) 
        
        ## #############################################
        ## Get info about discordant clusters
        ## #############################################

        clusters.info <- getDiscordantClustersInfo(y, id)
        disc.clusters <- with(clusters.info, clusters.info[which(disc), ])
        
        ## #############################################
        ## Create a data frame 
        ## #############################################
        
        yx.data <- data.frame( y, x.cent, id)

        ## Make sure the names are correct
        id.names <- names(yx.data)[length(names(yx.data))]
        names(yx.data) <- c(y.names, x.names, id.names)
        
        ## ##############################################
        ## Create a formula object to be used for clogit
        ## ##############################################
        
        oformula <- reformulate( c(x.names, "strata(id)"), y.names)

        ## ##############################################
        ## Fit using clogit
        ## ##############################################

        fit.clogit <- survival::clogit(oformula,
                                       data = yx.data, 
                                       method = "exact",
                                       subset = id %in% disc.clusters$id ) 
        
        beta.hat <- coef(fit.clogit)

        resids <- getResidsFromClogit(beta.hat,
                                      clusters.info,
                                      y,
                                      x.cent,
                                      discordant.only = TRUE)

        eq.res <- as.vector( resids$res )
        
        d.eq.res <- resids$dres

        naive.var <- vcov(fit.clogit)

        ## Use the old covariate names
        ## Will be different if the covariates contained interactions

        names(beta.hat) <- x.names
        colnames(d.eq.res) <- x.names
        dimnames(naive.var) <- list(x.names, x.names)

        return( list(coefficients = beta.hat,
                     res = eq.res,
                     d.res = d.eq.res,
                     eq.x = x.cent,
                     optim.object = NULL,
                     naive.var = naive.var
                     ) )
        
    }
