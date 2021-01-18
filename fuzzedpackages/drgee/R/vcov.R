## vcov.R

## Various functions to calculate the variance

## Calculate robust or cluster robust variance
## give residuals from estimating equations
## and the sum of their variance
robustVcov <- 
    function(U, d.U.sum, id = NULL) { 
 
        n.par <- ncol(d.U.sum) 
        n.obs <- nrow(U) 
 
        if( !is.null(id) ) { 
             
            ## If data is clustered we 
            ## divide by the number of clusters 
            n.clust <- length( unique(id) ) 
             
            ## and sum clusterwise 
            if(n.clust != n.obs){
                 
                U.dt <- as.data.table(cbind(id, U))
                
                ## We need unique variable names in the data.table
                names(U.dt) <- paste(c("id", 1:n.par))
                
                setkey(U.dt, id)
                
                U <- as.matrix( U.dt[, lapply(.SD, sum), by = id] )[, -1] 
            } 
             
        } else { 
             
            n.clust <- n.obs
             
        } 
 
        d.U <- d.U.sum / n.clust
        
        inv.d.U <- try( solve(d.U) )
 
        if (inherits(inv.d.U, "try-error")) { 
 
            vcov <- matrix(rep( NA, n.par^2 ), ncol = n.par ) 
 
        } else { 
             
            vcov <- inv.d.U %*% cov(U, U) %*% t(inv.d.U) / n.clust
             
        } 
         
        return( vcov ) 
    } 
 
## This is the old version
## only kept for backward compatibility
robVcov <-
    function(U, d.U, id = NULL) {

        n.obs <- nrow(U)

        n.par <- ncol(d.U)

        if(!is.null(id)) {
            n.clust <- length(unique(id))
        } else {
            n.clust <- 1
        }

        inv.d.U <- try(solve(d.U))

        if (inherits(inv.d.U, "try-error")) { 

            vcov <- matrix(rep( NA, n.par^2 ), ncol = n.par )

        } else {
            
            ## If data is not clustered we divide by
            ## the number of observations
            if (is.null(id) | n.clust == n.obs) {
                
                correction.term <- 1 / n.obs
                
            } else {
                
                ## If data is clustered we use the clustersums of u
                ## and divide by the number of clusters
                u.dt <- as.data.table(cbind(id, U))
                setkey(u.dt, id)
                U <- as.matrix( u.dt[, lapply(.SD, sum), by = id] )[, -1]

                correction.term <- n.clust / n.obs^2
            }

            vcov <- inv.d.U %*% cov(U, U) %*% t(inv.d.U) * correction.term
            
        }

        return( vcov )
    }

## A generic function to return the naive variance
## that is -solve( d.U ), which is correct
## for ML estimation for a fitted object
naiveVcov <- function(object){
    
    UseMethod("naiveVcov")
    
}

## A generic function to return the cluster robust variance.
## This function is useful when we want to obtain
## correct variance of the estimate for clustering
## when the cluster variable is not the same as the
## one used in the analysis
## for ML estimation
clusterRobustVcov <- function(object, clusterid = NULL){
    
    UseMethod("clusterRobustVcov")
    
}

