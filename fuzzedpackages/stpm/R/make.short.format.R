#' An internal function which construct short data format from a given long
#' @param x Dataset
#' @param col.id Column ID index
#' @param col.status Column status index
#' @param col.t1 Column t1 index
#' @param col.t2 Column t2 index
#' @param col.cov Column covariates indices
#' @return column index(es) in the provided dataset
make.short.format <- function(x, col.id=1, col.status=2, col.t1=3, col.t2=4, col.cov=5)
{
    splitted <- split(x, x[ , col.id])
    
    for(i in 1:length(splitted))
    {
        dims <- dim(splitted[[i]])
        row.to.bind <- splitted[[i]][dims[1], ]
        row.to.bind[, c(col.id, col.status, col.t1, col.t2)] <- row.to.bind[, c(col.id, col.status,col.t2,col.t1)]
        row.to.bind[, seq(col.cov,dims[2],2)] <- row.to.bind[, seq(col.cov+1, dims[2], 2)]
        splitted[[i]][dims[1], col.status] <- 0 # Make it 0 even if its 1 because next row will be attached
        splitted[[i]] <- rbind(splitted[[i]], row.to.bind)
    }
    
    dd <- do.call("rbind", splitted)
    names(dd)[col.t1] <- "t"
    rownames(dd) <- c()
    dd <- dd[, c(col.id, col.status, col.t1, seq(col.cov, dims[2], 2))]
    return(dd)
}