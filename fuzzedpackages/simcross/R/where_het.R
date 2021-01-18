## where_het.R

# where_het
#
#' Find heterozygous regions
#'
#' Find regions of heterozygosity in an individual
#'
#' @param ind An individual object, as output be
#' [create_parent()] or [cross()]
#'
#' @return A matrix with two columns; each row indicates the start and
#' end of a region where the individual is heterozygous
#'
#' @export
#' @seealso [sim_from_pedigree()],
#' [convert2geno()]
#' @examples
#' mom <- create_parent(100, 1:2)
#' dad <- create_parent(100, 1:2)
#' child <- cross(mom, dad)
#' where_het(child)
where_het <-
    function(ind)
{
    if(length(ind$mat$locations)==length(ind$pat$locations) &&
       all(ind$mat$locations == ind$pat$locations) &&
       all(ind$mat$alleles == ind$pat$alleles)) {
        het <- matrix(ncol=2, nrow=0)
        dimnames(het) <- list(NULL, c("left", "right"))
        return(het)
    }

    u <- c(0, sort(unique(c(ind$mat$locations,ind$pat$locations))))
    het <- NULL
    for(i in 2:length(u)) {
        toright <- which(ind$mat$locations >= u[i])
        mat <- ind$mat$alleles[toright[1]]

        toright <- which(ind$pat$locations >= u[i])
        pat <- ind$pat$alleles[toright[1]]

        if(mat!=pat) { # heterozygous
            if(is.null(het)) het <- cbind(u[i-1],u[i])
            else het <- rbind(het,c(u[i-1],u[i]))
        }
    }

    # clean up
    if(nrow(het) > 1) {
        keep <- rep(TRUE,nrow(het))
        for(j in 2:nrow(het)) {
            if(het[j,1] == het[j-1,2]) {
                het[j,1] <- het[j-1,1]
                keep[j-1] <- FALSE
            }
        }
        het <- het[keep,,drop=FALSE]
    }

    dimnames(het) <- list(1:nrow(het), c("left", "right"))
    het
}
