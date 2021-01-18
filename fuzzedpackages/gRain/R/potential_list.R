## Create potential list (rip, universe)
##
.mkArrayList <- function(rip.order, universe, values=1){
    cliques <- rip.order$cliques
    
    potlist  <- as.list(rep(NA, length(cliques)))
    
    for ( i in seq_along(cliques)){
        cq    <- cliques[[ i ]]
        vlab  <- universe$levels[cq]
        potlist[[ i ]] <- tabNew(cq, vlab, values)
    }
    potlist
}

.initArrayList <- function(x, values=NA){
    lapply(x, function(z) {
        z[] <- values             
        z
    } )
}

## Insert cpt's into potential list (cptlist, APlist)
##
.insertCPT <- function(cptlist, potlist, details=0)
{
    if (details>=1) cat(".Inserting cpt's in potential list [.insertCPT]\n")

    APnames <- lapply(potlist, function(x) names(dimnames(x)))
    CPnames <- unname(lapply(cptlist, function(x) varNames(x)))

    hosts    <- .get_hosts(CPnames, APnames)

    for (i in 1:length(cptlist)) {
            cptc <- cptlist[[ i ]]
            j    <- hosts[ i ]
            ##print(potlist[[j]])
            potlist[[ j ]] <- tableOp( potlist[[ j ]], cptc, "*" )
        }
    .infoPrint(details, 4, {cat("....potlist (after insertion):\n"); print(potlist) })
    potlist
}



## For each element (vector) x in xx.set, find the element (vector) y in
## yy.set such that x is contained in y
.get_hosts <- function(xx.set, yy.set){
    unlist(lapply(1:length(xx.set), function(i) which(is_inset(xx.set[[i]], yy.set, index=TRUE) > 0)[1]))
    ## Alternative:
    ## v <- lapply(xx.set, get_superset, yy.set, all=FALSE)
    ## v[lapply(v, length) == 0] <- NA
    ## unlist(v)
}






