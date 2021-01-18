## function for calculating tolerances using bisect method
bisectTols <- function(out, data, tols, ptol, ptollim = 0.1) {
    
    ## select initial candidate tolerances
    diffs <- apply(out, 1, function(x, data) {
        abs(x - data)
    }, data = data)
    if(is.null(nrow(diffs))) {
        diffs <- matrix(diffs, ncol = 1)
    } else {
        diffs <- t(diffs)
    }
    tolcand <- apply(diffs, 2, function(x, tolprop) {
        stats::quantile(x, probs = tolprop)
    }, tolprop = ptol)
    
    if(length(tolcand) == 1) {
        return(tolcand)
    }
    
    ## now check the proportion of particles in current generation that 
    ## would have been accepted under candidate tolerances
    ind <- apply(diffs, 1, function(x, tolcand) {
        all(x <= tolcand)
    }, tolcand = tolcand)
    ind <- sum(ind) / length(ind)
    
    ## check whether this is with acceptable range of target proportion
    ind <- ifelse(ind > (ptol - ptollim), 1, 0)
    while(ind == 0) {
        ## if not, then conduct bisect method to generate new candidate tolerances
        tolcand1 <- tolcand
        tolcand <- tolcand + ((tols - tolcand) / 2)
        ## now check the proportion of particles in current generation that 
        ## would have been accepted under candidate tolerances
        ind <- apply(diffs, 1, function(x, tolcand) {
            all(x <= tolcand)
        }, tolcand = tolcand)
        ind <- sum(ind) / length(ind)
        if(abs(ptol - ind) <= ptollim) {
            ind <- 1
        } else {
            if((ptol - ind) < 0) {
                tols <- tolcand
                tolcand <- tolcand1
            }
            ind <- 0
        }
    }
    ## set new tolerances
    tolcand
}