## function to simulate outputs
runRef <- function(i, priors, func, func_args) {
    
    if(is.character(priors[, 1])){
        ## set up vector
        pars <- rep(NA, nrow(priors))
        
        ## sample from prior
        for(j in 1:nrow(priors)) {
            pars[j] <- do.call(priors$dist[j], list(1, priors$p1[j], priors$p2[j]))
        }
    } else {
        pars <- unlist(priors[i, ])
    }
    
    ## simulate from model
    fargs <- c(func_args, list(pars = pars))
    out <- do.call("func", fargs)
    
    ## return output
    list(pars = pars, out = out)
}      
