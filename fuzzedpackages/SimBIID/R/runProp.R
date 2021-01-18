## function to generate valid proposal
runProp <- function(i, t, priors, prevWeights, prevPars, propCov, tols, data, u, func, func_args) {
    accrate <- 0
    valid <- 0
    while(valid == 0) {
        if(t == 1) {
            pars <- rep(NA, nrow(priors))
            ## sample from prior
            for(j in 1:nrow(priors)) {
                pars[j] <- do.call(priors$dist[j], list(1, priors$p1[j], priors$p2[j]))
            }
        } else {
            ## sample from previous generation
            k <- sample(1:length(prevWeights), 1, prob = prevWeights)
            pars <- mvtnorm::rmvnorm(1, 
                            mean = prevPars[k, ], 
                            sigma = propCov
            )
            pars <- as.vector(pars)
        }
        
        ## check valid proposals
        chk <- 0
        for(j in 1:nrow(priors)) {
            chk <- chk + ifelse(do.call(priors$ddist[j], list(pars[j], priors$p1[j], priors$p2[j])) == 0, 1, 0)
        }
        if(chk == 0) {
            ## simulate from model
            fargs <- c(func_args, list(pars = pars, data = data, tols = tols, u = u))
            out <- do.call("func", fargs)
            if(!is.na(out[1])) {
                valid <- 1
            }
        }
        
        ## update counter
        accrate <- accrate + 1
    }
    
    if(t == 1) {
        weightsNew <- 1
    } else {
        ## calculate unnormalised weight
        weightsNew <- rep(NA, nrow(priors))
        for(j in 1:nrow(priors)) {
            weightsNew[j] <- do.call(priors$ddist[j], list(pars[j], priors$p1[j], priors$p2[j]))
        }
        weightsNew <- prod(weightsNew) / sum(prevWeights * apply(prevPars, 1, function(x, pars, propCov) {
            mvtnorm::dmvnorm(pars, mean = x, sigma = propCov)
        }, pars = pars, propCov = propCov))
    }
    list(pars = pars, out = out, weightsNew = weightsNew, accrate = accrate)
}      
