#' @title Simulate character displacement data wrapper
#' @description Simulates phylogenetic trait data under a character displacement model (Clarke et al. 2017) in which traits interact inter-specifically, with competition between sympatric lineages driving trait values apart
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param n.steps Number of time steps the for the simulation (default = 1000 time steps).
#' @param n.sim Number of replications to simulate data
#' @param max.sigma The maximum value of Brownian variance in the simulation sampled from a U(0, max.sigma) distribution for each iteration
#' @param max.a The maximum value of the strength of competition between inter-specific lineages sampled from a U(0, max.a) distribution for each iteration
#' @param est.blomberg.k Logical. If TRUE, Blomberg's K is simultaneously estimated
#' @param ntraits Number of traits to be simulated
#' @param sympatry an optional matrix giving the time that each pair of species starts to interact
#' @param allopatry an optional matrix giving the times when species stop interacting
#' @param trait.lim an optional parameter that puts limits on the available trait-space, preventing trait values with magnitude greater than the value of lim
#' @param mc.cores Numeric. The number of parallel cores to be used in simulations. Only applicable on Linux and Mac systems
#' @return List containing the simulated data 'simulated.param': a matrix with each row represented an iteration, the sigma (Brownian variance) used in the iteration, the 'a' value used in each iteration, the mean and standard deviation between neighbouring trait values. The 'input.arguments' from the model, the 'input.phy' from the model, and the input 'sympatry' and 'allopatry' matrices. 
#' @useDynLib motmot
#' @importFrom Rcpp sourceCpp
#' @seealso \code{\link{chr.disp.sim}}, \code{\link{chr.disp.lrt}}
#' @references Clarke M, Thomas GH, & Freckleton RP. 2017. Trait evolution in adaptive radiations: modeling and measuring interspecific competition on phylogenies. The American Naturalist 189, 121-137.
#' @author Magnus Clarke and Mark Puttick
#' @examples
#' ## import finch data form Clarke et al. (2017)
#' data(finches)
#' ## simulate small amount of data 
#' ## (example only - many more datasets are required for accuracy)
#' param.simulation <- chr.disp.param(finch.tree, n.sim = 3, n.steps=100,  
#' max.sigma = 8, max.a = 8, ntraits=1, 
#' allopatry=as.matrix(allopatric.data), mc.cores = 1)
#' @export

chr.disp.param <- function (phy, n.sim = 100, n.steps=1000, max.sigma = 8, max.a = 8, est.blomberg.k = FALSE, ntraits=1, sympatry=NA, allopatry=NA, trait.lim=NA, mc.cores = 1) {
	
    sig <- runif(n.sim, 0, max.sigma)
    atry <- runif(n.sim, 0, max.a)
    if(!all(is.na(allopatry))) allopatry <- as.matrix(allopatry)
	if(!all(is.na(sympatry))) sympatry <- as.matrix(sympatry)    
    param.spl <- parallel::mclapply(1:n.sim, mc.cores = mc.cores, function(i) {
        d.out <- chr.disp.sim(phy = phy, n.steps=n.steps, a = atry[i], sigma = sig[i], ntraits=ntraits, sympatry=sympatry, allopatry=allopatry, trait.lim=trait.lim)$tval
        stats.sim <- summary_stats(phy = phy, y = as.matrix(d.out), 
            est.blomberg.k = est.blomberg.k)
        c(sig[i], atry[i], unlist(stats.sim))
    })
    
    if(!est.blomberg.k) {
    	    n.col <- 4
     	output.prm <- matrix(unlist(param.spl), ncol=n.col, byrow=TRUE)
    		colnames(output.prm) <- c("sigma", "a", "mean.gap", "sd.gap")
    	} else {
    		n.col <- 5
     	output.prm <- matrix(unlist(param.spl), ncol=n.col, byrow=TRUE)
    		colnames(output.prm) <- c("sigma", "a", "mean.gap", "sd.gap", "blomberg.k")
    	}
    input.pram <- data.frame(max.sigma = 8, max.a = 8, est.blomberg.k = as.logical(est.blomberg.k), n.sim = n.sim, n.traits=ntraits, trait.lim=trait.lim)
    return(list(simulated.param = output.prm, input.arguments = input.pram, 
        input.phy = phy, sympatry=sympatry, allopatry=allopatry))
}