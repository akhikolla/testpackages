#' @title Character displacement likelihood ratio test
#' @description Conducts a likelihood ratio test between empirical data (phylogeny and trait data), and simumlations from the function chr.disp.sim using an approximate Bayesian computation (ABC) approach (Clarke et al. 2017)
#' @param emp.tree An empirical phylogeny - a object of class \code{phylo} (see \pkg{ape}).
#' @param emp.data Continuous trait data matrix
#' @param param.out simulated data from the function \code{chr.disp.sim}
#' @param posteriorSize The number of samples to use in the likelihood-ratio test
#' @return List containing element of 'estimates' with the estimates of sigma and a, with the Brownian motion (a = 0) summarised in column one and the character displacement (a > 0) in column two. 'likelihood' contains the likelihood of the Brownian motion model and the character displacement model, and the likelihood ratio test estimate. If used, there is an estimate of Blomberg's K for the empirical and simulated data.
#' @useDynLib motmot
#' @importFrom Rcpp sourceCpp
#' @import ks
#' @seealso \code{\link{chr.disp.sim}}, \code{\link{chr.disp.param}}
#' @references Clarke M, Thomas GH, Freckleton RP. 2017. Trait evolution in adaptive radiations: modelling and measuring interspecific competition on phylogenies. The American Naturalist. 189, 121-137.
#' @author Magnus Clarke and Mark Puttick
#' @examples
#' ## import finch data form Clarke et al. (2017)
#' data(finches)
#' ## simulate small amount of data 
#' ## (example only - many more datasets are required for accuracy)
#' param.simulation <- chr.disp.param(finch.tree, n.sim = 100, n.steps=100,
#' max.sigma = 8, max.a = 8, ntraits=1, 
#' allopatry=as.matrix(allopatric.data), mc.cores = 1)
#' chr.disp.lrt(finch.tree, finch.data, param.simulation, 50)
#' @export

chr.disp.lrt <- function(emp.tree, emp.data, param.out, posteriorSize=500) {

    sig <- param.out[[1]][ ,1]
    atry <- param.out[[1]][ ,2]
    sstat <- param.out[[1]][ ,-(1:2)]
    
    input.match <- param.out$input.arguments
    max.sigma <- input.match$max.sigma
	max.a <- input.match$max.a
	est.blomberg.k <- input.match$est.blomberg.k
	n.sim <- input.match$n.sim
	n.traits <- input.match$n.traits
	trait.lim <- input.match$trait.lim
	# sympatry <- param.out$sympatry
	# allopatry <- param.out$allopatry
	

    # Get summary stats for the true data, and distance to sims
    tstat <- summary_stats(phy=emp.tree, est.blomberg.k=est.blomberg.k, y=emp.data)
    diff <- colSums(abs(t(sstat)[-3,] - unlist(tstat[-3])) ^ 2)

    # Get simulations from nth closest to closest 
    h.1_post <- order(diff)[1:posteriorSize]
    u.sig <- sig[h.1_post]
    u.atry <- atry[h.1_post]
    h.1_post <- matrix(ncol=2, nrow=length(u.sig))
    h.1_post[,1] <- u.sig
    h.1_post[,2] <- u.atry

    k.out <- ks::kde(h.1_post, xmin=c(0, 0), xmax=c(max.sigma, max.a), binned=FALSE)
    k.0.out <- ks::kde(h.1_post, xmin=c(0, 0), xmax=c(max.sigma, 0), binned=FALSE)

    # Use kernel smoothing to estimate likelihood maxima with and without competition.
    k.max.index <- which(k.out$estimate == max(k.out$estimate), arr.ind = TRUE)
    h.1.lik <- k.out$estimate[k.max.index[1], k.max.index[2]]
    h.1.est <- c(unlist(k.out$eval.points)[k.max.index[1]], unlist(k.out$eval.points)[length(k.out$estimate[,1]) + k.max.index[2]])

    k.0.max.index <- which(k.0.out$estimate == max(k.0.out$estimate), arr.ind = TRUE)
    h.0.lik <- k.0.out$estimate[k.0.max.index[1], k.0.max.index[2]]
    h.0.est <- c(unlist(k.0.out$eval.points)[k.0.max.index[1]], unlist(k.0.out$eval.points)[length(k.0.out$estimate[,1]) + k.0.max.index[2]])
    
    likelihood.ratio.test <- -2 * log(h.0.lik / h.1.lik )
    p.value <- pchisq(likelihood.ratio.test, 1)
   
    output <- list()
    output$estimates <- data.frame(h.0.est, h.1.est)
    rownames(output$estimates) <- c("sigma", "a")
    colnames(output$estimates) <- c("Brownian motion", "Character displacement model")
    output$likelihood <- data.frame(log(h.0.lik), log(h.1.lik), likelihood.ratio.test, p.value)
	if(est.blomberg.k) output$blomberg.k <- c("empirical.blomberg.k"=unlist(tstat[3]), "simulated.mean.blomberg.k"=mean(param.out[[1]][,4]))
	return(output)
	
}