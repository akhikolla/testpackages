#' @title Bayesian MCMC for models of trait evolution
#' @description Fits Bayesian models for various models of continuous character evolution using a Metropolis-Hastings Markov Chain Monte Carlo (MCMC) approach
#' @param y A matrix of trait values.
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param model The model of trait evolution (see details).
#' @param mcmc.iteration Integer - the number of generations for which to run the MCMC chain
#' @param hiddenSpeciation Logical. If TRUE the psi model will include nodes that are on the 'full.phy' but not the tree pruned of trait data
#' @param useMean Logical. Use the branch-based estimates of extinction of mean (TRUE, default) for the "psi" and "multispi" models only applicable if "hiddenSpeciation" = TRUE
#' @param full.phy The full phylogeny containing the species that do not contain trait data so are not included in 'phy'
#' @param burn.in The proportion of the chain (as given by mcmc.iteration) which to discard as 'burn-in'
#' @param lowerBound Minimum value for parameter estimates
#' @param upperBound Maximum value for parameter estimates
#' @param random.start Use a random starting value for the MCMC run (TRUE), or use the environment set.seed() value
#' @param meserr A vector (or matrix) of measurement error for each tip. This is only applicable to univariate analyses
#' @param covPIC Logical. For multivariate analyses, allow for co-variance between traits rates (TRUE) or no covariance in trait rates (FALSE). If FALSE, only the trait variances not co-variances are used.
#' @param lambdaEst Logical.  Estimate lambda alongside parameter estimates to reduce data noise. Only applicable for models "kappa", "delta", "OU", "psi", "multispi", and "ACDC". Default=FALSE.
#' @param nodeIDs Integer - ancestral nodes of clades applicable to rate heterogenous and nested models of evolution (see details)
#' @param branchLabels Branches on which different psi parameters are estimated in the "multipsi" model
#' @param acdcScalar Logical.  For nested EB rate model, simultaneously estimated a rate scalar alongside EB model. Default=FALSE.
#' @param sample.every Number specifying the every nth that is sampled in the MCMC chain (default = 1).
#' @details The method estimates posterior probabilities using a Metropolis-Hastings MCMC approach that places a prior bounded uniform distribution on all parameters with an independence sampler. These prior distributions can be altered by changing the upperBound and lowerBound arguments. The MCMC model will estimate the posterior probability for the following models:
#' \itemize{
#' \item {model="kappa"} {fits Pagel's kappa by raising all branch lengths to the power kappa. As kappa approaches zero, trait change becomes focused at branching events. For complete phylogenies, if kappa approaches zero this infers speciational trait change. Default bounds from ~0 - 1.}
#' \item {model="lambda"} {fits Pagel's lambda to estimate phylogenetic signal by multiplying all internal branches of the tree by lambda, leaving tip branches as their original length (root to tip distances are unchanged). Default bounds from ~0 - 1.}
#' \item {model="delta"} {fits Pagel's delta by raising all node depths to the power delta. If delta <1, trait evolution is concentrated early in the tree whereas if delta >1 trait evolution is concentrated towards the tips. Values of delta above one can be difficult to fit reliably. Default bounds from ~0 - 5.}
#' \item {model="OU"} {fits an Ornstein-Uhlenbeck model - a random walk with a central tendency proportional to alpha. High values of alpha can be interpreted as evidence of evolutionary constraints, stabilising selection or weak phylogenetic signal. It is often difficult to distinguish among these possibilities. Default bounds from ~0 - 10.}
#' \item {model="psi"} {fits a acceleration-deacceleration model to assess to the relative contributions of speciation and gradual evolution to a trait's evolutionary rate (Ingram 2010).}
#' \item {model="ACDC"} {fits a model to in which rates can exponentially increased or decrease through time (Blomberg et al. 2003). If the upper bound is < 0, the model is equivalent to the 'Early Burst' model of Harmon et al. 2010. Default rate parameter bounds from ln(1e-10) ~ ln(20) divided by the root age.}
#' }
#' @return median The median estimate of the posterior for the parameter 
#' @return 95.HPD The 95 percent Highest Posterior Density for the parameter
#' @return ESS Effective Sample Size for the posterior
#' @return acceptance.rate The ratio for which new proposals were accepted during the MCMC chain
#' @return mcmc.chain Full MCMC chain containing all iterations (including burn-in)
#' @author Mark Puttick, Gavin Thomas
#' @seealso \code{\link{transformPhylo.ML}}, \code{\link{transformPhylo.ll}}, \code{\link{transformPhylo}}
#' @import coda
#' @examples
#' data(anolis.tree)
#' data(anolis.data)
#' attach(anolis.data)
#' male.length <- matrix(Male_SVL, dimnames=list(rownames(anolis.data)))
#' sortedData <- sortTraitData(anolis.tree, male.length)
#' phy <- sortedData$phy
#' male.length <- sortedData$trait
#' phy.clade <- extract.clade(phy, 182)
#' male.length.clade <- as.matrix(male.length[match(phy.clade$tip.label, rownames(male.length)),])
#' ## please note, this model will be need to run for longer to achieve convergence
#' lambda.mcmc <- transformPhylo.MCMC(y=male.length.clade, phy=phy.clade, 
#' model="lambda", mcmc.iteration=100, burn.in=0.1)
#' @export 

transformPhylo.MCMC <- function(y, phy, model, mcmc.iteration=1000, burn.in=0.1, hiddenSpeciation = FALSE, full.phy=NULL, lowerBound = NULL, upperBound = NULL, useMean = FALSE, random.start=TRUE, meserr = NULL, covPIC=TRUE, lambdaEst=FALSE, nodeIDs = NULL, branchLabels = NULL, acdcScalar = FALSE, sample.every=10) {

	if(length(model) > 1) stop("please provide one model only")
	if(ncol(y) > 1 && !is.null(meserr)) stop("meserr only applicable to univariate analyses")
	
	controlList=c(fnscale=-1, maxit=100, factr=1e-7, pgtol=0, type=2, lmm=5)
	
	model <- tolower(model)
	all.models <- c("lambda", "delta", "kappa", "ou", "acdc", "psi")
	if(is.na((match(model, all.models)))) stop(paste(model, "not recognised - please provide one of", paste0(all.models, collapse=", ")))
	if(random.start) set.seed(as.numeric(Sys.time()))

	bounds <- matrix(c(1e-08, 1, 1e-08, 1, 1e-08, 5, 1e-08, 20, 0, 1, 1e-08, 10, 1e-10, 100000), 7, 2, byrow = TRUE)
	rownames(bounds) <- c("kappa", "lambda", "delta", "alpha", "psi", "rate", "acdcRate")

	if(model == "lambda") {
		if (is.null(lowerBound)) {
			lowerBound <- bounds["lambda", 1]
			}
		if (is.null(upperBound)) {
			upperBound <- bounds["lambda", 2]
		}
		
		input.value <- runif(1, lowerBound, upperBound)
		lik.model <- function(pram) {
			lambda.phy <- transformPhylo(phy, model="lambda", y=y, lambda=pram, meserr=meserr)
			return(likTraitPhylo(y, lambda.phy, covPIC=covPIC)[[2]])
		}
		name.param <- c("Lambda")
	}

	if(model == "delta") {
		if (is.null(nodeIDs)) {
        		nodeIDs <- Ntip(phy) + 1
    		} else {
        		nodeIDs <- nodeIDs
    			}
		if (is.null(lowerBound)) {
			lowerBound <- bounds["delta", 1]
			}
		if (is.null(upperBound)) {
			upperBound <- bounds["delta", 2]
			}
		if (lambdaEst) {
    			lowerBound <- c(lowerBound, bounds["lambda", 1])
    			upperBound <- c(upperBound, bounds["lambda", 2])
    		    	}
		
		input.value <- sapply(1:length(lowerBound), function(x) runif(1, lowerBound[x], upperBound[x]))
		if(length(input.value) == 2) input.value[2] = 1
		lik.model <- function(pram) {
			if (lambdaEst) lambda <- tail(pram, 1) else lambda <- 1
        		delta.est <- pram[1]
        		lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
			delta.phy <- transformPhylo(lambdaPhy, model="delta", y=y, delta=delta.est, meserr=meserr, nodeIDs = nodeIDs)
			return(likTraitPhylo(y, delta.phy, covPIC=covPIC)[[2]])
			}
		
		if (lambdaEst) {
				name.param <- c("Delta", "Lambda")
			} else {
				name.param <- c("Delta")
			}
			
		}

		if(model == "kappa") {
			if (is.null(nodeIDs)) {
        			nodeIDs <- Ntip(phy) + 1
      		} else {
        			nodeIDs <- nodeIDs
      			}
			if (is.null(lowerBound)) {
				lowerBound <- bounds["kappa", 1]
				}
			if (is.null(upperBound)) {
				upperBound <- bounds["kappa", 2]
				}
		
			if (lambdaEst) {
       			lowerBound <- c(lowerBound, bounds["lambda", 1])
       			upperBound <- c(upperBound, bounds["lambda", 2])
        			}
		
			input.value <- sapply(1:length(lowerBound), function(x) runif(1, lowerBound[x], upperBound[x]))
					
			lik.model <- function(pram) {
				if (lambdaEst) lambda <- tail(pram, 1) else lambda <- 1
         		kappa.est <- pram[1]
         		lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
				kappa.phy <- transformPhylo(phy, model="kappa", y=y, kappa=kappa.est, meserr=meserr, nodeIDs = nodeIDs)
				return(likTraitPhylo(y, kappa.phy, covPIC=covPIC)[[2]])
				}
		
			if (lambdaEst) {
				name.param <- c("Kappa", "Lambda")
			} else {
				name.param <- c("Kappa")
			}
		}
	
		if(model == "ou") {
			if (is.null(nodeIDs)) {
        			nodeIDs <- Ntip(phy) + 1
      		} else {
        			nodeIDs <- nodeIDs
      			}
			if (is.null(lowerBound)) {
				lowerBound <- bounds["alpha", 1]
				}
			if (is.null(upperBound)) {
				upperBound <- bounds["alpha", 2]
				}
			input.value <- runif(1, lowerBound , upperBound)
			if(!is.ultrametric(phy)) {
				stop("non-ultrametric OU model not available in transformPhylo.MCMC at current, sorry")
				}
			if (lambdaEst) {
       			lowerBound <- c(lowerBound, bounds["lambda", 1])
       			upperBound <- c(upperBound, bounds["lambda", 2])
        		}
			input.value <- upperBound[1] + 1
			while(!(input.value <= upperBound[1] && input.value >= lowerBound[1])) input.value <- rexp(1, upperBound[1])
			if(lambdaEst) input.value[2] <- runif(1, lowerBound[2], upperBound[2])
		
			lik.model <- function(pram) {
				if (lambdaEst) lambda <- tail(pram, 1) else lambda <- 1
				ou.est <- pram[1]
        			lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
        			alpha.phy <- transformPhylo(phy, model="ou", y=y, alpha=ou.est, meserr=meserr)
				return(likTraitPhylo(y, alpha.phy, covPIC=covPIC)[[2]])
				}
	
		if (lambdaEst) {
				name.param <- c("alpha", "Lambda")
			} else {
				name.param <- c("alpha")
				}	
		}

		if(model == "acdc") {
			if (is.null(nodeIDs)) {
        			nodeIDs <- Ntip(phy) + 1
      		} else {
        			nodeIDs <- nodeIDs
      			}
			rootBranchingTime <- nodeTimes(phy)[1,1]
			if (is.null(lowerBound)) {
				lowerBound <- log(bounds["acdcRate", 1]) / rootBranchingTime
				}
			if (is.null(upperBound)) {
				upperBound <- log(bounds["acdcRate", 2]) / rootBranchingTime
				}
			if (acdcScalar) {
				upperBound[1] <- -1e-6
      			if(is.na(lowerBound[2])) {
      				lowerBound[2] <- 1
      				}
      			if(is.na(upperBound[2])) {
      				upperBound[2] <- 5
      			}
       		}
			if (lambdaEst) {
       			lowerBound <- c(lowerBound, bounds["lambda", 1])
       			upperBound <- c(upperBound, bounds["lambda", 2])
        			} 
			input.value <- sapply(1:length(lowerBound), function(x) runif(1, lowerBound[x], upperBound[x]))
			if (acdcScalar) {
        			lik.model <- function(pram) {
        				if (lambdaEst) lambda <- tail(pram, 1) else lambda <- 1
          			acdc.est <- pram[1]
          			scalarRate <- pram[2]
          			lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
          			acdc.phy <- transformPhylo(lambdaPhy, model="acdc", acdcRate=acdc.est, meserr=meserr, nodeIDs = nodeIDs, cladeRates = scalarRate)
          			return(likTraitPhylo(y, acdc.phy, covPIC=covPIC)[[2]])
        				}
				} else {
					lik.model <- function(pram) {
					if (lambdaEst) lambda <- tail(pram, 1) else lambda <- 1
         			acdc.est <- pram[1]
         			scalarRate <- 1
         			lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
         			acdc.phy <- transformPhylo(lambdaPhy, model="acdc", y=y, acdcRate=acdc.est, meserr=meserr)
					return(likTraitPhylo(y, acdc.phy, covPIC=covPIC)[[2]])
					}
				}
		
			if (lambdaEst) {
				name.param <- c("ACDC.rate", "Lambda")
				if (acdcScalar) name.param <- c("ACDC.rate", "ACDC.scalar", "Lambda")
			} else {
				name.param <- c("ACDC.rate")
				if (acdcScalar) name.param <- c("ACDC.rate", "ACDC.scalar")
				}
		}

		if(model == "psi") {
			if (is.null(nodeIDs)) {
    		    		nodeIDs <- Ntip(phy) + 1
    		  	} else {
        			nodeIDs <- nodeIDs
      			}
			if (is.null(lowerBound)) {
				lowerBound <- bounds["psi", 1]
				}
			if (is.null(upperBound)) {
				upperBound <- bounds["psi", 2]
				}
			if (lambdaEst) {
       			lowerBound <- c(lowerBound, bounds["lambda", 1])
       			upperBound <- c(upperBound, bounds["lambda", 2])
        			}	
			if (hiddenSpeciation) {
				if (is.null(full.phy)) stop("please provide a full phylogeny")
				full.data.match <- match(full.phy$tip.label, rownames(y))
				tips.no.data <- full.phy$tip.label[which(is.na(full.data.match))]
				phy <- dropTipPartial(full.phy, tips.no.data)
				}
	 		if(is.ultrametric(phy)) {
     	 		phy.bd <- birthdeath(phy)
     		 } else {
     		 	phy.bd <- birthdeath_motmot(phy)
     		 	}
     	 
      		mu_over_lambda <- phy.bd[[4]][1]
      		lambda_minus_mu <- phy.bd[[4]][2]
      		lambda.sp <- as.numeric(lambda_minus_mu / (1 - mu_over_lambda))
      		mu.ext <- as.numeric(lambda_minus_mu / (1 / mu_over_lambda - 1))			
			if (mu.ext > 0) {
				phy <- sampleHiddenSp(phy, lambda.sp = lambda.sp, mu.ext = mu.ext, useMean=useMean)
			} else {
				phy$hidden.speciation <- NULL
			}
			
			input.value <- sapply(1:length(lowerBound), function(x) runif(1, lowerBound[x], upperBound[x]))

			lik.model <- function(pram) {
				if (lambdaEst) lambda <- tail(pram, 1) else lambda <- 1
         		psi.est <- pram[1]
         		lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
         		psi.phy <- transformPhylo(phy, model="psi", y=y, psi=psi.est, meserr=meserr, lambda.sp = lambda.sp, nodeIDs = nodeIDs)
				return(likTraitPhylo(y, psi.phy, covPIC=covPIC)[[2]])
				}

			if (lambdaEst) {
				name.param <- c("Psi", "Lambda")
			} else {
				name.param <- c("Psi")
			}
	}
		
		if(model == "multispi") {
   		   if (is.null(branchLabels)) stop("for 'multipsi' model must provide branchLabels giving state for each branch")
      		states <- levels(factor(branchLabels))
      		
      		if (is.null(lowerBound)) {
        			lowerBound <- bounds[rep("psi", length(states)), 1]
        			if (lambdaEst) lowerBound[length(states) + 1] <- bounds["lambda", 1]
      			}
      		if (is.null(upperBound)) {
        			upperBound <- bounds[rep("psi", length(states)), 2]
        			if (lambdaEst) upperBound[length(states) + 1] <- bounds["lambda", 2]
      			}
      		
      		input.value <- sapply(1:length(lowerBound), function(x) runif(1, lowerBound[x], upperBound[x]))

      		if (hiddenSpeciation) {
        			if (is.null(full.phy)) stop("please provide a full phylogeny")
        			full.data.match <- match(full.phy$tip.label, rownames(y))
        			tips.no.data <- full.phy$tip.label[which(is.na(full.data.match))]
        			phy <- dropTipPartial(full.phy, tips.no.data)
      			}
			
			if(is.ultrametric(phy)) {
	      		phy.bd <- birthdeath(phy)
	      	} else {
	      		phy.bd <- birthdeath_motmot(phy)
	      		}
      	
      		mu_over_lambda <- phy.bd[[4]][1]
      		lambda_minus_mu <- phy.bd[[4]][2]
      		lambda.sp <- as.numeric(lambda_minus_mu / (1 - mu_over_lambda))
      		mu.ext <- as.numeric(lambda_minus_mu / (1 / mu_over_lambda - 1))

      		if (mu.ext > 0) {
      			phy <- sampleHiddenSp(phy, lambda.sp = lambda.sp, mu.ext = mu.ext, useMean = useMean)
      		} else {
      			phy$hidden.speciation <- NULL
      			}
		
			var.funmultipsi <- function(param) {
        			all.param <- length(param)
        			if (lambdaEst) {
        				lambda <- param[all.param]
          			psi <- param[-all.param]
        			} else {
        		  		psi <- param
        		  		lambda <- 1
        			}
        			lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
        			return(transformPhylo.ll(y = y, phy = lambdaPhy, branchLabels = branchLabels, psi = psi, model = "multipsi", meserr = meserr, covPIC = covPIC, lambda.sp = lambda.sp)[[2]])
      		}
      	
      		vo <- optim(start, var.funmultipsi, method = "L-BFGS-B", lower = lowerBound, upper = upperBound, control = controlList)

			lik.model <- function(pram) {
				all.param <- length(pram)
				if (lambdaEst) {
		 			lambda <- pram[all.param]
		 			psi <- pram[-all.param]
		 		} else {
		 			psi <- pram
		 			lambda <- 1
		 			}
       			lambdaPhy <- transformPhylo(y = y, phy = phy, lambda = lambda, model = "lambda", meserr = meserr)
       			psi.phy <- transformPhylo(phy = lambdaPhy, branchLabels = branchLabels, psi = psi, model = "multipsi", meserr = meserr, lambda.sp = lambda.sp)
				return(likTraitPhylo(y, psi.phy, covPIC=covPIC)[[2]])
				}
	
			if (lambdaEst) {
				name.param <- c(rep("multipsi", length(states)), "Lambda")
				} else {
				name.param <- rep("multipsi", length(states))
				}
		}
	
	
		prior.calc <- function(pram) {
			prior.uni <- sapply(1:length(pram), function(x) dunif(pram[x], lowerBound[x], upperBound[x]))
			return(sum(prior.uni))
			}

		model.posterior <- function(pram) return (lik.model(pram) + prior.calc(pram))	
	
		motmot.mcmc <- function(input.value, iterations, silent=FALSE) {
			propose.mcmc <- function(pram) {
				return(runif(length(pram), min=lowerBound, max=upperBound))
				}
			
		mcmc.chain <- matrix(input.value, nrow=1)
    			for (i in 1:iterations) {
    				proposed.move <- propose.mcmc(mcmc.chain[i, ])
     	   		chain.prob <- exp(model.posterior(proposed.move) - model.posterior(mcmc.chain[i,]))
     	   		if (runif(1) < chain.prob) {
     	   			mcmc.chain <- rbind(mcmc.chain, proposed.move)
        			} else {
     	   			mcmc.chain <- rbind(mcmc.chain, mcmc.chain[i,])
     	   		}
			if(!silent) {
				cat("\r", paste0("MCMC progress (%): ", sprintf("%.4f", i/iterations * 100)))
				if(i/iterations == 1) {
					cat("\n", " ")
					cat("\r", "finished MCMC")
					cat("\r", "")
					}
				}	
    			}
    			if(!silent) cat("\n", " ")
    			return(mcmc.chain)
		}
	
		chain <- motmot.mcmc(input.value, iterations=mcmc.iteration)
		burnIn <- ceiling(mcmc.iteration * burn.in)
		post.burn.in <- as.matrix(chain[-c(1:burnIn), ])
		post.burn.in <- as.matrix(post.burn.in[seq(1, nrow(post.burn.in), by=sample.every), ])
		rownames(post.burn.in) <- NULL
		
		acceptance.1 <- 1 - mean(apply(as.matrix(post.burn.in), 2, function(x) duplicated(signif(x, 7))))
		
		if(dim(chain)[2] == 1) {
			ess.mcmc <- coda::effectiveSize(post.burn.in)
			median.mcmc <- median(post.burn.in)
			hpd.mcmc <- coda::HPDinterval(as.mcmc(post.burn.in))
			hpd.mcmc <- as.numeric(hpd.mcmc)
			names(ess.mcmc) <- names(median.mcmc) <- name.param
			names(hpd.mcmc) <- c("lower 95% HPD", "upper 95% HPD")
		} else {
			ess.mcmc <- apply(post.burn.in, 2, coda::effectiveSize)
			median.mcmc <- apply(chain[-c(1:burnIn), ], 2, median)
			hpd.mcmc <- apply(chain[-c(1:burnIn), ], 2, function(x) {
				class(x) <- "mcmc"
				as.numeric(coda::HPDinterval(x))
				}
			)
			names(ess.mcmc) <- names(median.mcmc) <- name.param
			colnames(hpd.mcmc) <- name.param
		}

		cat("\n")
		output.mcmc <- list(median.mcmc, hpd.mcmc, ess.mcmc, acceptance.1, post.burn.in)
		names(output.mcmc) <- c("median", "95.HPD", "ESS", "acceptance.rate", "mcmc.chain")
		class(output.mcmc) <- "motmot.mcmc"
		print(output.mcmc[1:4])
		invisible(output.mcmc)
	}