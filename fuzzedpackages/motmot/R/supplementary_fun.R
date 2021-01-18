# allCladeMembers (internal function)
#
# This is an internal function to generate an \code{allCladeMembersMatrix}
# phy phylogeny in \pkg{ape} \code{phylo} format
# cladeMembersMatrix
# Gavin Thomas

allCladeMembers <- function(phy) {
	
	nodeIDs <- c(1:Ntip(phy), (Ntip(phy)+2):(Ntip(phy) + Nnode(phy)))
	k <- length(nodeIDs)
   cladeMembersMatrix <-  sapply(nodeIDs, function(k) {
    	nodeShiftID <- c(k, node.descendents(x = k , phy = phy))
    	as.numeric(phy$edge[, 2] %in% nodeShiftID)
    	}
    )
    return(cladeMembersMatrix)
}

# Identify branches (including tips) descended from a node (internal function).
# description Internal function to get presence absence of descendent branches from a vector of node numbers. The descendents include the branch leading to the focal node (i.e. node defines the stem group no crown group
# phy An object of class \code{phylo} (see \pkg{ape}).
# nodeIDs Vector of node numbers (positive integers).
# cladeMembersObj Matrix of clade membership
# The function returns a matrix of unique presences given the selected node. If the selected nodes are nested then presences are only recorded for the least inclusive node.
# return matrix Matrix of unique presences for each node id
# author Gavin Thomas
# examples
# ## Read in phylogeny and data from Thomas et al. (2009)
# data(anolis.tree)
# data(anolis.data)
# cladeIdentityMatrix <- cladeIdentity(phy=anolis.tree, nodeIDs=170)

cladeIdentity <- function (phy, nodeIDs, cladeMembersObj=NULL) 
{
	
	k <- length(nodeIDs)
	
    if(is.null(cladeMembersObj)) {
		cladeMembers <- matrix(NA, ncol = k, nrow = length(phy$edge[, 
														   1]))
		for (i in 1:k) {
			nodeShiftID <- c(nodeIDs[i], node.descendents(x = nodeIDs[i], phy = phy))
			cladeMembers[, i] <- as.numeric(phy$edge[, 2] %in% nodeShiftID)
		}
	}
	
	if (is.null(cladeMembersObj)==FALSE) {
		allNodes <- c(1:Ntip(phy), (Ntip(phy)+2):(length(phy$edge.length)+1))   #######
		cladeMembers <- as.matrix(cladeMembersObj[,match(nodeIDs, allNodes)])
	}
	
    originalOrder <- colSums(cladeMembers)
    richnessOrder <- sort(originalOrder, decreasing = FALSE, 
						  index.return = TRUE)
    cladeMembersOrdered <- matrix(cladeMembers[, richnessOrder$ix], 
								  ncol = length(nodeIDs))
	
	if (k>1) {
		for (i in 2:k){
			if (i ==2 ) {	cladeMembersOrdered[,i] <- cladeMembersOrdered[,i] - cladeMembersOrdered[,1:(i-1)]}
			else {cladeMembersOrdered[,i] <- cladeMembersOrdered[,i] - rowSums(cladeMembersOrdered[,1:(i-1)])}
		}
	}
	
    cladeMembers <- cladeMembersOrdered[, sort(richnessOrder$ix, index.return = TRUE)$ix]
    cladeMembers <- matrix(cladeMembers, ncol = length(nodeIDs))
    return(cladeMembers)
}


# Create design matrix (internal function)
# This is an internal function to generate the design matrix required to define different means for each hypothesised rate category.
# y The response variable - typically a continuous trait.
# x The explanatory (discrete) variable used to define the hypothesised rate categories. Can be specified as a column number or column name.
# data A data frame containing (minimally) the x and y variables as columns with species names as rownames.
# common.mean a logical specififying whether each rate category should have its own mean (\code{common.mean=FALSE}) or all categories should have the same mean (\code{common.mean=FALSE}). See Thomas et al. (2009) for a discussion on the impact of assumptions about mean on rate estimates..
# author Gavin Thomas

make.anc <-
function(y, x, data=NULL, common.mean=FALSE) {
	
	if(is.factor(x) == "FALSE"){
		stop("The discrete trait must be a factor")	
	}
	
	m <- model.frame(y ~ x, data)
	x <- model.matrix(y ~ x, m)
	return(x)
}


#  Name check
# This is an internal function to check if names in a trait vector or matrix match a tree
# phy a phylogeny in \pkg{ape} phylo format
# data A matrix or vector of named trait data corresponding to species in phy
# Gavin Thomas

name.check <- function (phy, data) 
{
    if (is.vector(data)) {
        data.names <- names(data)
    }
    else {
        data.names <- rownames(data)
    }
    t <- phy$tip.label
    r1 <- t[is.na(match(t, data.names))]
    r2 <- data.names[is.na(match(data.names, t))]
    r <- list(sort(r1), sort(r2))
    names(r) <- cbind("tree_not_data", "data_not_tree")
    if (length(r1) == 0 && length(r2) == 0) 
	return("OK")
    else return(r)
}

# Phylogenetically independent contrasts (internal)
# Calculates phylogenetically independent contrasts.
# param x A matrix of trait values with taxon names as rownames.
# param phy An object of class \code{phylo} (see \pkg{ape}).
# details Extracts values for contrasts, expected variances using contrasts by calling \code{pic}.
# return contr A matrix with two columns containing raw contrasts in the first column and their expected variances in the second column.
# return root.v Expected variances of branches either side of root
# return V Expected variance at the root 
# references Felsenstein J. 1985. Phylogenies and the comparative method. American Naturalist 125, 1-15.
# author Gavin Thomas, Rob Freckleton, Emmanuel Paradis

pic.motmot <- function (x, phy) {
out <- pic(x, phy, rescaled.tree=TRUE, var.contrasts=TRUE, scaled=FALSE)
   contr <- out[[1]]	
	nb.tip <- length(phy$tip.label)
	idx <- which(out[[2]]$edge[,1] == (nb.tip+1) )
	
	root.v <- out[[2]]$edge.length[idx]
	V <- prod(root.v)/(sum(root.v))
  
  return(list(contr = contr, root.v = root.v, V = V))
}

# remove species occuring before time in the past (internal function) 
# removes tips and lineages after a time in the past
# param phy An object of class \code{phylo} (see \pkg{ape}).
# param traitData data associated with the species
# param keepByTime an age at which to keep preserve tips before the present (time = 0)
# param tol edge length precision in time cut (default = 1e-08)
# return a list with the prunedPhylogeny and a prunedData
# Mark Puttick

removeNonBin <- function(phy, traitData, keepByTime=0, tol=1e-08){

    all.times <- nodeTimes(phy)
    tip.time <- all.times[which(phy$edge[,2] <= Ntip(phy)) , 2]
	non.bin <- which(tip.time > (keepByTime - tol))
	traitLess <- traitData[ - non.bin]
   	phylo <- drop.tip(phy, non.bin)
	return(list(prunedPhy=phylo, prunedData=traitLess))
}



# Sample hidden speciation events along branches of a tree (internal function)
# Uses estimated speciation and extinction rates to sample the number of speciation events 'hidden' by subsequent extinction on each branch of a tree following Bokma (2008). For use with the \code{psi} and \code{multipsi} models. 
# param phy An object of class \code{phylo} (see \pkg{ape}).
# param lambda.sp Estimate of the rate of speciation "lambda"
# param mu.ext Estimate of the rate of extinction "mu"
# param useMean A logical indicating whether to output the average or expected number of hidden speciation events per branch, which may be non-integer (if \code{TRUE}), or to sample an integer number on each branch from a Poisson distribution (if \code{FALSE}, the default)
# details The expected number of hidden speciation events are calculated for each branch given its start and end times, and estimates of lambda and mu which are assumed to be constant across the tree. To properly account for uncertainty in the effect of extinction on the number of nodes affecting each branch of a tree, it may be appropriate to repeat model-fitting on many realizations of \code{Sobs} on the tree of interest (similar to evaluating phylogenetic uncertainty)
# eturn Phylogenetic tree in \code{phylo} format, with an added element \code{Sobs}, a vector of numbers of hidden speciation events per branch, in the same order as the branches in the \code{phylo} object
# author Travis Ingram

sampleHiddenSp <- function (phy, lambda.sp = NULL, mu.ext = NULL, useMean = FALSE) 
{
    if (is.null(lambda.sp) | is.null(mu.ext)) stop("Please provide values for lambda and mu ")  
    phy$node.label <- NULL
    edge.node.times <- nodeTimes(phy)
    to <- edge.node.times[ , 1]
    te <- edge.node.times[ , 2]
    te[te < 0] <- 0
    if (mu.ext > 0) {
        if (lambda.sp == mu.ext) {
            expSh <- 2 * lambda.sp * (to - te) + 2 * log((1 + lambda.sp * te) / (1 + lambda.sp * to))
        } else {
            expSh <- 2 * lambda.sp * (to - te) + 2 * log((lambda.sp * exp(te *  (lambda.sp - mu.ext)) - mu.ext) / (lambda.sp * exp(to * (lambda.sp - mu.ext)) - mu.ext))
        }
        if (useMean) {
            hidden.speciation <- expSh
        } else {
            hidden.speciation <- rpois(length(expSh), expSh)
        }
	} else {
        hidden.speciation <- rep(0, length(to))
    }
    phy$hidden.speciation <- hidden.speciation
    phy
}

# Slice tree (internal function)
# Split tree for the time slice function 
# phy Phylogenetic tree in phylo format
# splitTime Split time in the past
# Branch lengths of sliced tree
# Mark Puttick

sliceTree <- function(phy, splitTime) {
  branchLengths <- c()
  node.times.2 <- node.times <- nodeTimes(phy)
  startTime <- node.times[1, 1]
  split.these.times <- c(splitTime, startTime)
  treeBrLength <- phy$edge.length
  
  for (p in 1:length(split.these.times)) {
    startTime.2 <- startTime - split.these.times[p]
    br.l <- rep(0, dim(node.times)[1])
    these.br <- which(node.times[, 1] > startTime.2)
    br.l[these.br] <- node.times[these.br, 1] - startTime.2
    treeBrLength <- node.times[, 1] - br.l
    intNode <- which(node.times[these.br, 2] != 0)
    test.finish <-
      which(treeBrLength[these.br][intNode] <= node.times[these.br, 2][intNode])
    node.times[these.br, 1] <- treeBrLength[these.br]
    if (length(test.finish) > 0) {
      br.l[these.br][intNode][test.finish] <-
        node.times.2[these.br, 1][intNode][test.finish] - node.times[these.br, 2][intNode][test.finish]
      node.times[these.br, 1][intNode][test.finish] <- 0
    }
    node.times.2 <- node.times
    branchLengths <- cbind(branchLengths, br.l)
  }
  return(branchLengths)
}


# birth death from APE for non-ultrametric trees

birthdeath_motmot <- function(phy) {
	if (!inherits(phy, "phylo")) 
        stop("object \"phy\" is not of class \"phylo\"")
    N <- length(phy$tip.label)
    x <- c(NA, nodeTimes(phy)[match(unique(phy$edge[,1]),phy$edge[,1]),1])
    dev <- function(a, r) {
        if (r < 0 || a > 1) 
            return(1e+100)
        -2 * (lfactorial(N - 1) + (N - 2) * log(r) + r * sum(x[3:N]) + 
            N * log(1 - a) - 2 * sum(log(exp(r * x[2:N]) - a)))
    }
    out <- nlm(function(p) dev(p[1], p[2]), c(0.1, 0.2), hessian = TRUE)
    if (out$estimate[1] < 0) {
        out <- nlm(function(p) dev(0, p), 0.2, hessian = TRUE)
        para <- c(0, out$estimate)
        inv.hessian <- try(solve(out$hessian))
        se <- if (class(inv.hessian) == "try-error") 
            NA
        else sqrt(diag(inv.hessian))
        se <- c(0, se)
    } else {
        para <- out$estimate
        inv.hessian <- try(solve(out$hessian))
        se <- if (class(inv.hessian) == "try-error") 
            c(NA, NA)
        else sqrt(diag(inv.hessian))
    }
    Dev <- out$minimum
    foo <- function(which, s) {
        i <- 0.1
        if (which == 1) {
            p <- para[1] + s * i
            bar <- function() dev(p, para[2])
        } else {
            p <- para[2] + s * i
            bar <- function() dev(para[1], p)
        }
        while (i > 1e-09) {
            while (bar() < Dev + 3.84) p <- p + s * i
            p <- p - s * i
            i <- i/10
        }
        p
    }
    CI <- mapply(foo, c(1, 2, 1, 2), c(-1, -1, 1, 1))
    dim(CI) <- c(2, 2)
    names(para) <- names(se) <- rownames(CI) <- c("d/b", "b-d")
    colnames(CI) <- c("lo", "up")
    obj <- list(tree = deparse(substitute(phy)), N = N, dev = Dev, 
        para = para, se = se, CI = CI)
    class(obj) <- "birthdeath"
    obj
    }
    
    
    	# my function to get edge lengths from a vcv matrix
	vcv2edge.length <- function(vcv, tre) {
		mrca.phy <- mrca(tre)
		inner <- tre$edge[,1]
		outer <- tre$edge[,2]
		out.len <- c()
		for(i in 1:length(tre$edge[,1])) {
			vcv.in <- inner[i]
			vcv.out <- outer[i]
			if(vcv.in == (Ntip(tre) + 1)) {
				out.len[i] <- vcv[which(mrca.phy == vcv.out)][1]
			} else {
				out.len[i] <- vcv[which(mrca.phy == vcv.out)][1] - vcv[which(mrca.phy == vcv.in)][1]
				}
			}
		return(out.len)
	}
	

# @title Internal function
# @description Internal function. Constructor to allow fixing of rate parameters.
# @param rateData an object of class \code{rateData}
# @param fixed A vector stating whether each parameter should be allowed to vary (either \code{FALSE} which results in a start value of 1, or a numeric start value) or should be fixed (\code{TRUE}).
# @param common.mean a logical specififying whether each rate category should have its own mean (\code{common.mean=FALSE}) or all categories should have the same mean (\code{common.mean=FALSE}). See Thomas et al. (2009) for a discussion on the impact of assumptions about mean on rate estimates..
# @param lambda.est Logical - Logical. Fit Pagel's lambda.
# @param meserr an object of class "rateData"
# @return Returns a function to be passed to \code{optim.likRatePhylo}
# @author Gavin Thomas
# @export

make.likRatePhylo <-
function(rateData, fixed, common.mean=FALSE, lambda.est, meserr) {
	
	op <- c(fixed, !lambda.est, !meserr)
	
	function(params){
		
	
		op[c(!fixed, lambda.est, meserr)] <- params[c(!fixed, lambda.est, meserr)]
		

		y <- rateData$y
		x <- as.factor(rateData$x)
		if (common.mean==FALSE) {k <- nlevels(x)} else {k <- 1}
		
		V <- rateData$Vmat
		v1 <- V[[1]]
		nV <- length(rateData$Vmat)
		
		rateMats <- vector(mode="list", length = nV)
		retMat <- matrix(0, nrow = dim(v1)[1], ncol = dim(v1)[2])
		
		for(i in 1:nV) {
			rateMats[[i]] <- op[i] * V[[i]]  
			retMat <- retMat + rateMats[[i]]
		}
										
			x <- make.anc(y, x)
			
			if(common.mean==FALSE) {
				x <- x} else { x <- rep(1, length(x[,1]))}
									
				V <- retMat
				v.temp <- V
				diag(v.temp) <- rep(0, dim(V)[1])
			
				V.lam <- op[nV+1]*v.temp
				diag(V.lam) <- diag(V)
				V <- V.lam
				if (meserr==TRUE) {	diag(V) <- diag(V) + rateData$meserr/op[nV+2] } 

				logDetV <- determinant(V)$modulus
				
				iV <- solve(V)
				xVix <- crossprod(x, iV %*% x)
				xViy <- crossprod(x, iV %*% y)
				mu <- solve(xVix) %*% xViy 
				
				e <- y - x %*% mu
				s2 <- crossprod(e, iV %*% e)
				n <- length(y) 
				phylo.var <- ( s2 / (n - k) )
				
				n <- length(y)
				ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(phylo.var) - logDetV / 2.0 - (n - k)/2.0
				
				lik.RatePhylo <- ( list(ll = ll, mu = mu, phylo.var = phylo.var, lambda=op[nV+2]) )
				return(-1 * lik.RatePhylo$ll)	
		}
	
	}


# Confidence intervals for rate parameters
# Calculates approximate confidence intervals for all rate parameters. CIs are esimated for one rate parameters while fixing others at a given value (usually the maximum likelihood estimate).
# These are reliable (given the asympotic assumptions of the chi-square distribution) if only two groups are being compared but should be regarded only as a rough approximation for =>3 different rate categories. If the rates are correlated the CIs may be underestimated.
# param rateData an object of class \code{rateDate}
# param MLrate a vector of relative rate parameters. The length of the vector is equal to the number of rates being estimated. If \code{rate=NULL} then rates are equal. Normally these will be the maximum likelihood rate estimates.
# param fixed A vector stating whether each parameter should be allowed to vary (either \code{FALSE} which results in a start value of 1, or a numeric start value) or should be fixed (\code{TRUE}).
# param rateMIN Minimum value for the rate parameters
# param rateMAX Maximum value for the rate parameters
# param common.mean a logical specififying whether each rate category should have its own mean (\code{common.mean=FALSE}) or all categories should have the same mean (\code{common.mean=FALSE}). See Thomas et al. (2009) for a discussion on the impact of assumptions about mean on rate estimates.
# param lambda.est Logical. Estimate Pagel's lambda.
# return rateLci Lower confidence interval for rate estimate
# return rateUci Upper confidence interval for rate estimate
# Gavin Thomas, Rob Freckleton
# data(anolis.tree)
# data(anolis.data)
# 
# ## Convert data to class rateData with a rateMatrix object as input
# anolis.rateMatrix <- as.rateMatrix(phy=anolis.tree, x="geo_ecomorph", data=anolis.data)
# anolis.rateData <- as.rateData(y="Female_SVL", x="geo_ecomorph", 
# rateMatrix = anolis.rateMatrix, phy=NULL, data=anolis.data, log.y=TRUE)
#
# # A model with a different rate in each of the four groups. The 'fixed' command is used to determine
# # determine whether a particular rate is to be constrained or not. Use '1' to fix a group and
# # 'FALSE' to show that the parameter is not fixed and should be estimated. The values
# #  should be entered in the same order as the ranking of the groups. That is, group
# # 0 (small islands) takes position one in the fixed vector, group 1 (large island trunk crown 
# # and trunk ground) takes position 2 and so on. The default is to allow each group 
# # to take a different mean. 
#
#anole.ML <- optim.likRatePhylo(rateData=anolis.rateData, rate=NULL, 
#fixed=c(FALSE,FALSE,FALSE, FALSE), common.mean=FALSE, lambda.est=FALSE) # ML estimates
#
# # Confidence intervals for the first parameters
#
# RatePhylo.CI(rateData=anolis.rateData, MLrate = anole.ML$MLRate, 
# fixed=c(FALSE, TRUE, TRUE, TRUE), rateMIN = 0.001, rateMAX = 50, 
# common.mean = FALSE)

RatePhylo.CI <-
function(rateData, MLrate=NULL, fixed, rateMIN = 0.001, rateMAX = 50, common.mean=FALSE, lambda.est=TRUE) {

	if(is.null(MLrate))  { MLrate <- c(rep(1,length(rateData$Vmat))) } else { MLrate <- MLrate }	
	
	MLrate <- as.numeric(format(MLrate))

	ML <- likRatePhylo(rateData, MLrate, common.mean=common.mean, lambda.est=lambda.est)$ll
	
		fixed.rates <- data.frame(MLrate, fixed)
		use.rate.ind <- which(fixed.rates$fixed==FALSE)
		use.rate <- fixed.rates$MLrate[use.rate.ind]
	
		var.fun <- function(vary.rate) {
		
				fixed.rates$MLrate[use.rate.ind] <- vary.rate
				test.rate <- fixed.rates$MLrate
				ll <- likRatePhylo(rateData, test.rate, common.mean=common.mean, lambda.est=lambda.est)$ll
		
		return( ll - ML + 1.92)
		 }
		 
	if(var.fun(rateMIN) < 0) { 
			rateLci <- uniroot(var.fun, interval = c(rateMIN, use.rate))$root 
			}
	if(var.fun(rateMAX) < 0) {
			rateUci <- uniroot(var.fun, interval = c(use.rate, rateMAX))$root
			}
			return(c(rateLci=rateLci, rateUci=rateUci))
			}


# #' @title timeTravelPhy (internal function) 
# #'
# #' @description removes tips and lineages after a time in the past
# #' @param phy An object of class "phylo" (see ape package).
# #' @param node nodes arising more recently than the cut time
# #' @param nodeEstimate trait the number of descendants arising from the nodes
# #' @param timeCut position at which to cut the phylogeny
# #' @param traits Logical. Include trait values in the output tree
# #' @return the pruned phylogeny and a 'tipObject' of the number of lineages found in the pruned branches
# #' @author Mark Puttick

timeTravelPhy <- function(phy, node, nodeEstimate, timeCut, traits=TRUE){

	countDesAll <- c()
	node.in <- match(unique(phy$edge[,1]), phy$edge[, 1])
	originalNodeTimes <- nodeTimes(phy)[node.in, 1]
  	names(originalNodeTimes) <- unique(phy$edge[ ,1])
	internal <- phy$edge[, 1]
	external <- phy$edge[, 2]

	if(length(node) > 0) {

		for (i in 1:length(node)) {
			branchOne <- which(internal == node[i])
			countDes <- ext <- external[branchOne]
			isTerminal <- match(internal, ext)
			terminalCheck <- which(is.na(isTerminal) == FALSE)
			loopStop <- any(terminalCheck)
			
			while (loopStop == TRUE) {
					ext <- external[terminalCheck]
					countDes <- c(countDes, ext)
					isTerminal <- match(internal, ext)
					terminalCheck <- which(is.na(isTerminal) == FALSE)
					loopStop <- any(terminalCheck)
				}

			alreadyCounted <- which(is.na(match(countDes, countDesAll)) == FALSE)
			if (sum(alreadyCounted) > 0) countDes <- countDes[-alreadyCounted]
			countDesAll <- c(countDes, countDesAll)
			}

	tre2 <- phy
	tipsInSample <- countDesAll[which(countDesAll <= Ntip(phy))]
	desInMat <- match(countDesAll, phy$edge[,2])
	tre2$edge <- phy$edge[-desInMat, ]
	keepN <- match(tre2$edge[, 2], node)
	nodeToTips <- tre2$edge[complete.cases(keepN), 2]
	nodalPosition <- which(is.na(keepN) != TRUE)
	originalNodes <- unique(tre2$edge[, 1])
	for (k in 1:length(node)) tre2$edge[which(tre2$edge[ ,2] == node[k]), 2] <- NA
		
	naMat <- which(is.na(tre2$edge[, 2]))
	tipsRemain <- which(tre2$edge[, 2] <= Ntip(phy))
	keeperTips <- tre2$edge[tipsRemain, 2]

	
	tre2$edge[which(tre2$edge[, 2] <= Ntip(phy)), 2] <- NA
	tre2$edge.length <- phy$edge.length[-desInMat]
	internals <- tre2$edge[ ,1]
	tipsInTre2 <- length(which(is.na(tre2$edge[ , 2])))
	n_node <- tipsInTre2 + 1
	nodeTodal <- seq(tre2$edge[1,1], tre2$edge[1,1] + dim(tre2$edge)[1] / 2 - 1, by=1)
	tre2NodeTotal <- match(tre2$edge[,1], nodeTodal)
	nodeTodalNode <- match(nodeTodal, tre2$edge[,1])
	notFoundNode <- nodeTodal[which(is.na(nodeTodalNode))]
	uniqueTre2Int <- unique(tre2$edge[which(is.na(tre2NodeTotal)), 1])
	for(y in 1:length(uniqueTre2Int)) tre2$edge[which(tre2$edge[,1] == uniqueTre2Int[y]), 1] <- notFoundNode[y]
	
	takeValue <- tre2$edge[1, 1] - n_node
	tre2$edge[ ,1] <- tre2$edge[,1] - takeValue
	addOnNodes <- dim(tre2$edge)[1] / 2 - 1
	internalNewTree <- rep(seq(n_node, n_node + addOnNodes, 1), each=2)
	externalNew <- rep(unique(tre2$edge[ ,1]), each=2)
	internalNewMatch <- match(tre2$edge[ ,1], externalNew)
	tre2$edge[,1] <- internalNewTree[internalNewMatch]
	buildInt <- which(complete.cases(tre2$edge[ ,2]) == T)
	extNode <- n_node + 1
	desNodes <- seq(from=extNode, to=extNode + length(buildInt) - 1, 1)
	tre2$edge[buildInt, 2] <- desNodes
	tre2$edge[which(is.na(tre2$edge[,2])), 2] <- 1:tipsInTre2
	apeTree <- list(edge=tre2$edge, tip.label=c(1:tipsInTre2), Nnode=tipsInTre2-1, edge.length=tre2$edge.length)
	class(apeTree) <- "phylo"
	
	toKeep <- apeTree$edge[naMat, 2]
	tipsToKeep <- match(c(1:tipsInTre2), toKeep)
	
	new.tip <- 1:Ntip(apeTree)
	new.tip[-apeTree$edge[naMat, 2]] <- phy$tip.label[-tipsInSample]
	new.tip[apeTree$edge[naMat, 2]] <- paste0("node_", unique(phy$edge[match(nodeToTips, phy$edge[,2]), 2]))
	apeTree$tip.label <- new.tip
	
	if(traits == TRUE) {
		
		keepTipPhy <- match(keeperTips, phy$edge[,2])
		keepTipPhyValue <- nodeEstimate[keepTipPhy]
		newTips <- match(nodeToTips, phy$edge[,2])
		tipObject <- rep(NA, Ntip(apeTree))
		tipObject[which(is.na(tipsToKeep))] <- as.numeric(keepTipPhyValue)	
		needNewTips <- which(is.na(tipObject))
		ageOld <- match(nodeToTips, as.numeric(names(originalNodeTimes)))

		for(p in 1:length(needNewTips)){
			origNode_n <- originalNodeTimes[ageOld][p]		
			nodeN <- nodeToTips[p]
			pastNode <- phy$edge[which(phy$edge[,2] == nodeN), 1]
			pastLength <- phy$edge.length[which(phy$edge[,2] == nodeN)]
			ancAge <- origNode_n + pastLength			
			decAge <- origNode_n
			findNear <- which(c(ancAge - timeCut, - (decAge - timeCut)) == min(c(ancAge - timeCut, - (decAge - timeCut))))
			nearestNode <- which(phy$edge[,2] == c(pastNode, nodeN)[findNear])
			tipObject[needNewTips[p]] <- as.numeric(nodeEstimate[nearestNode])[1]
			}
		
		names(tipObject) <- NULL
		}
	
	} else {
	
	apeTree <- phy	
	tipObject <- nodeEstimate[which(phy$edge[,2] <= Ntip(phy))]
	
	}
	
	node.in.ape <- match(unique(apeTree$edge[,1]), apeTree$edge[,1])
	all.ages.ape <- nodeTimes(apeTree)
	nodeTimesApe <- all.ages.ape[node.in.ape,1]
  	names(nodeTimesApe) <- unique(apeTree$edge[,1])
	tipTimesApe <- all.ages.ape[which(apeTree$edge[,2] <= Ntip(apeTree)), 2] 
	names(tipTimesApe) <- apeTree$edge[which(apeTree$edge[,2] <= Ntip(apeTree)), 2]
	cutOffAge <- nodeTimes(phy)[1,1] - timeCut
	cutOff <- nodeTimes(apeTree)[1,1] - cutOffAge
	preCutOff <- which(tipTimesApe < cutOff)
	preCutOffTips <- match(as.numeric(names(preCutOff)), apeTree$edge[ ,2])
	tipsPreCut <- tipTimesApe[preCutOff] - cutOff
	tipsInApe <- match(preCutOff, apeTree$edge[,2])
	
	apeTree$edge.length[preCutOffTips] <- apeTree$edge.length[preCutOffTips] + tipsPreCut

	if(traits == TRUE) return(list(phy=apeTree, tipData=tipObject))
	if(traits == FALSE) return(phy=apeTree)
}


reml.lik <- function(node.state, bm.rate, vcv.phylo, y) {
	n <- ncol(vcv.phylo)
	s.b <- (t(y - rep(node.state, n)) %*% (solve(vcv.phylo)) %*% (y - rep(node.state, n)) / bm.rate)[1,1]
	ll.in <- -0.5 * ((n) * log(2 * pi * bm.rate) + log(det(vcv.phylo)) + s.b)
	return(ll.in )
}
      
mu.mean <- function(vcv.phy, y.in) {
	x <- rep(1, ncol(vcv.phy))
	iV <- solve(vcv.phy)
	xVix <- crossprod(x, iV %*% x)
	xViy <- crossprod(x, iV %*% y.in)
	mu <- solve(xVix) %*% xViy
	return(mu)
	}

sig.sq <- function(mu, vcv.phy, y) {
	n <- ncol(vcv.phy)
	return((1 / (n - 1)) * (t(y - rep(mu, n)) %*% (solve(vcv.phy)) %*% (y - rep(mu, n))))
	}
		



