#' @title Phylogenetic tree transformations
#' @description Transforms the branch lengths of a phylo object according to a model of trait evolution (see details).
#' @param y A matrix of trait values.
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param model The model of trait evolution (see details).
#' @param kappa Value of kappa transform.
#' @param lambda Value of lambda transform.
#' @param delta Value of delta transform.
#' @param alpha Value of alpha (OU) transform.
#' @param psi Value of psi transform.
#' @param acdcRate Value of ACDC transform.
#' @param lambda.sp Estimate of speciation (lambda) for the psi models
#' @param nodeIDs Integer - ancestral nodes of clades.
#' @param splitTime A split time (measured from the present, or most recent species) at which a shift in the rate occurs for the "timeSlice" model
#' @param timeRates The rates (from ancient to recent) for the timeSlice model
#' @param rateType If model="clade", a vector specifying if rate shift occurs in a clade ("clade") or on the single branch leading to a clade ("branch").
#' @param branchRates Numeric vector specifying relative rates for individual branches
#' @param cladeRates Numeric vector specifying telative rates for clades or logical to indicate scalar is included in the 'modeslice' model (the scalar is included in the mode.param argument with the 'modeslice' model).
#' @param branchLabels A vector of length equal to the number of internal branches, signifying the which "multiPsi" class it belongs to
#' @param meserr A vector (or matrix) of measurement error for each tip. This is only applicable to univariate analyses. Largely untested - please use cautiously
#' @param cophenetic.dist a cophenetic distance matrix showing the absolute distance between taxa - only applicable for OU models run on non-ultrmetric trees. If null will be calculated internally, but supplying the data can speed up run time
#' @param vcv.matrix a variance-covariance matrix - only applicable for OU models run on non-ultrmetric trees. If null will be calculated internally, but supplying the data can speed up run time
#' @param mode.order The order of modes for the 'modeslice' model. Any combination of 'BM', 'OU', 'acdc', and 'kappa'
#' @param rate.var Allows rate variation in BM modes in the 'modeslice' model
#' @param mode.param Parameters for the modes of evoluton in the 'modeslice' model
#' @details Transforms the branch lengths of a phylo object according to one of the following models:
#' \itemize{
#'  \item {model="bm"}- Brownian motion (constant rates random walk)
#'  \item {model="kappa"} - fits Pagel's kappa by raising all branch lengths to the power kappa. As kappa approaches zero, trait change becomes focused at branching events. For complete phylogenies, if kappa approaches zero this infers speciational trait change. 
#'  \item {model="lambda"} - fits Pagel's lambda to estimate phylogenetic signal by multiplying all internal branches of the tree by lambda, leaving tip branches as their original length (root to tip distances are unchanged);
#'  \item {model="delta"} - fits Pagel's delta by raising all node depths to the power delta. If delta <1, trait evolution is concentrated early in the tree whereas if delta >1 trait evolution is concentrated towards the tips. Values of delta above one can be difficult to fit reliably.
#'  \item {model="free"} - fits Mooer's et al's (1999) free model where each branch has its own rate of trait evolution. This can be a useful exploratory analysis but it is slow due to the number of parameters, particularly for large trees.
#'  \item {model="clade"} - fits a model where particular clades are a priori hypothesised to have different rates of trait evolution (see O'Meara et al. 2006; Thomas et al. 2006, 2009). Clades are specified using nodeIDs and are defined as the mrca node. Unique rates for each clade are specified using cladeRates. rateType specifies whether the rate shift occurs in the stem clade or on the single branch leading to the clade.
#'  \item {model="OU"} - fits an Ornstein-Uhlenbeck model - a random walk with a central tendency proportional to alpha. High values of alpha can be interpreted as evidence of evolutionary constraints, stabilising selection or weak phylogenetic signal.
#'  \item {model="psi"} - fits a model to assess to the relative contributions of speciation and gradual evolution to a trait's evolutionary rate (Ingram 2010). Note that 'original nodes' from the full phylogeny can be included as an element on the phylogeny (e.g., phy$orig.node) as well as estimates of 'hidden' speciation (e.g., phy$hidden.speciation) if estimates of extinction (mu) are > 0.
#' \item {model="multiPsi"} {fits a model to assess to the relative contributions of speciation and gradual evolution to a trait's evolutionary rate but allows seperate values of psi fitted to seperate branches (Ingram 2010; Ingram et al. 2016). Note that 'original nodes' from the full phylogeny can be included as an element on the phylogeny (e.g., phy$orig.node) as well as estimates of 'hidden' speciation (e.g., phy$hidden.speciation) if estimates of extinction (mu) are > 0.}
#' \item {model="ACDC"} {fits a model to in which rates can exponentially increased or decrease through time (Blomberg et al. 2003). If the upper bound is < 0, the model is equivalent to the 'Early Burst' model of Harmon et al. 2010. If a nodeIDs is supplied, the model will fit a ACDC model nested within a clade, with a BM fit to the rest of the tree.}
#' \item {model="timeSlice"} {A model in which all branch rates change at time(s) in the past.}
#' \item {model="modeSlice"} {A model in which all branch modes change at a time or times set a priori by the user.}
#'  }
#' @return phy  A phylo object with branch lengths scaled according to the given model and parameters
#' @seealso \code{\link{transformPhylo.ML}}, \code{\link{transformPhylo.ll}}, \code{\link{transformPhylo.MCMC}}
#' @references Ingram T. 2011. Speciation along a depth gradient in a marine adaptive radiation. Proc. Roy. Soc. B. 278, 613-618.
#' @references Ingram T,  Harrison AD, Mahler L, Castaneda MdR, Glor RE, Herrel A, Stuart YE, and Losos JB. Comparative tests of the role of dewlap size in Anolis lizard speciation. Proc. Roy. Soc. B. 283, 20162199.
#' @references Mooers AO, Vamosi S, & Schluter D. 1999. Using phylogenies to test macroevolutionary models of trait evolution: sexual selection and speciation in Cranes (Gruinae). American Naturalist 154, 249-259.
#' @references O'Meara BC, Ane C, Sanderson MJ & Wainwright PC. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution 60, 922-933
#' @references Pagel M. 1997. Inferring evolutionary processes from phylogenies. Zoologica Scripta 26, 331-348.
#' @references Pagel M. 1999 Inferring the historical patterns of biological evolution. Nature 401, 877-884.
#' @references Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.
#' @author Gavin Thomas, Mark Puttick
#' @examples
# Apply delta transform to anolis tree
#' data(anolis.tree)
#' anolis.treeDelta <- transformPhylo(phy=anolis.tree, model="delta", delta=0.5)
#' @export
#' @import caper

transformPhylo <- function (phy, model = NULL, y = NULL, meserr=NULL, kappa = NULL, lambda = NULL, delta = NULL, alpha = NULL, psi = NULL, lambda.sp = NULL, nodeIDs = NULL, rateType = NULL, branchRates = NULL, cladeRates = NULL,  splitTime = NULL, timeRates = NULL, acdcRate=NULL,  branchLabels = NULL, cophenetic.dist=NULL, vcv.matrix=NULL, mode.order=NULL, mode.param=NULL, rate.var=NULL) {

    n <- length(phy$tip.label)
	
   if (is.null(meserr) == FALSE) {
        if (dim(y)[2] > 1) {
            meserr <- NULL
            (stop("Measurement error can only be included for univariate models. Set meserr to NULL."))
        }
    }
    
    model <- tolower(model)
		
    switch(model, bm = {
		   if (is.null(meserr) == FALSE) {
		   height <- nodeTimes(phy)[1,1]
		   interns <- which(phy$edge[, 2] > n)
		   externs <- which(phy$edge[, 2] <= n)
		   phy$edge.length[externs] <- phy$edge.length[externs] + (meserr ^ 2) / (var(y) / height)[1]
		   } else {
		   phy <- phy
		   }
		   
		   }, kappa = {
		   	
		   	if(is.null(nodeIDs)) {
		   		node <- Ntip(phy) + 1	
		   		} else {
		   		node <- nodeIDs
		   		}
		   	
		   	if (is.null(meserr) == FALSE) {
		   		height <- nodeTimes(phy)[1,1]
		   		interns <- which(phy$edge[, 2] > n)
		   		externs <- which(phy$edge[, 2] <= n)
		   		}
		   		
		   	if(node == (Ntip(phy) + 1)) {
		   		relations_num <- 1:(Nnode(phy) + Ntip(phy))
		   	} else {
		   		relations_num <- c(node.descendents(node, phy))
		   	}
		   	
		   	branches <- match(relations_num, phy$edge[,2])
		   	if(any(is.na(branches))) branches <- branches[complete.cases(branches)]
		   	branches <- sort(branches)
		   	phy$edge.length[branches] <- phy$edge.length[branches] ^ kappa
		   	if (is.null(meserr) == FALSE) {
		   		phy$edge.length[externs] <- phy$edge.length[externs] + (meserr^2)/(var(y)/height)[1]
        		}
		   
		   }, lambda = {
		   if (is.null(meserr) == FALSE) {
		   height <- nodeTimes(phy)[1,1]
		   interns <- which(phy$edge[, 2] > n)
		   externs <- which(phy$edge[, 2] <= n)
		   }
		   if (is.ultrametric(phy)) {
		   rootOrig <- nodeTimes(phy)[1,1]
		   tips <- match(c(1:Ntip(phy)), phy$edge[, 2])
		   phy$edge.length <- phy$edge.length * lambda
		   phy$edge.length[tips] <- phy$edge.length[tips] + (rootOrig * (1 - lambda))
		   }
		   if (is.ultrametric(phy) == FALSE) {
		   tips <- match(c(1:Ntip(phy)), phy$edge[, 2])
		   cladeMat <- caper::clade.matrix(phy)
		   branchHeights <- sapply(1:Ntip(phy), function(x) sum(cladeMat$edge.length[cladeMat$clade.matrix[, x] == 1]))
		   phy$edge.length <- phy$edge.length * lambda
		   phy$edge.length[tips] <- phy$edge.length[tips] + (branchHeights * (1 - lambda))
		   }
		   if (is.null(meserr) == FALSE)	{
		   	phy$edge.length[externs] <- phy$edge.length[externs] +  (meserr^2)/(var(y)/height)[1]
		   	}
		   
		   
		   }, delta = {
		   
		   if(is.null(nodeIDs)) {
		   	node <- Ntip(phy) + 1
		   	} else { 
		   	node <- nodeIDs
		   	}
		   allTimes <- nodeTimes(phy)
		   if(node == phy$edge[1, 1]) {
		   	 height <- allTimes[1,1]
		   } else {
		   	height <- allTimes[which(node == phy$edge[,2]), 2]
		   }
		   	
	   if (is.null(meserr) == FALSE) {
		   	interns <- which(phy$edge[, 2] > n)
		   	externs <- which(phy$edge[, 2] <= n)
		   }
		   
		   times <- allTimes[match(c((Ntip(phy) + 1) : (Nnode(phy) + Ntip(phy))), phy$edge[,1]) , 1]
		   names(times) <- (Ntip(phy) + 1):(Ntip(phy) + Nnode(phy))
		   
		   if(node == (Ntip(phy) + 1)) {
		   		relations_num <- 1:(Nnode(phy) + Ntip(phy))
		   	} else {
		   		relations_num <- c(node.descendents(node, phy))
			}
		   	
		   originTime <- times[which(names(times) == node)][1]
		   branches <- match(relations_num, phy$edge[,2])
		   if(any(is.na(branches))) branches <- branches[complete.cases(branches)]
			
		   times <- originTime - times
		   tips <- Ntip(phy)
		   res <- phy
		   branches <- sort(branches)
		   		   
		   res$edge.length[branches] <- sapply(branches, function(x) {
				bl <- phy$edge.length[x]
		   		age <- times[phy$edge[x, 1] - tips]
		   		(age + bl)^delta - age^delta
		   		}
		   	)
		   	
		   phy <- res
		   if (is.null(meserr) == FALSE)	{
		   	phy$edge.length[externs] <- phy$edge.length[externs] + (meserr^2)/(var(y)/height)[1]
		   	}
		   	
		   # allTimes2 <- nodeTimes(phy)
		   # if(node == phy$edge[1, 1]) {
		   	# height2 <- allTimes2[1,1]
		   # } else {
		   	# height2 <- allTimes2[which(node == phy$edge[,2]), 2]
		   # }
		   # height2 <- allTimes[1,1]
		   # phy$edge.length[branches] <- phy$edge.length[branches] * (height/ height2)
		   		   
		   }, free = {
		   if (is.null(meserr) == FALSE) {
		   height <- nodeTimes(phy)[1,1]
		   interns <- which(phy$edge[, 2] > n)
		   externs <- which(phy$edge[, 2] <= n)
		   }
		   branchRates <- branchRates + (1 - min(branchRates))
		   phy$edge.length <- phy$edge.length * branchRates
		   if (is.null(meserr) == FALSE) {
		   	phy$edge.length[externs] <- phy$edge.length[externs] + (meserr^2)/(var(y)/height)[1]
		   	}
		      
		   }, clade = {
		   if (is.null(meserr) == FALSE) {
		   height <- nodeTimes(phy)[1,1]
		   interns <- which(phy$edge[, 2] > n)
		   externs <- which(phy$edge[, 2] <= n)
		   }
		   if (is.null(rateType)) {
		   rateType <- rep("clade", length(nodeIDs))
		   } else {
		   rateType <- rateType
		   }
		   branchShiftNms <- branchShiftID <- cladeShiftNms <- cladeShiftID <- NULL
		   cladeMembers <- matrix(0, ncol = length(nodeIDs), nrow = length(phy$edge[, 1]))
		   shiftType <- data.frame(rateType, nodeIDs, cladeRates)
		   colnms <- paste(shiftType[, 1], shiftType[, 2], sep = "")
		   shiftType <- data.frame(shiftType, colnms)
		   if (sum(shiftType[, 1] == "clade") > 0) {
		   cladeShiftID <- shiftType[shiftType[, 1] == "clade", 2]
		   cladeShiftNms <- as.character(shiftType[shiftType[, 1] == "clade", 4])
		   cladeMembers[, 1:length(cladeShiftID)] <- cladeIdentity(phy = phy, nodeIDs = cladeShiftID) #cladeMembersObj=cladeMembersObj)
		   }
		   if (sum(shiftType[, 1] == "branch") > 0) {
		   branchShiftID <- shiftType[shiftType[, 1] == "branch", 2]
		   branchShiftNms <- as.character(shiftType[shiftType[, 1] == "branch", 4])
		   }
		   if (is.null(branchShiftNms) == FALSE) { 
		   	colnames(cladeMembers) <- c(cladeShiftNms, rep(NA, length(branchShiftNms))) 
		   	} else  { 
		   	colnames(cladeMembers) <- cladeShiftNms
		   	}
		   	
		   if (sum(shiftType[, 1] == "branch") > 0) {
		   for (i in 1:length(branchShiftID)) {
			   branchID <- which(phy$edge[, 2] == branchShiftID[i])
			   cladeMembers[branchID, ] <- 0
			   cladeMembers[branchID, length(cladeShiftID) + i] <- 1
			   colnames(cladeMembers)[length(cladeShiftID) + i] <- branchShiftNms[i]
			   }
		   }
		   cladeMembers <- as.matrix(cladeMembers[, match(shiftType[, 4], colnames(cladeMembers))])
		   for (i in 1:length(cladeRates)) {
		   	phy$edge.length[cladeMembers[, i] == 1] <- phy$edge.length[cladeMembers[, i] == 1] * cladeRates[i]
		   	}
		   if (is.null(meserr) == FALSE) {
            phy$edge.length[externs] <- phy$edge.length[externs] +(meserr^2)/(var(y)/height)[1]
        	}
		   
		   }, ou = {
		   	
		   is.contemp <- is.ultrametric(phy)
		   
		   if(is.null(nodeIDs)) { 
		   	node <- Ntip(phy) + 1
		   	} else { 
		   	node <- nodeIDs
		   	}
		   	
		   if (is.null(meserr) == FALSE) {
		   		interns <- which(phy$edge[, 2] > n)
		   		externs <- which(phy$edge[, 2] <= n)
		   }
		   
		   times <- nodeTimes(phy)
   		   height <- times[1,1]
		   times <- times[match(c((Ntip(phy) + 1) : (Nnode(phy) + Ntip(phy))), phy$edge[,1]) , 1]
		   names(times) <- (Ntip(phy) + 1):(Ntip(phy) + Nnode(phy))
		   
			if(node == (Ntip(phy) + 1)) {
		   		relations_num <- 1:(Nnode(phy) + Ntip(phy))
		   	} else {
		   		relations_num <- c(node.descendents(node, phy))
		   	}

		   originTime <- times[which(names(times) == node)][1]
		   branches <- match(relations_num, phy$edge[,2])
		   branches <- sort(branches)
		   if(any(is.na(branches))) {
		   	branches <- branches[complete.cases(branches)]
		   	}
			times2 <- nodeTimes(phy)[,1]
		
			phy2 <- phy
			
			if(is.contemp) {
				edge.length.branches <- sapply(branches, function(i) {
					Tmax <- originTime
					bl <- phy$edge.length[i]
					age <- times2[i]
					t1 = originTime - age
					t2 = t1 + bl
					(1/(2 * alpha)) * exp(-2 * alpha * (Tmax - t2)) * (1 - exp(-2 * alpha * t2)) - (1/(2 * alpha)) * exp(-2 * alpha * (Tmax - t1)) * (1 - exp(-2 * alpha * t1))
					}
				)
				phy2$edge.length[branches] <- edge.length.branches
				phy <- phy2
				if (is.null(meserr) == FALSE) {
					phy$edge.length[externs] <- phy$edge.length[externs] + (meserr^2)/(var(y)/height)[1]
		   		}
		   		
			} else {
   	
   				if(is.null(cophenetic.dist)) {
   					t_ij <- cophenetic(phy)
	   				} else {
    				t_ij <- cophenetic.dist
    				}
    			if(is.null(vcv.matrix)) {
    				c_ij <- vcv(phy)
    				} else {
    				c_ij <- vcv.matrix
    				}
    			vcv.out <- (1 / (2 * alpha)) * exp(-alpha * t_ij) * (1 - exp(-2 * alpha * c_ij))
    			if (is.null(meserr) == FALSE) { 
    				diag(vcv.out) <- diag(vcv.out) + (meserr^2)/(var(y)/height)[1]
    			}
    			phy <- vcv.out
			}
		   
		 	
		   		   
		}, psi = {
        if (is.null(meserr) == FALSE) {
            height <- nodeTimes(phy)[1,1]
            interns <- which(phy$edge[, 2] > n)
            externs <- which(phy$edge[, 2] <= n)
        }
		if (is.null(phy$orig.node)) {
			orig.node <- rep(1, length(phy$edge.length))
		} else {
			orig.node <- phy$orig.node
		}
		if (!is.null(phy$hidden.speciation)) orig.node <- orig.node + phy$hidden.speciation
        phy2 <- phy
        phy2$edge.length <- (psi / ( 2* lambda.sp)) * orig.node + (1 - psi) * phy$edge.length
        phy <- phy2
        if (is.null(meserr) == FALSE) {
            phy$edge.length[externs] <- phy$edge.length[externs] + (meserr^2) / (var(y) / height)[1]
        }
        
        
         }, multipsi = {
        if (is.null(meserr) == FALSE) {
            height <- nodeTimes(phy)[1,1]
            interns <- which(phy$edge[, 2] > n)
            externs <- which(phy$edge[, 2] <= n)
        }
		if(is.null(phy$orig.node)){
			orig.node <- rep(1, length(phy$edge.length))
		}else{
			orig.node <- phy$orig.node
		}
		if(!is.null(phy$hidden.speciation)) orig.node <- orig.node + phy$hidden.speciation
        phy2 <- phy
		states <- levels(factor(branchLabels))
		for(i in 1:length(states)) {
            phy2$edge.length[branchLabels==states[i]] <- (psi[i]/(2*lambda.sp)) * orig.node[branchLabels==states[i]] + (1 - psi[i]) * phy$edge.length[branchLabels==states[i]]
		}
        phy <- phy2
        if (is.null(meserr) == FALSE) {
            phy$edge.length[externs] <- phy$edge.length[externs] + (meserr^2)/(var(y)/height)[1]
        	}
        }, 
        
        
        timeslice = {
		   if (is.null(meserr) == FALSE) {
		   height <- nodeTimes(phy)[1,1]
		   interns <- which(phy$edge[, 2] > n)
		   externs <- which(phy$edge[, 2] <= n)
		   }

		   splitTime <- nodeTimes(phy)[1,1] - sort(splitTime, TRUE)
		   sliceLengths <- sliceTree(phy, splitTime)
		   phy$edge.length <- as.vector(timeRates %*% t(sliceLengths))
		   if (is.null(meserr) == FALSE) {
		   	phy$edge.length[externs] <- phy$edge.length[externs] + (meserr ^2 )/(var(y) / height)[1]
		   	}
		   } ,
		   
		   
		modeslice = {
		  		   
		  		   
		  phy.in <- phy
		   if(!is.null(splitTime)) {
				inputSplitTime <- splitTime
				splitTime <- nodeTimes(phy)[1,1] - sort(splitTime, TRUE)
				sliceLengths <- sliceTree(phy, splitTime)
				slice.phy <- lapply(1:length(mode.order), function(sx) {
					phy2 <- phy
					phy2$edge.length <-  sliceLengths[,sx]
					phy2
					}
				)
				all.splits <- c(max(nodeTimes(phy)), inputSplitTime)	
		 	} else {
		 		slice.phy <- list()
		 		slice.phy[[1]] <- phy
		 		all.splits <- max(nodeTimes(phy))	
		 		
		 	}
		 
		 
		all.vcvs <- list()
		count <- 1 	
		for(x in 1:length(mode.order)) {
			
			#bm
			if(mode.order[x] == "bm") {
				if(rate.var) {
					phy2 <- slice.phy[[x]]
					all.vcvs[[x]] <- vcv(phy2) * mode.param[count]
					count <- count + 1
				} else {
					phy2 <- slice.phy[[x]]
					all.vcvs[[x]] <- vcv(phy2)
				}
				
			}
			
			# OU 
			if(mode.order[x] == "ou") {
				alpha <- mode.param[count]
				phy2 <- slice.phy[[x]]
				t_ij <- cophenetic(phy2)
				c_ij <- vcv(phy2)
				all.vcvs[[x]] <- (1 / (2 * alpha)) * exp(-alpha * t_ij) * (1 - exp(-2 * alpha * c_ij))							
				count <- count + 1
			}
			
			

			
			if(mode.order[x] == "kappa") {
				kappa <- mode.param[count]
				phy2 <- slice.phy[[x]]
				phy2$edge.length <- phy2$edge.length ^ kappa
				all.vcvs[[x]] <- vcv(phy2) 
				count <- count + 1
			}
			
			if(mode.order[x] == "acdc") {
				phy2 <- slice.phy[[x]]
				if(isTRUE(cladeRates)) {
					phy2$edge.length <- phy2$edge.length * mode.param[count]
					count <- count + 1
					}
				acdcRate <- mode.param[count]
				originTime <-  all.splits[x]
				if(is.null(originTime)) originTime <- nodeTimes(phy)[1,1]
				edges <- phy2$edge.length
				if(!is.null(splitTime)) {
					times <- sliceLengths[,x]
				} else {
					times <- nodeTimes(phy)[,1]
					}
				
				edgeACDC <- sapply(1:length(edges), function(i) {
					bl <- edges[i]
					age <- times[i]
					t1 <- originTime - age
					t2 <- t1 + bl
					(exp(acdcRate * t2) - exp(acdcRate * t1)) / (acdcRate)
					}
				)
				phy2$edge.length <- edgeACDC
				all.vcvs[[x]] <- vcv(phy2)
				count <- count + 1
				}
			}


		
			
		phy <- all.vcvs[[1]]
		if(length(all.vcvs) > 1) for(u in 2:length(all.vcvs)) phy <- phy + all.vcvs[[u]]
		
		if (is.null(meserr) == FALSE) {
			height <- nodeTimes(phy.in)[1,1]
		   	diag(phy) <- diag(phy) + (meserr ^2 )/(var(y) / height)[1]
		   	}
		   } ,
		   
		   
		acdc = {
		   	
			if(is.null(nodeIDs)) node <- Ntip(phy) + 1 else node <- nodeIDs
			times <- nodeTimes(phy)[,1]
			if (is.null(meserr) == FALSE) {
				height <- times[1]
				interns <- which(phy$edge[, 2] > n)
				externs <- which(phy$edge[, 2] <= n)
				}
			
			names(times) <- phy$edge[,1]

			if(node == (Ntip(phy) + 1)) {
				relations_num <- 1:(Nnode(phy) + Ntip(phy))
				} else {
				relations_num <- c(node, node.descendents(node, phy))
				}

			originTime <- times[which(names(times) == node)][1]
			branches <- match(relations_num, phy$edge[,2])
			if(any(is.na(branches))) branches <- branches[complete.cases(branches)]
			branches <- sort(branches)
			times <- nodeTimes(phy)[,1]
			
			if(is.null(cladeRates) == FALSE) phy$edge.length[branches] <- phy$edge.length[branches] * cladeRates
			edgeACDC <- sapply(branches, function(i) {
				bl <- phy$edge.length[i]
				age <- times[i]
				t1 <- originTime - age
				t2 <- t1 + bl
				(exp(acdcRate * t2) - exp(acdcRate * t1))/(acdcRate)
				})
				
			phy$edge.length[branches] <- edgeACDC	
			if (is.null(meserr) == FALSE)	 {
				phy$edge.length[externs] <- phy$edge.length[externs] + (meserr^2)/(var(y)/height)[1]
					}
				}
		)
		return(phy)
}
