# -----------------------------------------------------------------------------
# +++ recursive merging of clusters based on signal-to-noise-ratio (S2NR)
# -----------------------------------------------------------------------------

s2nr.merge <- function(bdm, k = 10, verbose = T, layer = 1)
{
	# input data
	X <- Xdata.get(bdm$data, bdm$Xdata$whiten, bdm$Xdata$input.dim, bdm$is.distance)$X
	# output data
	l <- c(1, 2) + (layer -1) *2
	Y <- bdm$ptsne$Y[, l]
	# density estimation
	Z <- bdm$pakde[[layer]]$z
	# grid parameters
	G <- bdm$wtt[[layer]]$grid
	# density landscape peaks
	P <- bdm$wtt[[layer]]$P
	# base cell cluster assignment
	C <- bdm$wtt[[layer]]$C
	L <- merge.labels(Y, C, G)

	# base clustering information
	if (verbose) cat('+++ computing clustering hierarchy ...\n')
	K <- s2nr.kList(C, L, X, Z, G, P)

	# base S2NR
	H <- s2nr.value(K)
	if (verbose) cat('+++ initial  S2N ratio ...', formatC(H, format='e', digits=4, width=12), formatC(length(K), width=4), '\n')

	# start recursive merging by minimum S2NR-loss
	step.num <- 1
	step.max <- ifelse(k < 0, abs(k), bdm$wtt[[layer]]$s -k)
	steps <- matrix(rep(0, step.max *2), nrow = step.max)
	H <- rep(0, step.max)

	while (step.num <= step.max)
	{
		# get current list of child clusters
		K.childs <- which(!sapply(K, function(k) k$k) %in% unique(sapply(K, function(k) k$fthr)))

		# !!! check blocked state (i.e. any remaining cluster is its own father)
		# might represent intermediate robust clustering situations (??)
		if (length(K.childs) == 0)
		{
			if (verbose) cat('+++ blocked hierarchy level reached at', length(K), 'clusters !!!\n')
			if (verbose) cat('+++ next level ... \n')
			# reset data clustering assignment
			L <- merge.labels(Y, C, G)
			# reset clustering information
			K <- s2nr.kList(C, L, X, Z, G, P)
			# recompute list of child clusters
			K.childs <- which(!sapply(K, function(k) k$k) %in% unique(sapply(K, function(k) k$fthr)))
			if (length(K.childs)==0) break
		}

		# compute s2nr for all possible merges
		H.step <- s2nr.step(K, K.childs)
		# get index of next clusters to be merged
		i <- K.childs[which.max(H.step)]
		j <- which(sapply(K, function(k) k$k)==K[[i]]$fthr)
		# merge clusters
		C[which(C==K[[i]]$k)] <- K[[j]]$k

		# save merge information
		steps[step.num,] <- c(K[[i]]$k, K[[j]]$k)
		H[step.num] <- max(H.step)
		# show merge information
		if (verbose) cat('+++ merging', formatC(K[[i]]$k, width=4), 'into', formatC(K[[j]]$k, width=4), formatC(max(H.step), format='e', digits=4, width=12), formatC(length(K)-1, width=4), '\n')

		# update status of merged clusters in K
		K[[j]] <- s2nr.update(K[[j]], K[[i]])
		K[[i]] <- NULL

		step.num <- step.num +1
	}

	if (step.num < step.max) return(K)
	return(list(C = C, steps = steps, H = H))
}

s2nr.kList <- function(C, L, X, Z, G, P)
{
	max.Var <- apply(X, 2, var)
	K <- lapply(unique(C), function(k)
	{
		k.data <- which(L == k)
		k.list <- list(k = k, fthr = s2nr.kfthr(C, Z, G, P, k), N = length(k.data), M = rep(0, ncol(X)), V = max.Var)
		if (length(k.data) > 1){
			k.list$M <- apply(X[k.data,], 2, mean)
			k.list$V <- apply(X[k.data,], 2, var)
		} else if (length(k.data) == 1){
			k.list$M <- X[k.data,]
		}
		return(k.list)
	})
	return(K)
}

s2nr.kfthr <- function(C, Z, G, P, k.merged)
{
	# get nodes corresponding to k.merged
	k.nodes <- which(C == k.merged)
	# get k.father assuming a down-hill labelling
	# (as if the merged cluster had no valid peak)
	# get internal boundering nodes
	k.bound <- lapply(k.nodes, function(n)
	{
		n.bound <- grid_bound(n -1, G) +1
		# filter out all nodes belonging to k.merged
		# (keep only external boundering neighbours)
		n.bound <- n.bound[which(C[n.bound] != k.merged)]
		# if it is a boundering node
		# return it along with its highest external neighbour
		if (length(n.bound) > 0)
			return(list(n.int = n, n.ext = n.bound[which.max(Z[n.bound])]))
	})
	# filter out NULL elements
	k.bound <- k.bound[lapply(k.bound, length) > 0]
	# get highest internal boundering node
	n.highest <- which.max(Z[sapply(k.bound, function(k) k$n.int)])
	# get the cluster of its highest external neighbour
	k.father <- C[k.bound[[n.highest]]$n.ext]
	# check that k.father's peak is higher
	if (Z[P[k.father]] < Z[P[k.merged]]) k.father <- k.merged
	return(k.father)
}

# update cluster statistics after merging
s2nr.update <- function(kj, ki)
{
	N <- (kj$N + ki$N)
	M <- (kj$M *kj$N + ki$M *ki$N) /N
	V <- ((kj$N *(kj$V +(kj$M -M)**2)) +(ki$N *(ki$V +(ki$M -M)**2))) /N
	return(list(k = kj$k, fthr = kj$fthr, N = N, M = M, V = V))
}

# !!! Remember: \sum_{ij}((x_i-x_j)^2) = 2n \sum_i((x_i-\hat{x})^2)
# i.e. ... the variance is equivalent to the sum of pairwise distances among datapoints
# i.e ... the sum of squared-distances between cluster centroids and global mean is equivalent to the sum of pairwise distances among cluster centroids.
s2nr.value <- function(K)
{
	k.N <- sapply(K, function(k) k$N)
	k.valid <- which(k.N > 0)
	# signal: inter-cluster (explained) variance
	exp.var <- t(sapply(K, function(k) k$M))[k.valid, ]
	exp.var <- k.N[k.valid] *(exp.var - apply(exp.var, 2, mean))**2
	# noise: intra-cluster (unexplained) variance
	unx.var <- t(sapply(K, function(k) k$V))[k.valid, ]
	return(sum(exp.var) /sum(unx.var) /sum(k.N))
}

# compute S2NR for all possible merges at current step
s2nr.step <- function(K, K.childs)
{
	H.step <- sapply(K.childs, function(i)
	{
		# get index of k.father of child i
		j <- which(sapply(K, function(k) k$k)==K[[i]]$fthr)
		# save current status of merged clusters
		k.merged <- K[[i]]
		k.father <- K[[j]]
		# update merged cluster ststistics
		K[[j]] <- s2nr.update(K[[j]], K[[i]])
		K[[i]]$N <- 0
		# compute resulting s2n ratio
		h <- s2nr.value(K)
		# restore status of merged clusters
		K[[i]] <- k.merged
		K[[j]] <- k.father
		return(h)
	})
	return(H.step)
}

s2nr.info <- function(H, B = NULL)
{
	# H goes from initial partition with s = bdm$wtt[[layer]]$s clusters to 2 clusters
	# thus, length(H) is (s -1)
	s <- length(H) +1
	dffH <- c(diff(H), 0)
	optk <- sort(sort(-dffH, decreasing = TRUE, index.return = TRUE)$ix[1:length(B)])
	# !!! dffH[1] is the decrease in s2nr that results from the first merge
	# thus, the corresponding optimum is at H[optk = 1], i.e. s-1+1 == s
	if (!is.null(B)) cat('+++ hierarchy levels .... : ', formatC(unlist(B), width=4), '\n')
	cat('+++ S2NR significant loss : ', formatC(s -optk +1, width = 4), '\n')
	return(list(heuristic = 's2nr', H = H, lvls = unlist(B), loss = s- optk +1))
}


# -----------------------------------------------------------------------------
# +++ find optimal number of lcusters (by recursive merging based on S2NR)
# -----------------------------------------------------------------------------

s2nr.optk <- function(bdm, verbose = T, layer = 1)
{
	# input data
	X <- Xdata.get(bdm$data, bdm$Xdata$whiten, bdm$Xdata$input.dim, bdm$is.distance)$X
	# output data
	l <- c(1, 2) + (layer -1) *2
	Y <- bdm$ptsne$Y[, l]
	# density estimation
	Z <- bdm$pakde[[layer]]$z
	# grid parameters
	G <- bdm$wtt[[layer]]$grid
	# density landscape peaks
	P <- bdm$wtt[[layer]]$P
	# base cell cluster assignment
	C <- bdm$wtt[[layer]]$C
	L <- merge.labels(Y, C, G)

	# base clustering information
	if (verbose) cat('+++ computing clustering hierarchy ...\n')
	K <- s2nr.kList(C, L, X, Z, G, P)

	# initialize s2n vector, blocked states list
	H <- rep(0, length(K) -1)
	B <- list()
	# base S2NR
	H[1] <- s2nr.value(K)
	# show base information
	if (verbose) cat('+++ initial S2N ratio ... ', formatC(H[1], format='e', digits=4, width=12), formatC(length(K), width=4), '\n')

	# start recursive merging by minimum S2NR-loss
	for (step in seq(2, length(K) -1))
	{
		# get current list of child clusters
		K.childs <- which(!sapply(K, function(k) k$k) %in% unique(sapply(K, function(k) k$fthr)))
		# !!! check blocked state (i.e. any remaining cluster is its own father)
		# might represent intermediate robust clustering situations (??)
		if (length(K.childs)==0)
		{
			if (verbose) cat('+++ blocked hierarchy level reached at', length(K), 'clusters !!!\n')
			B[[length(B) +1]] <- length(K)
			if (verbose) cat('+++ next level ... \n')
			# reset data clustering assignment
			L <- merge.labels(Y, C, G)
			# reset clustering information
			K <- s2nr.kList(C, L, X, Z, G, P)
			# recompute list of child clusters
			K.childs <- which(!sapply(K, function(k) k$k) %in% unique(sapply(K, function(k) k$fthr)))
			if (length(K.childs) == 0) break
		}

		# compute S2NR for all possible merges
		H.step <- s2nr.step(K, K.childs)
		# get index of clusters to be merged at this step
		i <- K.childs[which.max(H.step)]
		j <- which(sapply(K, function(k) k$k) == K[[i]]$fthr)
		# merge clusters
		C[which(C == K[[i]]$k)] <- K[[j]]$k
		# save s2n ratio
		H[step] <- max(H.step)

		# show merge information
		if (verbose) cat('+++ merging', formatC(K[[i]]$k, width=4), 'into', formatC(K[[j]]$k, width=4), formatC(H[step], format='e', digits=4, width=12), formatC(length(K)-1, width=4), '\n')

		# change status of merged clusters
		K[[j]] <- s2nr.update(K[[j]], K[[i]])
		K[[i]] <- NULL

		}

	optk <- s2nr.info(H, B)
	return(optk)

}
