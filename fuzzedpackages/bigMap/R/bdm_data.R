# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Process raw input data

Xdata.get <- function(X, whiten = 4, input.dim = NULL, is.distance = F, quiet = FALSE)
{
	if (!quiet) cat('+++ processing data \n')
	if (is.null(input.dim))	input.dim <- ncol(X)

	if (is.distance){
	    return(list(X = X, whiten = 0, input.dim = 1))
	}
	else if (whiten == 0){
	    return(list(X = X[, 1:input.dim], whiten = 0, input.dim = input.dim))
	}
	else {
		X <- data.get(X, whiten, input.dim)
		return(list(X = X, whiten = whiten, input.dim = input.dim))
	}
}


# -----------------------------------------------------------------------------
# +++ Preprocessing of input-data
# -----------------------------------------------------------------------------

data.get <- function(X, whiten, input.dim)
{
	# filter out all irrelevant features
	# (makes sense in datasets like the mnist.optical.digits where some features might have zeros for all observations)
	X <- X[, which(apply(X, 2, sum) != 0)]
	input.dim <- min(input.dim, ncol(X))

	if (whiten == 1)	# centering
	{
		X <- scale(X, center = T, scale = F)
	}
	else if (whiten == 2)	# centering & scaling
	{
		X <- scale(X, center = T, scale = T)
		if (any(is.na(X))) {
			return(message('+++ Error: scaling return NaNs !!!'))
		}
	}
	else if (whiten %in% c(3, 4))	# PCA/whitening
	{
		X <- scale(X, center = T, scale = F)
		# covariance matrix
		V <- t(X) %*% X / nrow(X)
		# singular value decomposition
		s <- La.svd(V)
		# whitening
		if (whiten == 4) s$u <- s$u %*% diag(1/sqrt(s$d))
		# take first input.dim dimensions
		X <- X[, 1:input.dim] %*% s$u[1:input.dim, ]
	}

	return(X)
}

# -----------------------------------------------------------------------------
# +++ Export input-data to workers
# -----------------------------------------------------------------------------

Xdata.exp <- function(cl, X, is.distance)
{
	clusterExport(cl, c('is.distance'), envir = environment())
	if (attr(cl[[1]], 'class') == 'SOCKnode')
	{
		# input-data big.matrix
		Xbm <- as.big.matrix(X, type = 'double')
		Xbm.dsc <- describe(Xbm)
		# export big matrix descriptor to workers
		clusterExport(cl, c('Xbm.dsc'), envir = environment())
		# attach big matrix to workers
		clusterEvalQ(cl, Xbm <- attach.big.matrix(Xbm.dsc))
	}
	else
	{
		f <- tName.get('X')
		Xbf <- as.big.matrix(X, type='double', backingpath = f$path, backingfile = f$bin, descriptorfile = f$desc)
		Xbf.dsc <- describe(Xbf)
		clusterExport(cl, c('Xbf.dsc'), envir = environment())
		# attach big.matrix backing.file to holders
		nulL <- clusterEvalQ(cl,
			if (thread.rank == thread.hldr) {
				Xhl <- attach.big.matrix(Xbf.dsc)
				Xbm <- as.big.matrix(as.matrix(Xhl[ , ]), type='double')
				rm(Xhl)
			})
		# get big.matrix backing.file descriptors from holders
		cl.Xdsc <- clusterEvalQ(cl,
			if (thread.rank == thread.hldr) {
				describe(Xbm)
			})
		# export shared-memory descriptors
		clusterExport(cl, c('cl.Xdsc'), envir = environment())
		# attach big matrix to workers
		nulL <- clusterEvalQ(cl,
			if (thread.rank != thread.hldr){
				Xbm <- attach.big.matrix(cl.Xdsc[[thread.hldr]])
			})
		# remove backing file
		unlink(paste(f$path, f$bin, sep = '/'))
		unlink(paste(f$path, f$desc, sep = '/'))
	}
}
