# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# -----------------------------------------------------------------------------
# +++ Perplexity-Adaptive Kernel Density Estimation (paKDE)
# -----------------------------------------------------------------------------

pakde.get <- function(cl, Y, ppx, itr, tol, g, g.exp)
{
	N <- nrow(Y)

	# get perplexity-based local Betas
	# Xbeta.get() controls for Xbm (say Y here) being attached to all workers
	Xbeta.list <- Xbeta.get(cl, Y, ppx = ppx, itr = itr, tol = tol, is.distance = FALSE)
	B <- Xbeta.list$B
	# export Betas to workers
	Xbeta.exp(cl, B)

	# get pakde grid
	pakde.grid <- pakde.grid.get(Y, B, g, g.exp)

	# export grid & beta to nodes
	clusterExport(cl, c('pakde.grid'), envir=environment())

	cat(' computing grid cell densities ... \n')
	z <- unlist(clusterCall(cl, thread.pakde))
	# cell-wise sum of density contributions from all chunks (i.e. all sample-points)
	z <- apply(matrix(z, nrow=g**2), 1, sum) * pakde.grid$d / pi / N
	cat('+++ cdf ', round(sum(z),4), '\n', sep='')

	layer.pakde <- list(ppx = ppx, beta = B, g.exp = g.exp, x = pakde.grid$x, y = pakde.grid$y, z = matrix(z, nrow=g, byrow=TRUE))

	return(layer.pakde)

}


# -----------------------------------------------------------------------------
# +++ paKDE - Auxiliary functions
# -----------------------------------------------------------------------------

# +++ kernel density grid
pakde.grid.get <- function(Y, beta, g, g.exp){
	# expand grid limits to avoid border effects;
	# enclose the kernel densities of the most extreme mapped-points up to g.exp
	idx <- which.min(Y[ ,1])
	x.min <- Y[idx, 1] - ifelse(beta[idx]>0, g.exp/sqrt(2*beta[idx]), 0)
	idx <- which.max(Y[, 1])
	x.max <- Y[idx, 1] + ifelse(beta[idx]>0, g.exp/sqrt(2*beta[idx]), 0)
	idx <- which.min(Y[, 2])
	y.min <- Y[idx, 2] - ifelse(beta[idx]>0, g.exp/sqrt(2*beta[idx]), 0)
	idx <- which.max(Y[, 2])
	y.max <- Y[idx, 2] + ifelse(beta[idx]>0, g.exp/sqrt(2*beta[idx]), 0)
	# pakde.grid
	pakde.grid <- list()
	pakde.grid$x <- seq(x.min, x.max, length.out=g)
	pakde.grid$y <- seq(y.min, y.max, length.out=g)
	pakde.grid$d <- (pakde.grid$x[2]-pakde.grid$x[1]) * (pakde.grid$y[2]-pakde.grid$y[1])
	return(pakde.grid)
}

# compute mapped densities
thread.pakde <- function()
{
	chnk.brks <- round(seq(1, nrow(Xbm)+1, length.out=(threads+2)), 0)
	brk <- thread.rank +1
	z <- rep(0, length(pakde.grid$x) * length(pakde.grid$y))
	for (i in chnk.brks[brk]:(chnk.brks[brk+1] -1))
	{
		# grid cell distances to sample-point-i
		sqdst2i <- as.numeric(t(outer((Xbm[i, 1] - pakde.grid$x)^2, (Xbm[i, 2] - pakde.grid$y)^2, '+')))
		# grid cell densities by kernel-i
		z <- z + Bbm[i] * exp( -Bbm[i] * sqdst2i)
	}
	return(z)
}
