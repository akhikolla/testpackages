# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Perplexity-based local betas

Xbeta.get <- function(cl, X, ppx = 100, itr = 100, tol = 1e-05, is.distance = F, quiet = FALSE)
{
	if (!quiet) cat('+++ Computing Betas, perplexity ', ppx, ' \n', sep='')
	if (!all(unlist(clusterEvalQ(cl, ('Xbm' %in% ls()))))) {
		Xdata.exp(cl, X, is.distance)
	}
	B <- beta.get(cl, ppx, itr, tol)
	# (beta is 1/(2*sigma**2) !!!)
	if (!quiet) print(summary(1/sqrt(2 *B)))
	return(list(B = B, ppx = ppx, itr = itr, tol = tol))
}

# -----------------------------------------------------------------------------
# +++ Parallelized computation of perplexity-based local betas (by chunks)
# -----------------------------------------------------------------------------

# Compute perplexity-based local betas. This function assumes that \var{Xbm} (input-data big.matrix) is attached to the workers of \var{cl}. Also \var{is.distance} must have been exported to the workers. The up-stream step \code{bdm.data()} works out both conditions.

# Att!!
# beta_i = 1/(2 *sigma**2) and sigma_i = 1/sqrt(2 *beta_i)
# thus, the BIVARIATE normal density function:
# pdf(dij) = 1/(2 *pi *sigma_i^2) exp(-1/2 * dij^2/sigma_i^2)
# is writen in terms of betak as:
# pdf(dij) = 1/pi *beta_i exp(- beta_i * dij^2)

beta.get <- function(cl, ppx, itr, tol)
{
	# export parameters
	clusterExport(cl, c('ppx', 'itr', 'tol'), envir=environment())
	# get perplexity-based local betas by chunks
	beta.list <- clusterCall(cl, thread.beta)
	return(unlist(beta.list))
}

# -----------------------------------------------------------------------------
# +++ worker function
# -----------------------------------------------------------------------------

thread.beta <- function(){
	if (thread.rank != 0) {
		zBeta(thread.rank, threads, Xbm@address, is.distance, ppx, tol, itr)
	}
}

# -----------------------------------------------------------------------------
# +++ Export beta to workers
# -----------------------------------------------------------------------------

Xbeta.exp <- function(cl, B)
{
	if (attr(cl[[1]], 'class') == 'SOCKnode')
	{
		# define Bbm big.matrix (betas)
		Bbm <- big.matrix(length(B), 1, type = 'double')
		Bbm[ ] <- B
		Bbm.dsc <- describe(Bbm)
		# export big matrix descriptor to workers
		clusterExport(cl, c('Bbm.dsc'), envir = environment())
		# attach big matrix to workers
		clusterEvalQ(cl, Bbm <- attach.big.matrix(Bbm.dsc))
	}
	else
	{
		f <- tName.get('B')
		Bbf <- big.matrix(length(B), 1, type='double', backingpath = f$path, backingfile = f$bin, descriptorfile = f$desc)
		Bbf.dsc <- describe(Bbf)
		Bbf[ ] <- B
		clusterExport(cl, c('Bbf.dsc'), envir = environment())
		# attach big.matrix backing.file to holders
		nulL <- clusterEvalQ(cl,
		if (thread.rank == thread.hldr) {
			Bbm <- attach.big.matrix(Bbf.dsc)
		})
		# get big.matrix backing.file descriptors from holders
		cl.Bdsc <- clusterEvalQ(cl,
		if (thread.rank == thread.hldr) {
			describe(Bbm)
		})
		# export sharedmemory descriptors
		clusterExport(cl, c('cl.Bdsc'), envir = environment())
		# attach big matrix to workers
		nulL <- clusterEvalQ(cl,
		if (thread.rank != thread.hldr){
			Bbm <- attach.big.matrix(cl.Bdsc[[thread.hldr]])
		})
		# remove backing file
		unlink(paste(f$path, f$bin, sep = '/'))
		unlink(paste(f$path, f$desc, sep = '/'))
	}
}
