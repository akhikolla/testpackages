# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Initialize parallel computing environment.
#
# @param threads The number of parallel threads (in principle only limited by hardware resources, \code{i.e.} number of cores and available memory)
#
# @param type The type of cluster: 'SOCK' (default) for intranode parallelization, 'MPI' for internode parallelization (\code{message passing interface} parallel environment).
#
# @return cl A cluster instance (as created by the snow::makeCluster() function).

cluster.start <- function(threads, type = 'SOCK')
{
	cl <- NULL
	if (type == 'MPI')
	{
		cl <- makeCluster()
		# cluster structure
		cl.ranks <- unlist(clusterEvalQ(cl, thread.rank <- Rmpi::mpi.comm.rank(comm = 0)))
		cl.nodes <- unlist(clusterEvalQ(cl, thread.node <- Rmpi::mpi.get.processor.name()))
		cl.hldrs <- sapply(unique(cl.nodes), function(node) min(cl.ranks[which(cl.nodes == node)]))
		# assign holders to workers
		clusterExport(cl, c('cl.hldrs'), envir = environment())
		nulL <- clusterEvalQ(cl, thread.hldr <- cl.hldrs[Rmpi::mpi.get.processor.name()])
	}
	else
	{
		cl <- makeCluster(threads+1, type = 'PSOCK')
		# cluster structure
		clusterApply(cl, seq_along(cl), function(i) thread.rank <<- (i-1))
	}
	if (!is.null(cl))
	{
		clusterExport(cl, c('threads'), envir = environment())
		cat('+++ running ', threads, ' threads \n', sep = '')
		# load workers environment
		clusterEvalQ(cl, library(bigmemory))
		clusterEvalQ(cl, library(bigMap))
	}
	if (substr(bdm.local(), 1, 7) == 'xxx.xxx'){
		cat('+++ WARNING: bdm.local() not set !! \n')
	}
	return(cl)
}

# Close parallel computing environment.
#
# @param cl A cluster instance (as created by the snow::makeCluster() function).
#
# @return None

cluster.stop <- function(cl)
{
	if (!is.null(cl)) stopCluster(cl)
}
