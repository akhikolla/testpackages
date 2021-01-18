# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# -----------------------------------------------------------------------------
# +++ distributed ptSNE implementation with shared memory
# -----------------------------------------------------------------------------

ptsne.get <- function(bdm, cl, Y.init = NULL, info = 0)
{
	cat('+++ computing ptSNE \n')
	if (attr(cl[[1]], 'class') == 'SOCKnode')
		sckt.ptsne(bdm, cl, Y.init = Y.init, progress = info)
	else
		mpi.ptsne(bdm, cl, Y.init = Y.init, progress = info)
}

sckt.ptsne <- function(bdm, cl, Y.init = NULL, progress = 0)
{
	# setup parameters
	threads <- bdm$ptsne$threads
	layers <- bdm$ptsne$layers
	rounds <- bdm$ptsne$rounds
	alpha <- bdm$ptsne$alpha

	chnk.brks <- round(seq(1, nrow(bdm$data)+1, length.out=(threads+1)), 0)
	thrd.size <- round(mean(diff(chnk.brks)) * layers, 0)
	iters <- ceiling(sqrt(thrd.size) *bdm$ptsne$boost)
	epochs <- floor(sqrt(nrow(bdm$data)) /bdm$ptsne$boost)

	# the next condition is meant only for paper pourposes
	# (must be deleted before building the package)
	if (!is.null(bdm$force)) epochs <- bdm$force
	rxEpochs <- rounds *epochs

	# export setup parameters
	clusterExport(cl, c('layers', 'iters', 'alpha'), envir = environment())
	# clusterExport(cl, c('progress'), envir = environment())

	# initial mapping (random circular embedding of radius 1)
	if (is.null(Y.init)) Y.init <- ptsne.init(nrow(bdm$data), layers)
	Ybm <- as.big.matrix(Y.init, type='double')
	Ybm.dsc <- describe(Ybm)

	# row sampling indexes (substract 1 to convert to C++ indexes !!)
	Ibm <- big.matrix(nrow(bdm$data), 1, type='integer')
	Ibm.dsc <- describe(Ibm)

	# epoch/thread cost matrix
	Cbm <- as.big.matrix(matrix(0, threads, (rxEpochs +1)), type='double')
	Cbm.dsc <- describe(Cbm)

	# epoch/thread size matrix
	Sbm <- as.big.matrix(matrix(0, layers, (rxEpochs +1)), type='double')
	Sbm.dsc <- describe(Sbm)

	# attach bigmatrices to workers
	clusterExport(cl, c('Ybm.dsc'), envir=environment())
	clusterEvalQ(cl, Ybm <- attach.big.matrix(Ybm.dsc))
	clusterExport(cl, c('Ibm.dsc'), envir=environment())
	clusterEvalQ(cl, Ibm <- attach.big.matrix(Ibm.dsc))
	clusterExport(cl, c('Cbm.dsc'), envir=environment())
	clusterEvalQ(cl, Cbm <- attach.big.matrix(Cbm.dsc))
	clusterExport(cl, c('Sbm.dsc'), envir=environment())
	clusterEvalQ(cl, Sbm <- attach.big.matrix(Sbm.dsc))

	# special initialization for thread.rank == 0
	clusterEvalQ(cl,
		if (thread.rank == 0) {
			w <- new.env()
			w$mapp.list <- list()
		})

	# compute initial cost
	iters <- 0
	clusterExport(cl, c('iters'), envir=environment())
	# resample dataset
	Ibm[ ] <- sample(seq(0, nrow(bdm$data)-1))
	# perform ptSNE
	nulL <- clusterCall(cl, sckt.ztsne, 0)

	# reset number of iterations
	iters <- ceiling(sqrt(thrd.size) *bdm$ptsne$boost)
	clusterExport(cl, c('iters'), envir=environment())

	# starting time
	t0 <- Sys.time()

	# report starting information
	avgCost <- mean(Cbm[, 1])
	avgSize <- mean(Sbm[, 1])
	if (progress >= 0) {
		ptsne.info(threads, thrd.size, rxEpochs, iters, 0, avgCost, avgSize, t0)
	}

	for (r in seq(rounds)){

		if (progress >= 0) {
			cat('--- round ', formatC(r, width=2, flag='0'), '/', formatC(rounds, width=2, flag='0'), sep='')
			cat('     ')
			cat('   <Cost>   ')
			cat('   <Size>   ')
			cat('   epSecs ')
			cat(' time2End ')
			cat(' eTime ')
			cat('\n')
		}

		for (e in seq((r-1)*epochs+1, (r*epochs))) {
			# epoch start-time
			te <- Sys.time()
			# resample dataset
			Ibm[ ] <- sample(seq(0, nrow(bdm$data)-1))
			# perform ptSNE
			nulL <- clusterCall(cl, sckt.ztsne, e)
			# security break control
			if (mean(Cbm[, e+1]) <0) break
			# report status
			if (progress >=0) {
				avgCost <- mean(Cbm[, e+1])
				avgSize <- mean(Sbm[, e+1])
				epoch.info(e, rxEpochs, avgCost, avgSize, t0, te)
			}

		}

		# save round
		if (progress >= 1 || r == rounds){
			bdm$ptsne$rounds <- r
			bdm$ptsne$Y <- as.matrix(Ybm[ , ])
			bdm$ptsne$cost <- as.matrix(Cbm[, 1:(r*epochs +1)])
			bdm$ptsne$size <- as.matrix(Sbm[, 1:(r*epochs +1)])
		}
		if (progress == 1) bdm.scp(bdm)

		# security break control
		if (mean(Cbm[, (r*epochs +1)]) <0 ) break

	}

	# if (progress == 2) {
	# 	bdm$progress <- clusterEvalQ(cl, if (thread.rank == 0) w$mapp.list)[[1]]
	# }

	# report status
	avgCost <- mean(Cbm[, e +1])
	avgSize <- mean(Sbm[, e +1])

	ptsne.info(threads, thrd.size, rxEpochs, iters, e, avgCost, avgSize, t0)

	return(bdm)

}

# -----------------------------------------------------------------------------
# +++ ptSNE SCKT thread functions
# -----------------------------------------------------------------------------

sckt.ztsne <- function(epoch)
{
	if (thread.rank == 0) {
		# reduce embedding size
		Sbm[, (epoch +1)] <- sapply(seq(layers), function(l) {
			lc <- (l-1)*2 +1
			sqrt(diff(range(Ybm[, lc]))**2 + diff(range(Ybm[, (lc+1)]))**2)
		})
		# # save current embedding to make movie
		# if (progress == 2) {
		# 	w$mapp.list[[(epoch +1)]] <- list(epoch = (epoch +1), Y = as.matrix(Ybm[ , 1:2]))
		# }
	}
	else {
		Cbm[thread.rank, (epoch +1)] <- sckt_zTSNE(thread.rank, threads, layers, Xbm@address, Bbm@address, Ybm@address, Ibm@address, iters, alpha, is.distance)
	}
}

# -----------------------------------------------------------------------------
# +++ ptSNE MPI
# -----------------------------------------------------------------------------

mpi.ptsne <- function(bdm, cl, Y.init = NULL, progress = 0)
{
	# setup parameters
	threads <- bdm$ptsne$threads
	layers <- bdm$ptsne$layers
	rounds <- bdm$ptsne$rounds
	alpha <- bdm$ptsne$alpha

	# export ptSNE setup parameters
	clusterExport(cl, c('alpha'), envir = environment())

	zSize <- round(nrow(bdm$data) /threads *layers, 0)
	iters <- ceiling(sqrt(zSize) *bdm$ptsne$boost)
	epochs <- floor(sqrt(nrow(bdm$data)) /bdm$ptsne$boost)

	# the next condition is meant only for paper pourposes
	# (must be deleted before building the package)
	if (!is.null(bdm$force)) epochs <- bdm$force
	rxEpochs <- rounds *epochs

	# chunking
	chnk.brks <- round(seq(1, nrow(bdm$data) +1, length.out = threads +1), 0)
	# thread segments: row/col indexes in Y
	thrd.brks <- lapply(seq(threads), function(z)
	{
		t(sapply(seq(layers), function(l) {
			a <- z + (l-1)
			if (a > threads) a <- a %% threads
			# Att!! C++ indexes
			c(chnk.brks[a], chnk.brks[a+1]) -1
		}))
	})

	# initialize Z.list to map I and Y chunks to workers
	Z.list <- lapply(seq(threads), function(z) {
		matrix(0, sum(apply(thrd.brks[[z]], 1, diff)), 3)
	})

	# initial mapping (random circular embedding of radius 1)
	if (is.null(Y.init)) Y.init <- ptsne.init(nrow(bdm$data), layers)
	Y <- Y.init

	# initialize cost&size functions
	e.cost <- matrix(0, nrow = threads, ncol = rxEpochs +1)
	e.size <- matrix(0, nrow = layers, ncol = rxEpochs +1)

	# special initialization for thread.rank != 0
	clusterEvalQ(cl,
		if (thread.rank != 0) {
			w <- new.env()
			w$zI <- numeric()
			w$zY <- numeric()
		})

	# +++ compute initial cost&size
	clusterEvalQ(cl, iters <- 0)
	# resample dataset. Att!! C++ indexes
	I <- sample(seq(nrow(bdm$data))) -1
	# get I&Y chunks
	zChnks(Z.list, Y, I, thrd.brks)
	# pool cost from workers
	# NULL returned by thread.rank == 0 is removed by unlist
	e.cost[, 1] <- unlist(clusterApply(cl, Z.list, mpi.ztsne))
	# compute layers' size
	e.size[ , 1] <- eSize(Y)
	# +++ set number of iterations on workers
	clusterExport(cl, c('iters'), envir = environment())

	# starting time
	t0 <- Sys.time()

	# report starting information
	avgCost <- mean(e.cost[ , 1])
	avgSize <- mean(e.size[ , 1])
	ptsne.info(threads, zSize, rxEpochs, iters, 0, avgCost, avgSize, t0)

	for (r in seq(rounds)){

		cat('--- round ', formatC(r, width=2, flag='0'), '/', formatC(rounds, width=2, flag='0'), sep='')
		cat('     ')
		cat('   <Qlty>   ')
		cat('   <Size>   ')
		cat(' <epSecs> ')
		cat(' time2End ')
		cat(' eTime ')
		cat('\n')

		for (e in seq((r-1)*epochs+1, (r*epochs))){

			# epoch starting time
			te <- Sys.time()

			# resample dataset. Att!! C++ indexes
			I <- sample(seq(nrow(bdm$data))) -1
			# get I&Y chunks
			zChnks(Z.list, Y, I, thrd.brks)

			# perform partial ptSNE pooling workers epoch-cost
			# NULL returned by thread.rank == 0 is removed by unlist
			e.cost[, (e +1)] <- unlist(clusterApply(cl, Z.list, mpi.ztsne))

			# pool partial mappings from workers
			zMap.list <- clusterEvalQ(cl, if (thread.rank != 0) w$zY)
			# restructure global mapping
			updateY(Y, I, zMap.list, thrd.brks)
			# compute epoch's embedding size
			e.size[, (e+1)] <- eSize(Y)

			# report status
			if (progress >= 0) {
				avgCost <- mean(e.cost[, e+1])
				avgSize <- mean(e.size[, e+1])
				epoch.info(e, rxEpochs, avgCost, avgSize, t0, te)
			}

			# security break control
			if (mean(e.cost[, e+1]) < 0) break
		}

		# save round
		if (progress >= 1 || r == rounds) {
			bdm$ptsne$rounds <- r
			bdm$ptsne$Y <- Y
			bdm$ptsne$cost <- as.matrix(e.cost[, 1:(r *epochs +1)])
			bdm$ptsne$size <- as.matrix(e.size[, 1:(r *epochs +1)])
		}
		if (progress == 1) bdm.scp(bdm)

	}

	# report status
	avgCost <- mean(e.cost[, e +1])
	avgSize <- mean(e.size[, e +1])
	ptsne.info(threads, zSize, rxEpochs, iters, e, avgCost, avgSize, t0)

	return(bdm)

}

# -----------------------------------------------------------------------------
# +++ ptSNE MPI worker functions
# -----------------------------------------------------------------------------

mpi.ztsne <- function(zChnk)
{
	if (thread.rank != 0) {
		w$zI <-  zChnk[, 1]
		w$zY <- zChnk[, 2:3]
		mpi_zTSNE(Xbm@address, Bbm@address, w$zY, w$zI, iters, alpha, is.distance)
	}
}

# -----------------------------------------------------------------------------
# +++ auxiliary functions for ptSNE
# -----------------------------------------------------------------------------

# +++ initial mapping in a random circular embedding of radius 1
ptsne.init <- function(N, layers)
{
	# relation square/circular embedding area (2*r)^2/(pi*r^2) = 4/pi
	# 5/pi ensures we have enough initial positions
	xN <- ceiling(N * 5/pi)
	Y <- matrix(runif(xN*2*layers, -1, 1), nrow = xN*layers)
	I <- which(apply(Y**2, 1, sum) < 1)[1:(N*layers)]
	return(matrix(unlist(t(Y[I, ])), nrow = N, byrow = T))
}

time.format <- function(time.secs)
{
	time.secs <- round(time.secs, 0)
	if (time.secs > 86400)
	{
		dd <- time.secs %/% 86400
		hh <- (time.secs - dd*86400) %/% 3600
		mm <- (time.secs - dd*86400 - hh*3600) %/% 60
		ft <- paste(formatC(dd, width=2, flag='0'), 'D', formatC(hh, width=2, flag='0'), ':', formatC(mm, width=2, flag='0'), sep='')
	} else
	{
		hh <- time.secs %/% 3600
		mm <- (time.secs - hh*3600) %/% 60
		ss <- (time.secs - hh*3600 - mm*60)
		ft <- paste(formatC(hh, width=2, flag='0'), ':', formatC(mm, width=2, flag='0'), ':', formatC(ss, width=2, flag='0'), sep='')
	}
	return(ft)
}

ptsne.info <- function(threads, zSize, rxEpochs, iters, e, avgCost, avgSize, t0)
{
	cat('...')
	cat(' threads ', threads, sep='')
	cat(', size ', zSize, sep='')
	cat(', epochs ', rxEpochs, sep='')
	cat(', iters ', iters, sep='')
	cat('\n')
	cat('+++ epoch ', formatC(e, width=4, flag='0'), '/', formatC(rxEpochs, width=4, flag='0'), sep='')
	cat(formatC(avgCost, format='e', digits=4, width=12))
	cat(formatC(avgSize, format='e', digits=4, width=12))
	if (e > 0){
		runTime <- as.numeric(difftime(Sys.time(), t0, units='secs'))
		cat('  ', formatC(runTime/rxEpochs, format='f', digits=4, width=8), sep='')
		cat('  ', time.format(runTime), sep='')

	}
	cat('\n')
}

epoch.info <- function(e, rxEpochs, avgCost, avgSize, t0, te)
{
	cat('+++ epoch ', formatC(e, width=4, flag='0'), '/', formatC(rxEpochs, width=4, flag='0'), sep='')
	cat(formatC(avgCost, format='e', digits=4, width=12))
	cat(formatC(avgSize, format='e', digits=4, width=12))
	epchSecs <- as.numeric(difftime(Sys.time(), te, units='secs'))
	cat('  ', formatC(epchSecs, format='f', digits=4, width=8), sep='')
	runTime <- as.numeric(difftime(Sys.time(), t0, units='secs'))
	tm2End <- runTime/e * (rxEpochs-e)
	cat('  ', time.format(tm2End), sep='')
	cat('  ', format(Sys.time()+tm2End, '%H:%M'), sep='')
	cat('\n')
}
