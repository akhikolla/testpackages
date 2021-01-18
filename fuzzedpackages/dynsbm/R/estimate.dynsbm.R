# This file is part of dynsbm.

# dysbm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# dynsbm is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with dynsbm.  If not, see <http://www.gnu.org/licenses/>
		
estimate.dynsbm <- function (Y, present=NULL, Q,
                             edge.type=c("binary","discrete","continuous"), K=-1,
                             directed=FALSE, self.loop=FALSE,
                             init.cluster=NULL,
                             nb.cores=1,
                             perturbation.rate=0., iter.max=20,
                             fixed.param=FALSE,
                             bipartition=NULL) {
    T <- dim(Y)[1]
    N <- dim(Y)[2]
    if (is.null(present)){
        present <- matrix(1L,dim(Y)[2],dim(Y)[1])
    } else{
        if (storage.mode(present)!= "integer") storage.mode(present) <- "integer"
    }
    is.bipartite = FALSE
    if(!is.null(bipartition))
        is.bipartite = TRUE
    ##---- initialization
    if(is.null(init.cluster)){
        ## aggregation
        Yaggr <- c()
        if(!any(present==0L)){ ## concat all time step if no absent nodes
            for (t in 1:T){
                if(directed){
                    Yaggr <- cbind(Yaggr,Y[t,,],t(Y[t,,]))
                } else{
                    Yaggr <- cbind(Yaggr,Y[t,,])
                }
            }
        } else{ ## aggregate all time step if absent nodes
            Yaggr <- Y[1,,]
            for (t in 2:T){
                Yaggr <- Yaggr+Y[t,,]
            }
            nbcopresent <- present %*% t(present)
            nbcopresent[nbcopresent==0] <- -1 # to avoid NaN
            Yaggr <- Yaggr / nbcopresent
            if(directed) Yaggr <- cbind(Yaggr,t(Yaggr))
        }
        if(is.bipartite){ 
            ## bipartite case       
            set1 <- which(bipartition==1)
            set2 <- which(bipartition==2)
            if(length(set1)>length(set2)){
                Q1 <- ceiling(Q/2.) 
                Q2 <- floor(Q/2.) 
            } else{
                Q1 <- floor(Q/2.) 
                Q2 <- ceiling(Q/2.) 
            }
            if (nb.cores==1){
                km1 <- kmeans(Yaggr[set1,], Q1, nstart=64)
                km2 <- kmeans(Yaggr[set2,], Q2, nstart=64)
                init.cluster <- c(km1$cluster, Q1+km2$cluster)
            } else{
                RNGkind("L'Ecuyer-CMRG")
                pkm1 <- mclapply(rep(64/nb.cores,nb.cores), function(nstart) kmeans(Yaggr[set1,], Q1, nstart=nstart), mc.cores=nb.cores)
                i1 <- sapply(pkm1, function(result) result$tot.withinss)
                pkm2 <- mclapply(rep(64/nb.cores,nb.cores), function(nstart) kmeans(Yaggr[set2,], Q2, nstart=nstart), mc.cores=nb.cores)
                i2 <- sapply(pkm2, function(result) result$tot.withinss)
                init.cluster <- c(pkm1[[which.min(i1)]]$cluster, Q1+pkm2[[which.min(i2)]]$cluster)
            } 
            ## perturbation
            if(perturbation.rate>0.){
                N1 <- length(set1)
                v <- sample(1:N1, perturbation.rate*N1)
                vs <- sample(v)
                init.cluster[set1[v]] <- init.cluster[set1[vs]]
                N2 <- length(set2)
                v <- sample(1:N2, perturbation.rate*N2)
                vs <- sample(v)
                init.cluster[set2[v]] <- init.cluster[set2[vs]]
            }
        } else{
            ## classical case
            if (nb.cores==1){
                km <- kmeans(Yaggr, Q, nstart=64)
                init.cluster <- km$cluster
            } else{
                RNGkind("L'Ecuyer-CMRG")
                pkm <- mclapply(rep(64/nb.cores,nb.cores), function(nstart) kmeans(Yaggr, Q, nstart=nstart), mc.cores=nb.cores)
                i <- sapply(pkm, function(result) result$tot.withinss)
                init.cluster <- pkm[[which.min(i)]]$cluster
            }
            ## perturbation
            if(perturbation.rate>0.){
                v <- sample(1:N, perturbation.rate*N)
                vs <- sample(v)
                init.cluster[v] <- init.cluster[vs]
            }
        }
    }
    ##---- estimation
    result <- list(
        init.cluster=init.cluster,
        dynsbm=dynsbmcore(T, N, Q, as.vector(Y), present, edge.type, K, init.cluster-1, iter.max, nb.cores, directed, self.loop, fixed.param) # -1 for compatibility with C++   
    )
    ##---- verifying bipartition
    if(is.bipartite){ 
        bipartite.check <- apply(result$dynsbm$membership, 2, function(mbt){
            crosstable <- table(mbt[which(mbt>0)],bipartition[which(mbt>0)])
            return(length(which(crosstable>0))>Q)
        })
        if(sum(bipartite.check)>0) warning("Model estimation is not coherent with bipartition: group membership is not bipartite\n")
    }
    
    ##---- correcting trans in the case of classical SBM
    if(T==1){
        result$dynsbm$trans <- matrix(0,Q,Q)
        diag(result$dynsbm$trans) <- 1
    }
    return(result)
}
