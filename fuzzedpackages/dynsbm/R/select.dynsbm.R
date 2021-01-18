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

select.dynsbm <- function(Y, present=NULL, Qmin, Qmax,
                          edge.type=c("binary","discrete","continuous"), K=-1,
                          directed=FALSE, self.loop=FALSE,
                          nb.cores=1,
                          iter.max=20, nstart=25, perturbation.rate=0.2,
                          fixed.param=FALSE,
                          bipartition=NULL,
                          plot=TRUE){
    if (is.null(present)){
        present <- matrix(0L,dim(Y)[2],dim(Y)[1])    
        for (t in 1:dim(Y)[1])
            present[union(which(apply(Y[t,,],1,FUN=function(v) sum(v>0))>0),which(apply(Y[t,,],2,FUN=function(v) sum(v>0))>0)),t] <- 1    		
    }
    never.present <- which(rowSums(present)==0L)
    if(length(never.present)){
        stop("Data format error: one or more nodes are never present.\nThis is not supported (see help about the present argument).\nPlease correct or remove the nodes from the adjacency matrices.")
    }
    if(!is.null(bipartition)){
        fixed.param <- TRUE
        if(length(bipartition)!=dim(Y)[2]){
            stop("Data format error: bipartition length is uncorrect.")
        }
    }
    list.dynsbm <- list()
    for (Q in Qmin:Qmax){
        results <- list()
        for (rep in 1:nstart){
            if (rep==1) this.perturbation.rate <- 0.0 else this.perturbation.rate <- perturbation.rate
            results[[length(results)+1]] <- estimate.dynsbm(Y=Y, present=present, Q=Q, directed=directed,
                                                            self.loop=self.loop,
                                                            edge.type=edge.type, K=K,
                                                            nb.cores=nb.cores, init.cluster=NULL,
                                                            perturbation.rate=this.perturbation.rate,
                                                            iter.max=5, fixed.param=fixed.param,
                                                            bipartition=bipartition)
        }
        best.result <- which.max(sapply(results, FUN=function(result) result$dynsbm$loglikelihood))
        best.init.cluster <- results[[best.result]]$init.cluster
        dynsbm <- estimate.dynsbm(Y=Y, present=present, Q=Q, directed=directed,
                                  self.loop=self.loop,
                                  edge.type=edge.type, K=K,
                                  nb.cores=nb.cores,
                                  init.cluster=best.init.cluster,
                                  iter.max=iter.max, fixed.param=fixed.param,
                                  bipartition=bipartition)$dynsbm
        class(dynsbm) <- c("list","dynsbm")
        list.dynsbm[[length(list.dynsbm)+1]] <- dynsbm
    }
    if(plot) plot.icl(list.dynsbm)
    list.dynsbm
}

plot.icl <- function(list.dynsbm){
    Qmin <- ncol(list.dynsbm[[1]]$trans)
    Qmax <- ncol(list.dynsbm[[length(list.dynsbm)]]$trans)
    logl <- sapply(list.dynsbm, FUN=function(model) model$loglikelihood)
    plot(Qmin:Qmax, logl, type='b', ylab="", xlab="Number of groups", yaxt='n')
    par(new=TRUE)
    if("gamma" %in% names(list.dynsbm[[1]])){
        legend("topright",legend=c("Loglikelihood", "ICL\n(not available)"), col=1:2, pch=19, lwd=4)
    } else{
        plot(Qmin:Qmax, sapply(list.dynsbm, compute.icl), col=2, type='b', ylab="", xlab="Number of groups", yaxt='n')
        legend("topright",legend=c("Loglikelihood","ICL"), col=1:2, pch=19, lwd=4)
    }
}

compute.icl <- function(dynsbm){    
    T <- ncol(dynsbm$membership)
    Q <- nrow(dynsbm$trans)
    N <- nrow(dynsbm$membership)
    pen <- 0.5*Q*log(N*(N-1)*T/2) + 0.25*Q*(Q-1)*T*log(N*(N-1)/2) # binary case
    if ("sigma" %in% names(dynsbm)) pen <- 2*pen # continuous case
    return(dynsbm$loglikelihood - ifelse(T>1,0.5*Q*(Q-1)*log(N*(T-1)),0) - pen)    
}
