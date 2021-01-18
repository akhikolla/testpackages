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


alluvial.plot <- function(dynsbm, timestep.abbrev="T", minimal.flow=1, only.present=FALSE){
                                        #library(riverplot)
    if (! "dynsbm" %in% class(dynsbm)) stop("Object must be of class 'dynsbm'")
    T <- ncol(dynsbm$membership)
    membership <- dynsbm$membership
    if(only.present) membership[membership==0] <- NA
    edges <- c()
    for (t in 1:(T-1)){
        groups1 <- membership[,t]
        groups2 <- membership[,t+1]
        tb <- table(groups1,groups2)
        groupname1 <- dimnames(tb)$groups1
        groupname2 <- dimnames(tb)$groups2
        for (q in 1:dim(tb)[1])
            for (l in 1:dim(tb)[2]){
                if (tb[q,l]>=minimal.flow){
                    if(groupname1[q]=="0" && groupname2[l]=="0"){ # SKIP 0->0
                    } else{
                        edges <- rbind(edges, c(paste(timestep.abbrev,t,"-",groupname1[q],sep=""), paste(timestep.abbrev,t+1,"-",groupname2[l],sep=""), tb[q,l]))
                    }
                }
            }
    }
    edges <- data.frame(N1 = as.character(edges[,1]),
			N2 = as.character(edges[,2]),
			Value=as.numeric(edges[,3]),
			stringsAsFactors = F)
    nodes <- data.frame(ID = unique(c(edges$N1, edges$N2)), stringsAsFactors = F)
    vID <- strsplit(nodes$ID,"-")
    nodes$x <- as.numeric(substr(sapply(vID, FUN=function(v){return(v[1])}),2,5))
    nodes$y <- as.numeric(sapply(vID, FUN=function(v){return(v[2])}))
    rownames(nodes) <- nodes$ID
    
    styles <- lapply(nodes$y, function(yy) {
        list(col = ifelse(yy==0,colors()[1],colors()[200]), lty = 0, textcol = "black")
    })
    names(styles) <- nodes$ID
    rp <- list(nodes = nodes, edges = edges, styles = styles)
    class(rp) <- c(class(rp), "riverplot")
    plot(rp, nodewidth = .5, plot_area = 0.95)
}

gradient.polygon <- function(x1,x2,y1,y2,col1,col2,ylim.low){
    gap <- (x2-x1)/20
    a <- (y2-y1)/(x2-x1)
    b <- (x2*y1-x1*y2)/(x2-x1)
    xx <- x1
    colfunc <- colorRampPalette(c(col1, col2))
    my.col <- colfunc(20)
    for(i in 1:20){
        xx <- x1+(i-1)*gap
        yy1 <- a*xx+b
        yy2 <- a*(xx+gap)+b
        polygon(c(xx,xx,xx+gap,xx+gap), c(ylim.low,yy1,yy2,ylim.low),                 
                col=my.col[i],
                border=NA)
    }
}

connectivity.plot <- function(dynsbm, Y){
    if (! "dynsbm" %in% class(dynsbm)) stop("Object must be of class 'dynsbm'")
                                        #library(RColorBrewer)
    directed <- dynsbm$directed
    T <- ncol(dynsbm$membership)
    Q <- nrow(dynsbm$trans)
    membership <- dynsbm$membership
    old.mar <- par()$mar
    par(mar = rep(0.33, 4))
    par(mfrow=c(Q,Q))
    ## BETA
    beta <- array(0, c(T,Q,Q))
    for (q in 1:Q){
        for (l in 1:Q){
            beta[,q,l] <- sapply(1:T, FUN=function(t){
                subYt <- Y[t, which(membership[,t]==q),which(membership[,t]==l)]
                if("matrix" %in% class(subYt)){
                    d1 <- nrow(subYt)
                    d2 <- ncol(subYt)
                    if(q==l){ #subYt is square
                        if(dynsbm$self.loop){
                            this.beta <- sum(subYt>0) / (d1*d1)
                        } else{
                            diag(subYt) <- NA
                            this.beta <- sum(subYt>0, na.rm=T) / (d1*(d1-1))
                        }
                    } else{
                        this.beta <- sum(subYt>0) / (d1*d2)
                    }
                } else{ ## any(class(subYt) %in% c("numeric","integer"))
                    this.beta <- sum(subYt>0) / length(subYt)
                }
                return(this.beta)
            })
        }
    }
    beta[is.na(beta)] <- 0 
    if(!("gamma" %in% names(dynsbm)) && !("mu" %in% names(dynsbm))){
        ## binary case  
        for (q in 1:Q){
            for (l in 1:Q){
                y <- beta[,q,l]
                par(mgp=c(0,-1.4,0)) 
                plot(1:T, y, ylim=c(0,1), type="n", xlab="",ylab="", xaxt="n", yaxt="n")
                polygon(c(1,1:T,T), c(0,y,0), col=brewer.pal(5,"Blues")[2])
                lines(1:T, y, type="b", pch = 20)  
                axis(1,at=1:T, tck = 0.02, labels=rep("",T))
                axis(2,at=c(0,1), labels=c("0","1"),tck = 0.02)         
            }
        }
    }
    if("gamma" %in% names(dynsbm)){
        ## discrete case
        K <- dim(dynsbm$gamma)[4]
        for (q in 1:Q){
            for (l in 1:Q){
                y <- beta[,q,l]
                par(mgp=c(0,-1.4,0)) 
                plot(1:T, y, ylim=c(0,1), type="n", xlab="",ylab="", xaxt="n", yaxt="n")
                polygon(c(1,1:T,T), c(0,y,0), col=brewer.pal(K,"Blues")[K])
                z <- sapply(1:T, function(t){
                    subYt <- Y[t, which(membership[,t]==q),which(membership[,t]==l)]
                    if (length(subYt)==0){
                                        #subYt is empty
                        this.t <- rep(0,K+1)
                    } else{
                        if(("matrix" %in% class(subYt)) & q==l){ #subYt is square
                            if(!dynsbm$self.loop) diag(subYt) <- NA
                        } else{ ## ("integer" %in% class(subYt)==) #subYt is a vector
                            if(!dynsbm$self.loop) subYt <- 0L
                        }
                        this.t <- (table(c(subYt,0:K))-1)/length(subYt[!is.na(subYt)])
                        this.t <- this.t[order(as.integer(names(this.t)))][c(2:(K+1),1)]
                    }
                    cumsum(this.t)
                })
                for(k in ((K-1):1)){
                    polygon(c(1,1:T,T), c(0,z[k,],0), col=brewer.pal(K,"Blues")[k])
                    lines(1:T, z[k,], type="b", pch = 20)
                }
                lines(1:T, y, type="b", pch = 20)
                axis(1,at=1:T, tck = 0.02, labels=rep("",T))
                axis(2,at=c(0,1), labels=c("0","1"),tck = 0.02)
            }
        }
    }
    if("mu" %in% names(dynsbm)){
        connec <- array(0, c(T,Q,Q))
        for (q in 1:Q){
            for (l in 1:Q){
                connec[,q,l] <- sapply(1:T, FUN=function(t){
                    subYt <- Y[t, which(membership[,t]==q),which(membership[,t]==l)]
                    if(("matrix" %in% class(subYt)) & q==l){ #subYt is square
                        if(!dynsbm$self.loop) diag(subYt) <- NA
                    } else{ # ("integer" %in% class(subYt)==) #subYt is a vector
                        if(!dynsbm$self.loop) subYt <- 0
                    }
                    if (length(subYt[subYt>0])){
                        mean(subYt[subYt>0],na.rm=T)
                    } else 0         
                })
            }
        }
        for (q in 1:Q){
            for (l in 1:Q){
                y <- beta[,q,l]
                par(mgp=c(0,-1.4,0)) 
                plot(1:T, y, ylim=c(0,1), type="n", xlab="",ylab="", xaxt="n", yaxt="n")
                for(t in 1:(T-1)){
                    col1 <- (colorRampPalette(brewer.pal(9,"Blues"))(12))[3:15][1+floor((connec[t,q,l]-min(connec))/(max(connec)-min(connec))*9)]
                    col2 <- (colorRampPalette(brewer.pal(9,"Blues"))(12))[3:15][1+floor((connec[t+1,q,l]-min(connec))/(max(connec)-min(connec))*9)]
                    col1 <- brewer.pal(9,"Blues")[1+floor((connec[t,q,l]-min(connec))/(max(connec)-min(connec))*8)]
                    col2 <- brewer.pal(9,"Blues")[1+floor((connec[t+1,q,l]-min(connec))/(max(connec)-min(connec))*8)]
                    gradient.polygon(t,t+1,y[t],y[t+1],col1,col2,0)
                }   
                lines(1:T, y, type="b", pch = 20)
                axis(1,at=1:T, tck = 0.02, labels=rep("",T))
                axis(2,at=c(0,1), labels=c("0","1"),tck = 0.02)
            }
        }
    }
    par(mar=old.mar)
    par(mfrow=c(1,1))
}

adjacency.plot <- function(dynsbm, Y, present=NULL, col=heat.colors(9)){
    if (is.null(present)){
        present <- matrix(0L,dim(Y)[2],dim(Y)[1])    
        for (t in 1:dim(Y)[1])
            present[union(which(apply(Y[t,,],1,FUN=function(v) sum(v>0))>0),which(apply(Y[t,,],2,FUN=function(v) sum(v>0))>0)),t] <- 1    		
    }
    ## Alternative: col <- gray.colors(10, start=0.8, end=0.)
    if(!("gamma" %in% names(dynsbm)) && !("mu" %in% names(dynsbm))){
        col=col[1]
    }
    par(mfrow=c(min(dim(Y)[1],3), ceiling(dim(Y)[1]/(min(dim(Y)[1],3)))))
    par(mar = rep(0.5,4))
    for(t in 1:dim(Y)[1]){
        Yt <- Y[t,,]
        index <- which(dynsbm$membership[,t]>0)
        subYt <- Yt[index,index]
        if(is.matrix(subYt)){
            membershipt <- dynsbm$membership[,t]
            membershipt <- membershipt[which(membershipt>0)]
            image(t(subYt[rev(order(membershipt)), order(membershipt)]), xaxt = "n", xlab = "", yaxt = "n", ylab = "",
                  col=c("white",col))            
            nb <- cumsum(rev(table(membershipt)))
            abline(h=(-1/nb[length(nb)]/2+(1+1/nb[length(nb)])*nb[-length(nb)]/nrow(subYt)), col=2)
            nb <- cumsum(table(membershipt))
            abline(v=(-1/nb[length(nb)]/2+(1+1/nb[length(nb)])*nb[-length(nb)]/nrow(subYt)), col=2)
        } else{
            image(matrix(1))
        }
    }
}
