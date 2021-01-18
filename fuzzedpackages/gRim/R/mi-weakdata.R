###
### weak.marginal for data
###



weak_marginal_data_list <- function(CGstats, Ad.list, Ac.list, type="ghk", details=0){

    .infoPrint(details,1,"weak_marginal_data_list: Finding weak (empirical) marginals for each generator\n")
    
    ans <- vector("list", length(Ad.list))
    for (ii in 1:length(Ad.list)){
        EE.mm <- weak_marginal_data(CGstats, disc=Ad.list[[ii]], cont=Ac.list[[ii]],
                                  type=type, details=details)
        ans[[ii]] <- EE.mm
    }
    ans
}


### cgdata: Heterogeneous CGstats
weak_marginal_data <- function(cgdata, disc=NULL, cont=NULL, type="pms", details=2){
    ## cat("weak_marginal_data - entry\n"); #str(cgdata)
    .infoPrint(details,5,
               cat("Finding weak marginal (data)",
                   "     disc:",disc, "cont:", cont, "type:", .genType(disc,cont),"\n"))
    ##str(list(dist=disc, cont=cont))
    
    ans <- switch(.genType(disc, cont),
                  "discrete"   ={.weak.datamarg.disc(cgdata, disc,       details)},
                  "mixed"      ={.weak.datamarg.mix (cgdata, disc, cont, details)},
                  "continuous" ={.weak.datamarg.cont(cgdata,       cont, details)})

    if (type == "ghk") parm_pms2ghk(ans) else ans
}

.weak.datamarg.disc <- function(cgdata, Ad.idx, details=2){
  .infoPrint(details,5, "Finding weak marginal (data-discrete):   Ad.idx: ", Ad.idx, "\n")
  
  p   <- tabMarg(cgdata$n.obs, Ad.idx)  
  res <- list(p=p / sum(p), mu=NULL, Sigma=NULL, gentype="discrete", N=cgdata[['N']], Ad.idx=Ad.idx)
  ##class(res) <- c("pms", "MIparms")
  res
}

.weak.datamarg.mix <- function(cgdata, Ad.idx, Ac.idx, details=5){
    .infoPrint(details, 5,"Finding weak marginal (data-mixed):  Ad.idx: ",
               Ad.idx, "  Ac.idx :", Ac.idx,"\n")
    #print(cgdata$n.obs)    
    ##cat(".weak.datamarg.mix\n"); print(cgdata);
    #print(.disc.levels(cgdata))

    ##n.obs <- cgdata[['n.obs']]
###n.vec <- as.numeric(cgdata[['n.obs']])
    n.vec <- c(cgdata[['n.obs']])
    n.tot <- sum(cgdata[['n.obs']])

##    print("HERE")
    p.tab <- cgdata[['n.obs']] / n.tot ## FIXME : Goes wrong here
    ## print("HERE2")
    ## print(.disc.levels(cgdata)); print(cgdata$n.obs)
        
    p.A   <- tabMarg(p.tab, Ad.idx)
    
    flevels   <- .disc.levels(cgdata)
    A.levels  <- flevels[Ad.idx]
    A.dim     <- prod(A.levels)  
    len.Ac    <- length(Ac.idx)
    n.cont    <- length(cgdata$cont.names)
    mu.A      <- cgdata$center[Ac.idx,,drop=FALSE]

    ##str(list(A.dim=A.dim, A.levels=A.levels, flevels=flevels))
    
###Sigma.A   <- cgdata$cov[as.numeric(.rowcol2idx(n.cont,Ac.idx)),,drop=FALSE]
    Sigma.A   <- cgdata$cov[c(.rowcol2idx(n.cont,Ac.idx)), , drop=FALSE]
    
### Allocate space for results
    mu.A.marg      <- matrix(0, ncol=A.dim,      nrow=len.Ac)
    Sigma.A.marg   <- matrix(0, ncol=len.Ac,     nrow=len.Ac)
    QQ             <- rep(0, len.Ac^2)  
    ia             <- rep(1, length(Ad.idx)) ## The first cell (1,1,...,1)
    ## Iteration goes 
    ##   cat(sprintf("flevels=%s, Ad.idx=%s A.dim=%d\n",
    ##               .toString(flevels), .toString(Ad.idx), A.dim))
    
##    print("HERE3"); print(cgdata$n.obs)

    ##str(list(A.dim=A.dim))
    for (ii in 1:A.dim){
##        print("HERE6"); print(.disc.levels(cgdata)); #print(cgdata$n.obs)
##        aa <- list(ia=ia, Ad.idx=Ad.idx, flevels=flevels)
##        str(aa)

        jia          <- slice2entry(ia, Ad.idx, flevels)

        ## print(jia)
        ## print("HERE7"); print(.disc.levels(cgdata)); ##print(cgdata$n.obs)
        ## counts n(ia)
        n.jia        <- n.vec[jia]
        n.ia         <- sum(n.jia)
        
        ## means \bar y (ia)
        mu.jia <- mu.A[,jia, drop=FALSE]
        mu.ia    <- rowSumsPrim(.colmult(n.jia/n.ia, mu.jia))
        mu.A.marg[,ii] <- mu.ia
        
        ## SSD(ia) 
        S.jia  <- Sigma.A[,jia,drop=FALSE]
        vvv1   <- rowSumsPrim(.colmult(n.jia , S.jia))
        sum.ssd.j  <- matrix(vvv1, nrow=length(Ac.idx)) ## sum_{j:ja=ia} SSD(j)
        
        mu.dif  <- mu.jia - mu.ia
        quad    <- .vMMt(n.jia, mu.dif)
        QQ      <- QQ + (sum.ssd.j + quad)/n.tot
        ia <- next_cell(ia, A.levels) ## FIXME Think that A.levels and flevels are always same, but do check!
    }
    
    rownames(mu.A.marg) <- rownames(Sigma.A.marg)
    ans <- list(p=p.A, mu=mu.A.marg, Sigma=QQ, 
                gentype="mixed", N=cgdata[['N']],
                Ad.idx=Ad.idx, Ac.idx=Ac.idx)
    ##class(ans) <- c("pms", "MIparms")


    ## print("HERE4");  print(cgdata$n.obs)    
    ans
}

.weak.datamarg.cont <- function(cgdata, Ac.idx, details=2){
  .infoPrint(details, 5, "Finding weak marginal (data-cont) :  Ac.idx: ", Ac.idx,"\n")

  n.cont   <- length(cgdata$cont.names)
  n.i      <- as.numeric(cgdata$n.obs)
  N        <- sum(n.i)
  
  mu.i     <- cgdata$center[Ac.idx,,drop=FALSE]
  Sigma.i  <- cgdata$cov[as.numeric(.rowcol2idx(n.cont,Ac.idx)),, drop=FALSE]
  p.i      <- n.i / sum(n.i)

  ssdA.i   <- .colmult(n.i, Sigma.i)
  sum.ssdA <- matrix(rowSumsPrim(ssdA.i), nrow=length(Ac.idx))

  mu.A     <- rowSumsPrim(.colmult(p.i, mu.i))
  mu.dif   <- mu.i - mu.A
  quad     <- .vMMt(n.i, mu.dif)
  Sigma.A  <- (sum.ssdA + quad) / N

  dim(mu.A) <- c(dim(mu.i)[1],1)
  rownames(Sigma.A) <- colnames(Sigma.A) <- rownames(mu.A) <- rownames(mu.i) 
  
  ans <- list(p=1, mu=mu.A, Sigma=Sigma.A, gentype="continuous", N=cgdata[['N']],  Ac.idx = Ac.idx)
  ##class(ans) <- c("pms", "MIparms")
  ans             
}





