
GenomicRel = function(M){
#       markers must be 0, 1, 2 for homozygous, heterozygous and other homozygous

    M= M+1  ## since M is -1,0,1
    p1=round((apply(M,2,sum)+nrow(M))/(nrow(M)*2),3)
    p=2*(p1-.5)
    P = matrix(p,byrow=T,nrow=nrow(M),ncol=ncol(M))
    Z = as.matrix(M-P)

    b=1-p1
    c=p1*b
    d=2*(sum(c))

    ZZt = Z %*% t(Z)
    G = (ZZt/d)
    #invG=solve(G)

    return(G)
}



