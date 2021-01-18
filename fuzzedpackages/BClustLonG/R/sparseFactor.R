#' A Dirichlet process mixture model for clustering longitudinal gene expression data.
#' @param data Data list with three elements: Y (gene expression data with
#' each column being one gene), ID, and years.
#' (The names of the elements have to be matached exactly.
#' See the data in the example section more info)
#' @param iter Number of iterations (excluding the thinning).
#' @param thin Number of thinnings.
#' @param savePara Logical variable indicating if all the parameters needed to be saved.
#' Default value is FALSE, in which case only the membership indicators are saved.
#' @param infoVar Either "both" (using both intercepts and slopes for clustering)
#' or "int" (using only intercepts for clustering)
#' @param factor Logical variable indicating whether factor analysis model is wanted.
#' @param hyperPara A list of hyperparameters with default values.
#' @return returns a list with following objects.
#' \item{e.mat}{Membership indicators from all iterations.}
#' \item{All other parameters}{only returned when savePara=TRUE.}
#' @references Jiehuan Sun, Jose D. Herazo-Maya, Naftali Kaminski, Hongyu Zhao,
#' and Joshua L. Warren. "A Dirichlet process mixture model for clustering
#' longitudinal gene expression data." Statistics in Medicine 36, No. 22 (2017): 3495-3506.
#' @examples
#' data(data)
#' ## increase the number of iterations
#' ## to ensure convergence of the algorithm
#'
#' res = BClustLonG(data, iter=20, thin=2,savePara=FALSE,
#' infoVar="both",factor=TRUE)
#' ## discard the first 10 burn-ins in the e.mat
#' ## and calculate similarity matrix
#' ## the number of burn-ins has be chosen s.t. the algorithm is converged.
#' mat = calSim(t(res$e.mat[,11:20]))
#' clust = maxpear(mat)$cl ## the clustering results.
#'\dontrun{
#' ## if only want to include intercepts for clustering
#' ## set infoVar="int"
#' res = BClustLonG(data, iter=10, thin=2,savePara=FALSE,
#' infoVar="int",factor=TRUE)
#'
#' ## if no factor analysis model is wanted
#' ## set factor=FALSE
#' res = BClustLonG(data, iter=10, thin=2,savePara=FALSE,
#' infoVar="int",factor=TRUE)
#' }
BClustLonG <- function(data=NULL, iter=20000, thin = 2, savePara=FALSE,
                       infoVar=c("both","int")[1], factor=TRUE,
                       hyperPara=list(v1=0.1, v2=0.1, v=1.5, c=1, a=0, b=10,
                                      cd=1, aa1=2, aa2=1, alpha0= -1,
                                      alpha1 = -1e-4, cutoff = 1e-4, h =100)
){
    if(infoVar=="both" & factor){
        cat("clustering on both intercepts and slopes...\n")
        cat("use factor analysis model to approximate the covariance matrix...\n")
        res = sparseBoth(data, iter, thin , savePara, hyperPara)
    }else if(infoVar=="int" & factor){
        cat("clustering only on intercepts... \n")
        cat("use factor analysis model to approximate the covariance matrix...\n")
        res = sparseInt(data, iter, thin , savePara, hyperPara)
    }else if(infoVar=="both" & (!factor) ){
        cat("clustering on both intercepts and slopes... \n")
        cat("assume diagonal covariance matrix...\n")
        res = sparseBothNof(data, iter, thin , savePara, hyperPara)
    }else if (infoVar=="int" & (!factor)){
        cat("clustering only on intercepts... \n")
        cat("assume diagonal covariance matrix...\n")
        res = sparseIntNof(data, iter, thin , savePara, hyperPara)
    }else{
        stop("At least one of the variables 'infoVar' or
             'factor' are not correctly specified!" )
    }

    res
    }




#' Internal function to sample c
#' @noRd
samCv2 <- function(a,b,e,c,cd,n){
    d <- length(unique(e))
    v <- log((c-a)/(b-c))
    vp <- rnorm(1,v,cd)
    cp <- (a+exp(vp)*b)/(exp(vp)+1)
    pp <- d*log(cp)+lgamma(cp) - lgamma(cp+n) + log(b-a) + vp - 2*log(exp(vp)+1) -
        ( d*log(c)+lgamma(c) - lgamma(c+n) + log(b-a) + v - 2*log(exp(v)+1) )
    pp <- as.numeric(runif(1) <= exp(pp))
    return(c*(pp==0)+cp*(pp==1))
}

#' Internal function to sample A
#' @noRd
samA <- function(Y,years,sig,B,muA,ID,e,sigmaInv){

    Y <- Y - B[ID,]*years
    sigma <- diag(1/sig)
    A <- matrix(0,nrow=max(ID),ncol=length(sig))
    tmp <- muA[,e]
    for(i in 1:max(ID)){
        tmp.sig <- solve(sigmaInv + sum(ID==i)*sigma)
        tmp.mu <- tmp.sig%*%(sigmaInv%*%(tmp[,i]) + sigma%*%apply(Y*(ID==i),2,sum) )
        A[i,] <- mvrnorm(1,tmp.mu,tmp.sig)
    }
    return(A)
}

#' Internal function to sample mean of A
#' @noRd
samMuA <- function(e,A,muA0,sigma0Inv,sigmaInv){

    K <- max(e)
    mat <- matrix(0,nrow=nrow(sigma0Inv),ncol=K)
    M <- sapply(seq_len(K),function(k){as.numeric(e==k)})


    tmp.sig <- lapply(seq_len(K), function(x){
        solve(sigma0Inv + sum(M[,x])*sigmaInv)
    })
    tmp.mu <- lapply(seq_len(K),function(k){
        tmp.sig[[k]]%*%(sigmaInv%*%t(A)%*%M[,k] + sigma0Inv%*%muA0) })
    mat[,seq_len(K)] <- sapply(seq_len(K),function(k){
        mvrnorm(1,tmp.mu[[k]],tmp.sig[[k]])})
    mat
}

#' Internal function to sample hyperparameters for mean of A
#' @noRd
samMuA0 <- function(sigma0Inv,muA,h){

    I <- diag(rep(1,ncol(sigma0Inv)))/h
    tmp.sig <- solve(I + ncol(muA)*sigma0Inv)
    tmp.mu <- tmp.sig%*%sigma0Inv%*%(apply(muA,1,sum))

    return(mvrnorm(1,tmp.mu,tmp.sig))
}

#' Internal function to sample eta in factor analysis
#' @noRd
samEtaA <- function(siga,A,lamA,muA,e,sig0a){

    I <- diag(rep(1,ncol(lamA)))/sig0a
    sigma <- diag(1/siga)
    A <- A - t(muA[,e])
    tmp.sig <- solve(I + t(lamA)%*%sigma%*%lamA)
    tmp.mu <- tmp.sig%*%t(lamA)%*%sigma%*%t(A)

    return( t(mvrnorm(nrow(A),rep(0,ncol(lamA)),tmp.sig)) + tmp.mu )
}


#' Internal function to sample diagonal elements of covariance matrix
#' @noRd
samSigA <- function(v1,v2,muA,e,A,lamA,etaA){
    v1 <- v1 + nrow(A)/2
    A <- A - t(muA[,e])
    A <- A  - t(lamA%*%etaA)
    v2 <- v2 + apply(A^2,2,sum)/2
    return(1/rgamma(length(v2),v1,v2))
}


#' Internal function to sample B
#' @noRd
samB <- function(Y,years,sig,sigmaInv,A,muB,e,ID){

    Y <- Y - A[ID,]
    sigma <- diag(1/sig)
    B <- matrix(0,nrow=max(ID),ncol=length(sig))
    tmp <- muB[,e]
    for(i in 1:max(ID)){
        tmp.sig <- solve(sigmaInv + sum(years[ID==i]^2)*sigma)
        tmp.mu <- tmp.sig%*%(sigmaInv%*%tmp[,i] + sigma%*%apply(Y*years*(ID==i),2,sum) )
        B[i,] <- mvrnorm(1,tmp.mu,tmp.sig)
    }
    return(B)
}

#' Internal function to sample mean of B
#' @noRd
samMuB <- function(e,B,muB0,sigma0Inv,sigmaInv){

    K <- max(e)
    mat <- matrix(0,nrow=nrow(sigma0Inv),ncol=K)
    M <- sapply(seq_len(K),function(k){as.numeric(e==k)})


    tmp.sig <- lapply(seq_len(K), function(x){
        solve(sigma0Inv + sum(M[,x])*sigmaInv)
    })
    tmp.mu <- lapply(seq_len(K),function(k){
        tmp.sig[[k]]%*%(sigmaInv%*%t(B)%*%M[,k] + sigma0Inv%*%muB0) })
    mat[,seq_len(K)] <- sapply(seq_len(K),function(k){
        mvrnorm(1,tmp.mu[[k]],tmp.sig[[k]])})
    mat
}


#' Internal function to sample hyperparameters for mean of B
#' @noRd
samMuB0 <- function(sigma0Inv,muB,h){

    I <- diag(rep(1,ncol(sigma0Inv)))/h
    tmp.sig <- solve(I + ncol(muB)*sigma0Inv)
    tmp.mu <- tmp.sig%*%sigma0Inv%*%(apply(muB,1,sum))

    return(mvrnorm(1,tmp.mu,tmp.sig))
}

#' Internal function to sample diagonal elements of covariance matrix
#' @noRd
samSigB <- function(v1,v2,muB,e,B,lamB,etaB){
    v1 <- v1 + nrow(B)/2
    B <- B - t(muB[,e])
    B <- B  - t(lamB%*%etaB)
    v2 <- v2 + apply(B^2,2,sum)/2
    return(1/rgamma(length(v2),v1,v2))
}

#' Internal function to sample eta in the factor analysis
#' @noRd
samEtaB <- function(sigb,B,lamB,muB,e,sig0b){

    I <- diag(rep(1,ncol(lamB)))/sig0b
    sigma <- diag(1/sigb)
    B <- B - t(muB[,e])
    tmp.sig <- solve(I + t(lamB)%*%sigma%*%lamB)
    tmp.mu <- tmp.sig%*%t(lamB)%*%sigma%*%t(B)

    return( t(mvrnorm(nrow(B),rep(0,ncol(lamB)),tmp.sig)) + tmp.mu )
}

#' Internal function to sample B for sparse int
#' @noRd
samBInt <- function(Y,years,sig,sigb,A,b0,ID){

    Y <- Y - A[ID,]
    sigma <- diag(1/sig)
    sigmab <- diag(1/sigb)
    B <- matrix(0,nrow=max(ID),ncol=length(sig))
    for(i in 1:max(ID)){
        tmp.sig <- solve(sigmab + sum(years[ID==i]^2)*sigma)
        tmp.mu <- tmp.sig%*%(sigmab%*%b0 + sigma%*%apply(Y*years*(ID==i),2,sum) )
        B[i,] <- mvrnorm(1,tmp.mu,tmp.sig)
    }
    return(B)
}

#' Internal function to sample mean of B for sparse int
#' @noRd
SamB0Int <- function(sigb,h,B){
    G <- length(sigb)
    sigma <- diag(1/sigb)
    I <- diag(rep(1/h,G))
    tmp.sig <- solve(I + nrow(B)*sigma)
    tmp.mu <- tmp.sig%*%sigma%*%apply(B,2,sum)
    mvrnorm(1,tmp.mu,tmp.sig)
}

#' Internal function to sample diagonal elements of covariance matrix for sparse int
#' @noRd
samSigBInt <- function(v1,v2,b0,B){
    v1 <- v1 + nrow(B)/2
    B <- B - matrix(rep(b0,nrow(B)), nrow=nrow(B),byrow=T)
    v2 <- v2 + apply(B^2,2,sum)/2
    return(1/rgamma(length(v2),v1,v2))
}


#' Internal function to sample gene specific variances
#' @noRd
samSig <- function(v1,v2,Y,years,ID,B,A){
    v1 <- v1 + nrow(Y)/2
    Y <- Y - A[ID,]
    Y <- Y - B[ID,]*years
    v2 <- v2 + apply(Y^2,2,sum)/2
    return(1/rgamma(length(v2),v1,v2))
}

#' Internal function to sample hyperparameters for mean of A
#' @noRd
samSig0 <- function(v1,v2,muA0,muA){
    v1 <- v1 + nrow(muA)*ncol(muA)/2
    res <- muA - muA0
    v2 <- v2 + sum(res^2)/2
    return(1/rgamma(1,v1,v2))
}

#' Internal function to sample hyperparameters in factor analysis
#' @noRd
samaV2 <- function(a,delta,vv=1,vv2=2,vv1=1){

    loglik <- sum( dgamma(delta,a,1,log=T) ) + dgamma(a,vv2,vv1,log=T)

    v <- log(a)
    vp <- rnorm(1,v,vv)
    ap <- exp(vp)

    loglikp <- sum( dgamma(delta,ap,1,log=T) ) + dgamma(ap,vv2,vv1,log=T)

    loglik <- sum(loglik)  + v
    loglikp <- sum(loglikp) + vp

    dif <- exp(loglikp - loglik)
    if(runif(1) < dif){a <- ap}
    return(a)
}

#' Internal function to sample hyperparameters in factor analysis
#' @noRd
samPhi <- function(v1,v2,lamA,tau){
    v1 = v1 + 1/2
    v2 = v2 + as.vector(matrix(t(lamA)^2,nrow=ncol(lamA))*tau)/2
    t(matrix(rgamma(length(v2),v1,v2),nrow=ncol(lamA)))
}

#' Internal function to sample hyperparameters in factor analysis
#' @noRd
samDelta <- function(a1,a2,phi,lamA,tau,delta){
    tmp <- apply(phi*(lamA^2),2,sum)
    p1 <- sapply(1:ncol(lamA),function(i){
        nrow(lamA)*(ncol(lamA)-i+1)/2
    })
    p2 <- sapply(1:ncol(lamA),function(i){
        sum(tau[i:ncol(lamA)]*tmp[i:ncol(lamA)]/2/delta[i])
    })
    p1 <- p1 + c(a1,rep(a2,length(p1)-1))
    p2 <- p2 + 1
    rgamma(length(p2),p1,p2)
}




#' Internal function to perform clustering using both intercepts and slopes
#' and using factor analysis model.
#' @noRd
sparseBoth <- function(data=NULL, iter=100, thin = 2, savePara=FALSE,
                       hyperPara=list(v1=0.1, v2=0.1, v=1.5, c=1, a=0, b=10,
                                      cd=1, aa1=2, aa2=1, alpha0= -1,
                                      alpha1 = -1e-4, cutoff = 1e-4, h =100)
){
    ## data section
    Y = data$Y
    ID = data$ID
    ID = match(ID,ID[!duplicated(ID)])
    years = data$years

    Y <- scale(Y,scale=TRUE)
    G <- ncol(Y)
    n <- max(ID)
    X <- cbind(1,years) ## years has unit scale

    ## hyperparameters
    v1 = hyperPara$v1; v2 = hyperPara$v2; v = hyperPara$v
    c = hyperPara$c; a = hyperPara$a; b = hyperPara$b
    cd = hyperPara$cd; aa1=2; aa2=1
    alpha0 = hyperPara$alpha0; alpha1 = hyperPara$alpha1
    cutoff = hyperPara$cutoff; h = hyperPara$h

    ## initiations

    a.mat <- matrix(0,nrow=n,ncol=G)
    b.mat <- matrix(0,nrow=n,ncol=G)
    sig <- rep(1,G)
    for(g in 1:G){
        #print(i)
        fit <- tryCatch(lmer(Y[,g]~(1+years|ID)))
        a.mat[,g] <- coef(fit)$ID[,2]
        b.mat[,g] <- coef(fit)$ID[,1]
        sig[g] <- summary(fit)$sigma^2
    }

    siga <- apply(a.mat,2,var)
    sigb <- apply(b.mat,2,var)

    sig[sig<1e-4] <- 0.01
    siga[siga<1e-4] <- 0.01
    sigb[sigb<1e-4] <- 0.01

    M <- 10

    e <- sample(1:as.integer(n/5),n,replace = TRUE)

    c1.a = 10
    c2.a = 10
    delta.a <- c(rgamma(1,c1.a,1),rgamma(M-1,c2.a,1))
    tau.a <- cumprod(delta.a)
    phi.a <- matrix(rgamma(G*M,v,v),ncol=M)

    lamA <- sapply(1:M,function(i){
        rnorm(G,0,1/sqrt(phi.a[,i]*tau.a[i]) )
    })
    etaA <- matrix(rnorm(n*M),nrow=M)

    c1.b = 10
    c2.b = 10
    delta.b <- c(rgamma(1,c1.b,1),rgamma(M-1,c2.b,1))
    tau.b <- cumprod(delta.b)
    phi.b <- matrix(rgamma(G*M,v,v),ncol=M)

    lamB <- sapply(1:M,function(i){
        rnorm(G,0,1/sqrt(phi.b[,i]*tau.b[i]) )
    })
    etaB <- matrix(rnorm(n*M),nrow=M)

    B <- b.mat
    A <- a.mat

    muB0 <- apply(b.mat,2,mean)
    muA0 <- apply(a.mat,2,mean)

    sig0b = 1
    sig0a = 1

    sigmaA0 <- diag(rep(sig0a,G))
    sigmaA0Inv <- diag(rep(1/sig0a,G))
    sigmaB0 <- diag(rep(sig0b,G))
    sigmaB0Inv <- diag(rep(1/sig0b,G))

    sigmaB = lamB%*%t(lamB) + diag(sigb)
    sigmaBInv <- chol2inv(chol(sigmaB))

    sigmaA = lamA%*%t(lamA) + diag(siga)
    sigmaAInv <- chol2inv(chol(sigmaA))

    pt.a = exp(alpha0+alpha1*(1:iter))
    ut.a <- runif(iter)
    pt.b = exp(alpha0+alpha1*(1:iter))
    ut.b <- runif(iter)

    ## construct data to store results
    e.mat <- matrix(0,n,iter)
    if(savePara){
        c.vec <- rep(0,iter)
        sigmaA.vec <- matrix(0,G,iter)
        sigmaB.vec <- matrix(0,G,iter)

        sig0b.vec <- rep(0,iter)
        sig0a.vec <- rep(0,iter)
        sig.vec <- matrix(0,G,iter)

        muB.vec <- matrix(0,G*n,iter)
        muB0.vec <- matrix(0,G,iter)
        muA.vec <- matrix(0,G*n,iter)
        muA0.vec <- matrix(0,G,iter)

        A.vec <- matrix(0,n*G,iter)
        B.vec <- matrix(0,n*G,iter)
    }

    ### running MCMC
    for( i in 1:iter){
        print(paste("iterations ",i))

        for(thinning in 1:thin){

            e <- polyurncppBoth(e,A, muA0, sigmaA,  sigmaAInv, sigmaA0,sigmaA0Inv,
                                B, muB0, sigmaB, sigmaBInv, sigmaB0,sigmaB0Inv, c)
            e <- match(e,sort(unique(e)))
            c <- samCv2(a,b,e,c,cd,n)

            ## DP for B
            muB <- samMuB(e,B,muB0,sigmaB0Inv,sigmaBInv)
            muB0 <- samMuB0(sigmaB0Inv,muB,h)

            ## DP for A
            muA <- samMuA(e,A,muA0,sigmaA0Inv,sigmaAInv)
            muA0 <- samMuA0(sigmaA0Inv,muA,h)

            ## sampling lamA,etaA
            lamA <- samLamV2Cpp(A - t(muA[,e]),  etaA, siga, lamA, phi.a, tau.a)

            if(ut.a[i]<=pt.a[i]){
                ind <- apply(lamA,2,function(x){all(abs(x)<cutoff)})
                if(sum(ind)>0){
                    lamA <- lamA[,!ind]
                    phi.a <- phi.a[,!ind]
                    delta.a <- delta.a[!ind]
                    tau.a <- cumprod(delta.a)
                }else{
                    delta.a <- c(delta.a,rgamma(1,c2.a,1))
                    tau.a <- cumprod(delta.a)
                    phi.a <- cbind(phi.a,rgamma(G,v,v))
                    lamA <- cbind(lamA,rnorm(G,0,1/sqrt(phi.a[,ncol(phi.a)]*tau.a[ncol(phi.a)]) ) )
                }
            }
            etaA <- samEtaA(siga,A,lamA,muA,e,1)
            etaA <- t(scale(t(etaA)))
            phi.a <- samPhi(v,v,lamA,tau.a)
            delta.a <- samDelta(c1.a,c2.a,phi.a,lamA,tau.a,delta.a)
            tau.a <- cumprod(delta.a)
            c1.a <- samaV2(c1.a,delta.a[1], vv=1,vv2=aa1,vv1=aa2)
            c2.a <- samaV2(c2.a,delta.a[-1], vv=1,vv2=aa1,vv1=aa2)

            #### lamB and etaB
            lamB <- samLamV2Cpp(B - t(muB[,e]),  etaB, sigb, lamB, phi.b, tau.b)

            if(ut.b[i]<=pt.b[i]){
                ind <- apply(lamB,2,function(x){all(abs(x)<cutoff)})
                if(sum(ind)>0){
                    lamB <- lamB[,!ind]
                    phi.b <- phi.b[,!ind]
                    delta.b <- delta.b[!ind]
                    tau.b <- cumprod(delta.b)
                }else{
                    delta.b <- c(delta.b,rgamma(1,c2.b,1))
                    tau.b <- cumprod(delta.b)
                    phi.b <- cbind(phi.b,rgamma(G,v,v))
                    lamB <- cbind(lamB,rnorm(G,0,1/sqrt(phi.b[,ncol(phi.b)]*tau.b[ncol(phi.b)]) ) )
                }
            }
            etaB <- samEtaB(sigb,B,lamB,muB,e,1)
            etaB <- t(scale(t(etaB)))
            phi.b <- samPhi(v,v,lamB,tau.b)
            delta.b <- samDelta(c1.b,c2.b,phi.b,lamB,tau.b,delta.b)
            tau.b <- cumprod(delta.b)
            c1.b <- samaV2(c1.b,delta.b[1], vv=1,vv2=aa1,vv1=aa2)
            c2.b <- samaV2(c2.b,delta.b[-1], vv=1,vv2=aa1,vv1=aa2)

            ### sampling A
            A <- samA(Y,years,sig,B,muA,ID,e,sigmaAInv)
            siga <- samSigA(v1,v2,muA,e,A,lamA,etaA)
            sig0a <- samSig0(v1,v2[1],muA0,muA)

            ### sampling B
            B <- samB(Y,years,sig,sigmaBInv,A,muB,e,ID)
            sigb <- samSigB(v1,v2,muB,e,B,lamB,etaB)
            sig0b <- samSig0(v1,v2[1],muB0,muB)

            sig <- samSig(v1,v2,Y,years,ID,B,A)

            sigmaB = lamB%*%t(lamB) + diag(sigb)
            sigmaBInv <- chol2inv(chol(sigmaB))

            sigmaA = lamA%*%t(lamA) + diag(siga)
            sigmaAInv <- chol2inv(chol(sigmaA))

            sigmaA0 <- diag(rep(sig0a,G))
            sigmaA0Inv <- diag(rep(1/sig0a,G))

            sigmaB0 <- diag(rep(sig0b,G))
            sigmaB0Inv <- diag(rep(1/sig0b,G))
        }

        print(e)
        if(savePara){
            c.vec[i] <- c

            muB.vec[,i] <- as.vector(muB[,e])
            muB0.vec[,i] <- muB0
            muA.vec[,i] <- as.vector(muA[,e])
            muA0.vec[,i] <- muA0

            sig0b.vec[i] <- sig0b
            sig0a.vec[i] <- sig0a
            sigmaA.vec[,i] <- diag(sigmaA)
            sigmaB.vec[,i] <- diag(sigmaB)

            sig.vec[,i] <- sig
            A.vec[,i] <- as.vector(A)
            B.vec[,i] <- as.vector(B)
        }
        e.mat[,i] <- e
    }

    ## return results
    if(savePara){
        list(e.mat=e.mat,muA.vec=muA.vec, muA0.vec=muA0.vec,sig0a.vec=sig0a.vec,
             muB.vec=muB.vec, muB0.vec=muB0.vec,sig0b.vec=sig0b.vec,sig.vec=sig.vec,
             A.vec=A.vec,B.vec=B.vec,sigmaA.vec=sigmaA.vec,sigmaB.vec=sigmaB.vec,
             c.vec=c.vec)
    }else{
        list(e.mat=e.mat)
    }

}

#' Internal function to perform clustering using only intercepts
#' and using factor analysis model.
#' @noRd
sparseInt <- function(data=NULL, iter=100, thin = 2, savePara=FALSE,
                      hyperPara=list(v1=0.1, v2=0.1, v=1.5, c=1, a=0, b=10,
                                     cd=1, aa1=2, aa2=1, alpha0= -1,
                                     alpha1 = -1e-4, cutoff = 1e-4, h =100)
){
    ## data section
    Y = data$Y
    ID = data$ID
    ID = match(ID,ID[!duplicated(ID)])
    years = data$years

    Y <- scale(Y,scale=TRUE)
    G <- ncol(Y)
    n <- max(ID)
    X <- cbind(1,years) ## years has unit scale

    ## hyperparameters
    v1 = hyperPara$v1; v2 = hyperPara$v2; v = hyperPara$v
    c = hyperPara$c; a = hyperPara$a; b = hyperPara$b
    cd = hyperPara$cd; aa1=2; aa2=1
    alpha0 = hyperPara$alpha0; alpha1 = hyperPara$alpha1
    cutoff = hyperPara$cutoff; h = hyperPara$h


    ## initiations

    a.mat <- matrix(0,nrow=n,ncol=G)
    b.mat <- matrix(0,nrow=n,ncol=G)
    sig <- rep(1,G)
    for(g in 1:G){
        #print(i)
        fit <- tryCatch(lmer(Y[,g]~(1+years|ID)))
        a.mat[,g] <- coef(fit)$ID[,2]
        b.mat[,g] <- coef(fit)$ID[,1]
        sig[g] <- summary(fit)$sigma^2
    }

    siga <- apply(a.mat,2,var)
    sigb <- apply(b.mat,2,var)
    sig[sig<1e-4] <- 0.01
    siga[siga<1e-4] <- 0.01
    sigb[sigb<1e-4] <- 0.01

    M = 10
    e <- sample(1:as.integer(n/5),n,replace = TRUE)


    c1 = 10
    c2 = 10
    delta <- c(rgamma(1,c1,1),rgamma(M-1,c2,1))
    tau <- cumprod(delta)
    phi <- matrix(rgamma(G*M,v,v),ncol=M)

    lamA <- sapply(1:M,function(i){
        rnorm(G,0,1/sqrt(phi[,i]*tau[i]) )
    })
    etaA <- matrix(rnorm(n*M),nrow=M)

    B <- b.mat
    A <- a.mat

    muA0 <- apply(a.mat,2,mean)

    sig0a = 1
    sigma0 <- diag(rep(sig0a,G))
    sigma0Inv <- diag(rep(1/sig0a,G))
    b0 <- apply(b.mat,2,mean)

    sigma = lamA%*%t(lamA) + diag(siga)
    sigmaInv <- chol2inv(chol(sigma))

    pt = exp(alpha0+alpha1*(1:iter))
    ut <- runif(iter)

    ## construct data to store results
    e.mat <- matrix(0,n,iter)
    if(savePara){
        c.vec <- rep(0,iter)

        sig0a.vec <- rep(0,iter)
        muA.vec <- matrix(0,G*n,iter)
        muA0.vec <- matrix(0,G,iter)
        b0.vec <- matrix(0,G,iter)

        sigb.vec <- matrix(0,G,iter)
        sig.vec <- matrix(0,G,iter)
        A.vec <- matrix(0,n*G,iter)
        B.vec <- matrix(0,n*G,iter)
        sigma.vec <- matrix(0,G,iter)
    }

    for( i in 1:iter){
        print(paste("iterations ",i))

        for(thinning in 1:thin){

            e <- polyurncppInt(e,muA0,sigma0, A, sigma, sigmaInv,sigma0Inv,c)
            e <- match(e,sort(unique(e)))

            muA <- samMuA(e,A,muA0,sigma0Inv,sigmaInv)
            muA0 <- samMuA0(sigma0Inv,muA,h)

            cat("1\n")
            c <- samCv2(a,b,e,c,cd,n)

            ### sampling lamA and etaA
            lamA <- samLamV2Cpp(A - t(muA[,e]),  etaA, siga, lamA, phi, tau)

            if(ut[i]<=pt[i]){
                ind <- apply(lamA,2,function(x){all(abs(x)<cutoff)})
                if(sum(ind)>0){
                    lamA <- lamA[,!ind]
                    phi <- phi[,!ind]
                    delta <- delta[!ind]
                    tau <- cumprod(delta)
                }else{
                    delta <- c(delta,rgamma(1,c2,1))
                    tau <- cumprod(delta)
                    phi <- cbind(phi,rgamma(G,v,v))
                    lamA <- cbind(lamA,rnorm(G,0,1/sqrt(phi[,ncol(phi)]*tau[ncol(phi)]) ) )
                }
            }
            etaA <- samEtaA(siga,A,lamA,muA,e,1)
            etaA <- t(scale(t(etaA)))
            phi <- samPhi(v,v,lamA,tau)
            delta <- samDelta(c1,c2,phi,lamA,tau,delta)
            tau <- cumprod(delta)
            c1 <- samaV2(c1,delta[1])
            c2 <- samaV2(c2,delta[-1])

            ### sampling A
            A <- samA(Y,years,sig,B,muA,ID,e,sigmaInv)
            siga <- samSigA(v1,v2,muA,e,A,lamA,etaA)


            ### sampling B
            B <- samBInt(Y,years,sig,sigb,A,b0,ID)
            b0 <- SamB0Int(sigb,h,B)
            sigb <- samSigBInt(v1,v2,b0,B)

            ## sampling \Sigma
            sig <- samSig(v1,v2,Y,years,ID,B,A)
            sig0a <- samSig0(v1,v2[1],muA0,muA)

            sigma = lamA%*%t(lamA) + diag(siga)
            sigmaInv <- chol2inv(chol(sigma))
            sigma0 <- diag(rep(sig0a,G))
            sigma0Inv <- diag(rep(1/sig0a,G))

        }

        print(e)
        e.mat[,i] <- e
        if(savePara){
            c.vec[i] <- c
            muA.vec[,i] <- as.vector(muA[,e])
            muA0.vec[,i] <- muA0
            sig0a.vec[i] <- sig0a

            b0.vec[,i] <- as.vector(b0)
            sigb.vec[,i] <- sigb
            sigma.vec[,i] <- diag(sigma)
            sig.vec[,i] <- sig
            A.vec[,i] <- as.vector(A)
            B.vec[,i] <- as.vector(B)

        }

    }


    ## return results
    if(savePara){
        list(e.mat=e.mat,muA.vec=muA.vec, muA0.vec=muA0.vec,sig0a.vec=sig0a.vec,
             sig.vec=sig.vec, b0.vec=b0.vec, sigb.vec=sigb.vec, A.vec=A.vec,
             B.vec=B.vec,sigma.vec=sigma.vec,c.vec=c.vec)
    }else{
        list(e.mat=e.mat)
    }


}




#' Internal function to perform clustering using both intercepts and slopes
#' but not using factor analysis model.
#' @noRd
sparseBothNof <- function(data=NULL, iter=100, thin = 2, savePara=FALSE,
                          hyperPara=list(v1=0.1, v2=0.1, v=1.5, c=1, a=0, b=10,
                                         cd=1, aa1=2, aa2=1, alpha0= -1,
                                         alpha1 = -1e-4, cutoff = 1e-4, h =100)
){
    ## data section
    Y = data$Y
    ID = data$ID
    ID = match(ID,ID[!duplicated(ID)])
    years = data$years

    Y <- scale(Y,scale=TRUE)
    G <- ncol(Y)
    n <- max(ID)
    X <- cbind(1,years) ## years has unit scale

    ## hyperparameters
    v1 = hyperPara$v1; v2 = hyperPara$v2; v = hyperPara$v
    c = hyperPara$c; a = hyperPara$a; b = hyperPara$b
    cd = hyperPara$cd; aa1=2; aa2=1
    alpha0= hyperPara$alpha0; alpha1 = hyperPara$alpha1
    cutoff = hyperPara$cutoff; h = hyperPara$h

    ## initiations

    a.mat = matrix(0,nrow=n,ncol=G)
    b.mat = matrix(0,nrow=n,ncol=G)
    sig = rep(1,G)
    for(g in 1:G){
        #print(i)
        fit = tryCatch(lmer(Y[,g]~(1+years|ID)))
        a.mat[,g] = coef(fit)$ID[,2]
        b.mat[,g] = coef(fit)$ID[,1]
        sig[g] = summary(fit)$sigma^2
    }

    siga = apply(a.mat,2,var)
    sigb = apply(b.mat,2,var)

    sig[sig<1e-4] = 0.01
    siga[siga<1e-4] = 0.01
    sigb[sigb<1e-4] = 0.01

    e <- sample(1:as.integer(n/5),n,replace = TRUE)

    B = b.mat
    A = a.mat

    muB0 = apply(b.mat,2,mean)
    muA0 = apply(a.mat,2,mean)

    sig0b = 1
    sig0a = 1


    sigmaA0 = diag(rep(sig0a,G))
    sigmaA0Inv = diag(rep(1/sig0a,G))
    sigmaB0 = diag(rep(sig0b,G))
    sigmaB0Inv = diag(rep(1/sig0b,G))


    sigmaB = diag(sigb)
    sigmaBInv = diag(1/sigb)
    sigmaA =  diag(siga)
    sigmaAInv = diag(1/siga)

    ## construct data to store results
    e.mat = matrix(0,n,iter)
    if(savePara){

        c.vec <- rep(0,iter)
        sigmaA.vec <- matrix(0,G,iter)
        sigmaB.vec <- matrix(0,G,iter)

        sig0b.vec <- rep(0,iter)
        sig0a.vec <- rep(0,iter)
        sig.vec <- matrix(0,G,iter)

        muB.vec <- matrix(0,G*n,iter)
        muB0.vec <- matrix(0,G,iter)
        muA.vec <- matrix(0,G*n,iter)
        muA0.vec <- matrix(0,G,iter)

        A.vec <- matrix(0,n*G,iter)
        B.vec <- matrix(0,n*G,iter)
    }

    ### running MCMC

    for( i in 1:iter){
        print(paste("iterations ",i))

        for(thinning in 1:thin){

            e <- polyurncppBoth(e,A, muA0, sigmaA,  sigmaAInv, sigmaA0,sigmaA0Inv,
                                B, muB0, sigmaB, sigmaBInv, sigmaB0,sigmaB0Inv, c)
            e <- match(e,sort(unique(e)))
            c <- samCv2(a,b,e,c,cd,n)

            ## DP for B
            muB <- samMuB(e,B,muB0,sigmaB0Inv,sigmaBInv)
            muB0 <- samMuB0(sigmaB0Inv,muB,h)


            #DP for A
            muA <- samMuA(e,A,muA0,sigmaA0Inv,sigmaAInv)
            muA0 <- samMuA0(sigmaA0Inv,muA,h)

            ### sampling A
            A <- samA(Y,years,sig,B,muA,ID,e,sigmaAInv)
            siga <- samSigA(v1,v2,muA,e,A,lamA=matrix(0,nrow=G,ncol=1),
                            etaA=matrix(0,nrow=1,ncol=n))
            sig0a <- samSig0(v1,v2[1],muA0,muA)

            ### sampling B
            B <- samB(Y,years,sig,sigmaBInv,A,muB,e,ID)
            sigb <- samSigB(v1,v2,muB,e,B,lamB=matrix(0,nrow=G,ncol=1),
                            etaB=matrix(0,nrow=1,ncol=n))
            sig0b <- samSig0(v1,v2[1],muB0,muB)

            sig <- samSig(v1,v2,Y,years,ID,B,A)


            sigmaB = diag(sigb)
            sigmaBInv <- diag(1/sigb)

            sigmaA =  diag(siga)
            sigmaAInv <- diag(1/siga)

            sigmaA0 <- diag(rep(sig0a,G))
            sigmaA0Inv <- diag(rep(1/sig0a,G))

            sigmaB0 <- diag(rep(sig0b,G))
            sigmaB0Inv <- diag(rep(1/sig0b,G))
        }


        if(savePara){
            c.vec[i] <- c

            muB.vec[,i] <- as.vector(muB[,e])
            muB0.vec[,i] <- muB0
            muA.vec[,i] <- as.vector(muA[,e])
            muA0.vec[,i] <- muA0

            sig0b.vec[i] <- sig0b
            sig0a.vec[i] <- sig0a
            sigmaA.vec[,i] <- diag(sigmaA)
            sigmaB.vec[,i] <- diag(sigmaB)

            sig.vec[,i] <- sig
            A.vec[,i] <- as.vector(A)
            B.vec[,i] <- as.vector(B)
        }

        print(e)
        e.mat[,i] <- e


    }



    ## return results
    if(savePara){
        list(e.mat=e.mat, muA.vec=muA.vec, muA0.vec=muA0.vec, sig0a.vec=sig0a.vec,
             muB.vec=muB.vec, muB0.vec=muB0.vec, sig0b.vec=sig0b.vec, sig.vec=sig.vec,
             A.vec=A.vec, B.vec=B.vec, sigmaA.vec=sigmaA.vec, sigmaB.vec=sigmaB.vec,
             c.vec=c.vec)
    }else{
        list(e.mat=e.mat)
    }



}



#' Internal function to perform clustering using only intercepts
#' but not using factor analysis model.
#' @noRd
sparseIntNof <- function(data=NULL, iter=100, thin = 2, savePara=FALSE,
                         hyperPara=list(v1=0.1, v2=0.1, v=1.5, c=1, a=0, b=10,
                                        cd=1, aa1=2, aa2=1, alpha0= -1,
                                        alpha1 = -1e-4, cutoff = 1e-4, h =100)
){
    ## data section
    Y = data$Y
    ID = data$ID
    ID = match(ID,ID[!duplicated(ID)])
    years = data$years

    Y <- scale(Y,scale=TRUE)
    G <- ncol(Y)
    n <- max(ID)
    X <- cbind(1,years) ## years has unit scale

    ## hyperparameters
    v1 = hyperPara$v1; v2 = hyperPara$v2; v = hyperPara$v
    c = hyperPara$c; a = hyperPara$a; b = hyperPara$b
    cd = hyperPara$cd; aa1=2; aa2=1
    alpha0 = hyperPara$alpha0; alpha1 = hyperPara$alpha1
    cutoff = hyperPara$cutoff; h = hyperPara$h

    ## initiations

    a.mat <- matrix(0,nrow=n,ncol=G)
    b.mat <- matrix(0,nrow=n,ncol=G)
    sig <- rep(1,G)
    for(g in 1:G){
        #print(i)
        fit <- tryCatch(lmer(Y[,g]~(1+years|ID)))
        a.mat[,g] <- coef(fit)$ID[,2]
        b.mat[,g] <- coef(fit)$ID[,1]
        sig[g] <- summary(fit)$sigma^2
    }

    siga <- apply(a.mat,2,var)
    sigb <- apply(b.mat,2,var)
    sig[sig<1e-4] <- 0.1
    siga[siga<1e-4] <- 0.1
    sigb[sigb<1e-4] <- 0.1

    e <- sample(1:as.integer(n/5),n,replace = TRUE)
    B <- b.mat
    A <- a.mat

    muA0 <- apply(a.mat,2,mean)
    sig0a = 1
    sigma0 <- diag(rep(sig0a,G))
    sigma0Inv <- diag(rep(1/sig0a,G))
    b0 <- apply(b.mat,2,mean)
    sigma =  diag(siga)
    sigmaInv <-  diag(1/siga)

    ## construct data to store results
    e.mat <- matrix(0,n,iter)
    if(savePara){
        c.vec <- rep(0,iter)
        sig0a.vec <- rep(0,iter)
        muA.vec <- matrix(0,G*n,iter)
        muA0.vec <- matrix(0,G,iter)
        b0.vec <- matrix(0,G,iter)

        sigb.vec <- matrix(0,G,iter)
        sig.vec <- matrix(0,G,iter)
        A.vec <- matrix(0,n*G,iter)
        B.vec <- matrix(0,n*G,iter)
        sigma.vec <- matrix(0,G,iter)
    }

    for( i in 1:iter){
        print(paste("iterations ",i))

        for(thinning in 1:thin){

            e <- polyurncppInt(e,muA0,sigma0, A, sigma, sigmaInv,sigma0Inv,c)
            e <- match(e,sort(unique(e)))

            muA <- samMuA(e,A,muA0,sigma0Inv,sigmaInv)
            muA0 <- samMuA0(sigma0Inv,muA,h)

            c <- samCv2(a,b,e,c,cd,n)

            ### sampling A
            A <- samA(Y,years,sig,B,muA,ID,e,sigmaInv)
            siga <- samSigA(v1,v2,muA,e,A,lamA=matrix(0,nrow=G,ncol=1),etaA=matrix(0,nrow=1,ncol=n))


            ### sampling B
            B <- samBInt(Y,years,sig,sigb,A,b0,ID)
            b0 <- SamB0Int(sigb,h,B)
            sigb <- samSigBInt(v1,v2,b0,B)

            ## sampling \Sigma
            sig <- samSig(v1,v2,Y,years,ID,B,A)
            sig0a <- samSig0(v1,v2[1],muA0,muA)

            sigma =  diag(siga)
            sigmaInv <- diag(1/siga)
            sigma0 <- diag(rep(sig0a,G))
            sigma0Inv <- diag(rep(1/sig0a,G))

        }

        print(e)
        e.mat[,i] <- e

        if(savePara){
            c.vec[i] <- c
            muA.vec[,i] <- as.vector(muA[,e])
            muA0.vec[,i] <- muA0
            sig0a.vec[i] <- sig0a

            b0.vec[,i] <- as.vector(b0)
            sigb.vec[,i] <- sigb
            sigma.vec[,i] <- diag(sigma)

            sig.vec[,i] <- sig
            A.vec[,i] <- as.vector(A)
            B.vec[,i] <- as.vector(B)

        }

    }

    ## return results
    if(savePara){
        list(e.mat=e.mat,muA.vec=muA.vec, muA0.vec=muA0.vec,sig0a.vec=sig0a.vec,
             sig.vec=sig.vec, b0.vec=b0.vec, sigb.vec= sigb.vec, A.vec=A.vec,
             B.vec=B.vec, sigma.vec=sigma.vec,c.vec=c.vec)
    }else{
        list(e.mat=e.mat)
    }


}


#' Simulated dataset for testing the algorithm
#' @docType data
#' @usage data(data)
#' @keywords datasets
#' @examples
#' data(data)
#' ## this is the required data input format
#' head(data.frame(ID=data$ID,years=data$years,data$Y))
"data"


