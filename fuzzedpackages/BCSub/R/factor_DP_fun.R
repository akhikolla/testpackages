
#' Internal function to sample cluster-specific means of eta
#' @keywords internal
samMu <- function(e,eta,mu0,sigma0Inv,sigmaInv){

    K <- max(e)
    mat <- matrix(0,nrow=nrow(sigma0Inv),ncol=K)
    M <- sapply(seq_len(K),function(k){as.numeric(e==k)})


    tmp.sig <- lapply(seq_len(K), function(x){
        solve(sigma0Inv + sum(M[,x])*sigmaInv)
    })
    tmp.mu <- lapply(seq_len(K),function(k){
        tmp.sig[[k]]%*%(sigmaInv%*%eta%*%M[,k] + sigma0Inv%*%mu0) })
    mat[,seq_len(K)] <- sapply(seq_len(K),function(k){
        mvrnorm(1,tmp.mu[[k]],tmp.sig[[k]])})
    mat
}


#' Internal function to sample eta
#' @keywords internal
samEta <- function(sig,A,lam,mu,e,sigmaInv){

    #I <- diag(rep(1,ncol(lam)))/sig0a
    sigma.a <- diag(1/sig)

    ##
    tmp.sig <- solve(sigmaInv + t(lam)%*%sigma.a%*%lam)
    tmp.mu <- tmp.sig%*%(t(lam)%*%sigma.a%*%t(A) + sigmaInv%*%mu[,e])

    return( t(mvrnorm(nrow(A),rep(0,ncol(lam)),tmp.sig)) + tmp.mu )
}

#' Internal function to sample gene-specific variances
#' @keywords internal
samSig <- function(v1,v2,A,lam,eta){
    v1 <- v1 + nrow(A)/2
    A <- A - t(lam%*%eta)
    v2 <- v2 + apply(A^2,2,sum)/2
    return(1/rgamma(length(v2),v1,v2))
}


#' Internal function to sample rho
#' @keywords internal
samRho2 <- function(mu,rho2,vv=1,up=1,lo=0){

    loglik <- sum( dnorm(as.vector(mu),0,sqrt(rho2),log=T) )

    v <- log((rho2-lo)/(up-rho2))
    vp <- rnorm(1,v,vv)
    rho2p <- (lo+up*exp(vp))/(exp(vp)+1)

    loglikp <- sum( dnorm(as.vector(mu),0,sqrt(rho2p),log=T) )

    loglik <- sum(loglik)  + v -2*log(exp(v)+1)
    loglikp <- sum(loglikp) + vp -2*log(exp(vp)+1)

    dif <- exp(loglikp - loglik)
    if(runif(1) < dif){rho2 <- rho2p}
    return(rho2)
}


#' Internal function to sample variances for eta
#' @keywords internal
samSige <- function(v1,v2,eta,mu,e){
    v1 <- v1 + ncol(eta)/2
    eta <- eta - mu[,e]
    v2 <- v2 + apply(eta^2,1,sum)/2
    return(1/rgamma(length(v2),v1,v2))
}


#' A Bayesian semiparametric factor analysis model for subtype identification
#'  (Clustering).
#' @param A Data matrix with rows being subjects and columns being genes.
#' @param iter Total number of iterations (including burn-in period).
#' @param seq Posterior samples used for inference of cluster structure.
#' @param M Number of factors.
#' @return returns a list with following objects.
#' \item{CL}{Inferred cluster strucutre based on the posterior samples.}
#' \item{E}{A matrix with each column being the cluster structre at each iteration.}
#' @references A Bayesian Semiparametric Factor Analysis Model for Subtype
#'             Identification. Jiehuan Sun, Joshua L. Warren, and Hongyu Zhao.
#' @export
#' @examples
#' set.seed(1)
#' n = 100 ## number of subjects
#' G = 200 ## number of genes
#' SNR = 0 ## ratio of noise genes
#' ## loading matrix with four factors
#' lam = matrix(0,G,4)
#' lam[1:(G/4),1] = runif(G/4,-3,3)
#' lam[(G/4+1):(G/2),2] = runif(G/4,-3,3)
#' lam[(G/2+1):(3*G/4),3] = runif(G/4,-3,3)
#' lam[(3*G/4+1):(G),4] = runif(G/4,-3,3)
#' ## generate low-rank covariance matrix
#' sigma <- lam%*%t(lam) + diag(rep(1,G))
#' sigma <- cov2cor(sigma)
#' ## true cluster structure ##
#' e.true = c(rep(1,n/2),rep(2,n/2))
#'
#' ## generate data matrix ##
#' mu1 = rep(1,G)
#' mu1[sample(1:G,SNR*G)] = 0
#' mu2 <- rep(0,G)
#' A = rbind(mvrnorm(n/2,mu1,sigma),mvrnorm(n/2,mu2,sigma))
#'
#' ## factor analysis to decide the number of factors
#' \dontrun{
#' ev = eigen(cor(A))
#' ap = parallel(subject=nrow(A),var=ncol(A),rep=100,cent=.05)
#' nS = nScree(x=ev$values, aparallel=ap$eigen$qevpea)
#' M = nS$Components[1,3] ## number of factors
#' }
#' M = 4
#' ## run BCSub for clustering
#' iters = 1000 ## total number of iterations
#' seq = 600:1000 ## posterior samples used for inference
#' system.time(res <- BCSub(A,iter=iters,seq=seq,M=M))
#' res$CL ## inferred cluster structure
#'
#' ## calculate and plot similarity matrix
#' sim = calSim(t(res$E[,seq]))
#'
#' ## plot similarity matrix
#' x <- rep(1:n,times=n)
#' y <- rep(1:n,each=n)
#' z <- as.vector(sim)
#' levelplot(z~x*y,col.regions=rev(gray.colors(n^2)), xlab = "Subject ID",ylab = "Subject ID")
BCSub <- function(A=NULL,iter=1000,seq=200:1000,M=5){

    n <- nrow(A)
    G <- ncol(A)
    hcc <- hclust(as.dist(1-abs(cor(A))))

    tree <- cutree(hcc,k=M)
    ord <- sapply(1:M, function(i){sample(which(tree==i))[1]})
    ordd <- setdiff(1:G,ord)
    A <- A[,c(ord,ordd)]
    # print(ord)
    ### hyperparameters
    v1 <- 0.01
    v2 <- rep(0.01,G)
    c <- 1
    a <- 0
    b <- 10
    cd <- 1
    h <- 100
    sigl <- 1

    ### run the model
    #e <- sample(1:n)
    e <- rep(1,n)
    mu0 <- rep(0,M)

    lam <- matrix(rnorm(G*M),nrow=G,ncol=M)
    lam[upper.tri(lam)] <- 0
    diag(lam) <- 1

    eta <- matrix(rnorm(M*n),nrow=M)
    eta <- t(scale(t(eta),scale=F))

    sig <-  rep(1,G)
    sige <- rep(1,M)
    rho1 <- 1
    rho2 <- 1
    sigma = diag(M)
    sigmaInv = diag(M)
    sigma0 = diag(M)
    sigma0Inv = diag(M)

    diag(sigma) <- rho1
    diag(sigmaInv) <- 1/rho1
    diag(sigma0) <- rho2
    diag(sigma0Inv) <- 1/rho2
    ###
    e.mat <- matrix(1,n,iter)
    sigma.vec <- matrix(0,M*M,iter)
    eta.vec <- matrix(0,M*n,iter)
    LAM = matrix(0,nrow=G,ncol=M)
    ETA = matrix(0,nrow=M,ncol=n)
    SIG = rep(0,G)
    for( i in 1:iter){
        e <- polyurncpp(e,mu0,sigma0, t(eta), sigma, sigmaInv,sigma0Inv,c)
        e <- match(e,sort(unique(e)))
        mu <- samMu(e,eta,mu0,sigma0Inv,sigmaInv)
        eta <- samEta(sig,A,lam,mu,e,sigmaInv)
        lam <- samLamV3Cpp(A,eta,sig,sigl,lam)
        sig <- samSig(v1,v2,A,lam,eta)
        rho2 <- samRho2(mu,rho2,vv=1,up=2,lo=0)

        sige <- samSige(v1,v2[1],eta,mu,e)
        diag(sigma) <- sige
        diag(sigmaInv) <- 1/sige

        diag(sigma0) <- rho2
        diag(sigma0Inv) <- 1/rho2

        e.mat[,i] <- e
    }
    CL = maxpear(calSim(t(e.mat[,seq])),method="avg")$cl
    list(E=e.mat,CL=CL)
}

