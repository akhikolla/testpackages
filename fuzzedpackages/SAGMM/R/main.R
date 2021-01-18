#'SAGMM: A package for Clustering via Stochastic Approximation and Gaussian Mixture Models.
#'
#'The SAGMM package allows for computation of gaussian mixture models using stochastic approximation to increase efficiency with large data sets.
#'The primary function \code{SAGMMFit} allows this to be performed in a relative flexible manner.
#'@author Andrew T. Jones and Hien D. Nguyen
#'@references Nguyen & Jones (2018). Big Data-Appropriate Clustering via Stochastic Approximation and Gaussian Mixture Models. In Data Analytics (pp. 79-96). CRC Press.
#'@docType package
#'@name SAGMM
NULL

#' @import stats MixSim mclust lowmemtkmeans
NULL

#' Return Gamma, a sequence of gain factors
#' 
#' @description  Generate a series of gain factors.
#' @param Number Number of values required.
#' @param Burnin Number of 'Burnin' values at the beginning of sequence.
#' @return Gamma, a vector of gain factors.
#' @examples
#' g<-gainFactors(10^4, 2*10^3)
#' @export
gainFactors <- function(Number, Burnin) {
    # Make a sequence of gain factors
    Gamma <- c(rep(log(Number)/(Number),round(Number/Burnin)),
               rep(1/Number,Number-round(Number/Burnin))) 
    return(Gamma)
}

#' Generate data for simulations to test the SAGMM package..
#' 
#' @description  This function is primarily a convienence wrapper for MixSim.
#' @param ngroups Number of mixture components. Default 5.
#' @param Dimensions number of Dimensions. Default 5.
#' @param Number number of samples. Default 10^4.
#' @return List of results: X, Y, simobject.
#' @examples
#' sims<-generateSimData(ngroups=10, Dimensions=10, Number=10^4)
#' sims<-generateSimData()
#'@export
generateSimData<-function(ngroups=5, Dimensions=5, Number=10^4){
    MS <- MixSim::MixSim(BarOmega=0.01,K=ngroups,p=Dimensions,PiLow=(0.1/ngroups)) # Simulation code
    Z <- MixSim::simdataset(n=Number,Pi=MS$Pi,Mu=MS$Mu,S=MS$S) # More simulation code, look at package vignette
    X <- Z[[1]] # Extract features
    Y <- Z[[2]] # Extract Labels
    SAMPLE <- sample(1:Number) # Randomize
    X <- X[SAMPLE,] # Randomize
    Y <- Y[SAMPLE] # Randomize
    return(list(X=X, Y=Y, MS=MS))
}


#' Clustering via Stochastic Approximation and Gaussian Mixture Models (GMM)
#' 
#' @description Fit a GMM via Stochastic Approximation. See Reference.
#' @param X numeric matrix of the data.
#' @param Y Group membership (if known). Where groups are integers in 1:ngroups. If provided ngroups can 
#' @param Burnin Ratio of observations to use as a burn in before algorithm begins.
#' @param ngroups Number of mixture components. If Y is provided, and groups is not then is overridden by Y.
#' @param kstart number of kmeans starts to initialise.
#' @param plot If TRUE generates a plot of the clustering.
#'@return A list containing
#'\item{Cluster}{The clustering of each observation.}
#'\item{plot}{A plot of the clustering (if requested).}
#'\item{l2}{Estimate of Lambda^2}
#'\item{ARI1}{Adjusted Rand Index 1 - using k-means}
#'\item{ARI2}{Adjusted Rand Index 2 - using GMM Clusters}
#'\item{ARI3}{Adjusted Rand Index 3 - using intialiation k-means}
#'\item{KM}{Initial K-means clustering of the data.}
#'\item{pi}{The cluster proportions (vector of length ngroups)}
#'\item{tau}{tau matrix of conditional probabilities.}
#'\item{fit}{Full output details from inner C++ loop.}
#' @examples
#' sims<-generateSimData(ngroups=10, Dimensions=10, Number=10^4)
#' res1<-SAGMMFit(sims$X, sims$Y)
#' res2<-SAGMMFit(sims$X, ngroups=5)
#'@author Andrew T. Jones and Hien D. Nguyen
#'@references Nguyen & Jones (2018). Big Data-Appropriate Clustering via Stochastic Approximation and Gaussian Mixture Models. In Data Analytics (pp. 79-96). CRC Press.
#'@export
SAGMMFit<-function(X, Y=NULL, Burnin=5, ngroups=5, kstart=10, plot=FALSE){

    Number<-nrow(X) # N observations
    Dimensions <-ncol(X) #dim of data
    
    if(length(Y)>0){
      if(ngroups==0 |class(ngroups)!="numeric"){
        stop("At least one of ngroups or Y must be provided")
      }
      if(length(Y)!=nrow(X)){
        stop("Y length is not equal number of rows of X.")
      }
    }else{
      if(ngroups==0 |class(ngroups)!="numeric"){
        ngroups <- max(Y)
      }
    }

    
    ### Initialize Algorithm
    KM <- suppressWarnings(stats::kmeans(X[1:round(Number/Burnin),],ngroups,nstart=kstart)) # Use K means on burnin sample
    MU <- KM[[2]] # Get initial MU
    LAMBDA <- rep(max(sqrt(KM$withinss/KM$size)),ngroups) # Get initial lambda
    SIGMA <- list() # Make sigma from Lambda
    SIGMA_C<-array(0,c(Dimensions,Dimensions,ngroups))
    
    for (gg in 1:ngroups) {
      SIGMA[[gg]] <- diag(Dimensions)*LAMBDA[gg]^2/2
      SIGMA_C[,,gg]<-SIGMA[[gg]]
    }
    
    PI <- rep(1/ngroups,ngroups) # Get initial PI
    PISTAR <- rep(0,ngroups) # Solve for initial Pistar
    for (it in 1:100) {
      for (gg in 1:(ngroups-1)) {
        PISTAR[gg] <- log((1/PI[gg]-1)^-1*sum(exp(PISTAR[-gg])))
      }
    }

    GAMMA<-gainFactors(Number, Burnin)

    ### Old stuff
    LAMBDA_O <- LAMBDA # old value of lambda
    MU_O <- MU # Old value of Mu
    PISTAR_O <- PISTAR # Old value of Pistar
  
    #MAIN ACTION HERE
    results<-main_loop(X, Dimensions, Number, ngroups, MU_O, LAMBDA_O,PISTAR_O, SIGMA_C, GAMMA)
    
    TauMAT<-results$Tau
    PI <-results$PI
    MU <-results$MU
    
    SIGMA <-list()
    for (gg in 1:ngroups) {
        SIGMA[[gg]] <- diag(Dimensions)*results$LAMBDA[gg]^2/2
    }
    
  
    Cluster <- apply(TauMAT,1,which.max)
    if(plot &(length(Y)>0)){
        p1<-plot(as.data.frame(X),col=Cluster,pch=Y)
    }else{
      if(plot){
        p1<-plot(as.data.frame(X),col=Cluster)  
      }else{
        p1<-NA
      }
        
    }
    
    l2<-LAMBDA^2
    pi <- sort(PI)
    KM <- kmeans(X,ngroups,nstart=10)
    
    if(length(Y)>0){
      ARI1<-adjustedRandIndex(lowmemtkmeans::nearest_cluster(X,KM$centers),Y)
      ARI2<-adjustedRandIndex(Cluster,Y)
      ARI3<-adjustedRandIndex(KM$cluster,Y)
    }else{
      ARI1<-NA
      ARI2<-NA
      ARI3<-NA
    }
    retList <-list(Cluster=Cluster, plot=p1, l2=l2, ARI1 = ARI1, ARI2 = ARI2, ARI3=ARI3, KM=KM, pi=pi, tau=TauMAT, fit=results)
    
    return(retList)
}