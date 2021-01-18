#' Evaluate the adequacy of the state-space model to reproduce the observed data
#'
#' @details The function implements a simulated-based method for assessing the adequacy of the model to reproduce the observed data. In particular,
#' the function provides two type of model adequacy, i.e. overall (PA_ov) and by-subject (PA_sbj). In the overall case the function provides the total amount of data reconstruction 
#' based on the I x J x N matrix Y of observed data. By contrast, in the second case the function provides the adequacy of the model to reconstruct the individual-based set of data
#' as it works on the matrix J x N over i=1,...,I. Both the indices are in the range 0\% - 100\%, with 100\% indicating perfect fit. In addition, the function returns a by-subject distance-based index (Dynamic Timw Warp distance) between observed and reproduced trajectories using \code{\link[dtw:dtw]{dtw}} function.
#'   
#' @param ssmfit (list) output of \code{\link{run_ssm}} function
#' @param M (integer) number of replications
#' @param plotx (boolean) if \code{plotx=TRUE} the function returns a graphical representation for the fit indices
#' @return a datalist containing the adequacy indices
#' @export
#' @examples
#' 
#' \dontrun{
#' ## Fit a state-space model using simulated data 
#' # Generate mouse-tracking data for an univariate experimental design with K = 3 categorical levels, 
#' ## J = 12 trials, I = 5 subjects
#' X1 <- generate_data(I=5,J=12,K=3,Z.formula="~Z1")
#' iid <- 23 # keep just one dataset from the simulated set of datasets
#' # Run the state-space model on the chosen dataset
#' X1_fit <- run_ssm(N = X1$N,I = X1$I,J = X1$J,Y = X1$data$Y[iid,,],D = X1$data$D[iid,,],
#' Z = X1$data$Z,niter=100,nwarmup=25)
#' # Evaluate the state-space model
#' evaluate_ssm(ssmfit = X1_fit,M = 10,plotx=FALSE)
#' }

evaluate_ssm <- function(ssmfit=NULL,M=100,plotx=TRUE){
  if(is.null(ssmfit))
    stop("The ssmfit object must be provided")
  if(M<0)
    stop("Positive integers should be provided for M")
  
  I <- ssmfit$I; J <- ssmfit$J; N <- ssmfit$N; Z <- ssmfit$data$Z; Y <- ssmfit$data$Y; gfunction <- ssmfit$Gfunction
  
  if(gfunction=="logistic"){
    gfunc <- function(lb,r,b,x){lb + (r / ( 1 + exp(b - x)))}
  }else if(gfunction=="gompertz"){
    gfunc <- function(lb,r,b,x){lb + (r * exp(-b * exp(-x)))}
  }
  
  Y_m <- array(0,c(N,I*J,M))
  PA_ov <- rep(NA,M)
  PA_sbj <- matrix(NA,M,I)
  dtws <- matrix(NA,M,I*J)
  
  iid <- rep(seq(from=1,to=I),each=J) #indices for subjects
  lb <- 0.1; bnds <- matrix(1,I*J,1)%*%matrix(c(lb,pi-lb,(pi-lb)-lb),1,3)
  
  if(M < dim(ssmfit$params$gamma)[1]){
    iidm <- sample(1:dim(ssmfit$params$gamma)[1],M,FALSE)
  }else{
    iidm <- 1:dim(ssmfit$params$gamma)[1]
  }
  
  for(m in 1:length(iidm)){ # outer loop over m=1...M
    
    mu_m=matrix(NA,N,I*J)
    b <- Z%*%t(ssmfit$params$gamma[iidm[m],])
    Y_m[1,,m] <- Y[1,] # do not consider the step n=0 in the evaluation
    for(n in 2:N){ # innter loop over n=1...N
      xv <- rep(ssmfit$data$X[iidm[m],n,],each=J)
      mu_m[n,] <- gfunc(bnds[,1],bnds[,3],b,xv)
      expx <- exp(ssmfit$params$lambda*ssmfit$data$D[n,])
      kappa1 <- ssmfit$params$kappa_bnds[1] + ((expx-min(expx)) / (max(expx)-min(expx))) * (ssmfit$params$kappa_bnds[2]-ssmfit$params$kappa_bnds[1])
      Y_m[n,,m] <- mapply(function(k){CircStats::rvm(n=1,mean = mu_m[n,k],k = kappa1[k])},1:(I*J))
    }
    
    PA_ov[m] <- max(1-norm(Y_m[,,m]-Y)^2/norm(Y)^2,0) # overall percentage of reconstruction
    PA_sbj[m,] <- max(mapply(function(i){1-norm(Y_m[,iid==i,m]-Y[,iid==i])^2/norm(Y[,iid==i])^2},1:I),0) # by subject percentage of reconstruction
    dtws[m,] <- max(mapply(function(i){dtw::dtw(Y_m[,i,m],Y[,i])$normalizedDistance},1:(I*J)),0) # dtw distance
    
  }
  
  dataout <- list(dist = list(PA_ov = PA_ov,PA_sbj = PA_sbj, DTW = dtws), 
                  indices = list(PA_ov = mean(PA_ov),PA_sbj = mean(PA_sbj), DTW = mean(dtws)))
  
  if(plotx==TRUE){
    data_plot=data.frame(y=as.vector(PA_sbj),x=rep(1:I,each=M))
    g1 = ggplot(data=data_plot,aes_string(x = "x",y = "y",group="x")) + geom_violin() + stat_summary(aes_string(group="x"),fun.y=mean, geom="point", color="black", size=1) +
      geom_hline(yintercept =mean(PA_ov),linetype=2,alpha=0.6,show.legend = FALSE) + xlab("subject index") + ylab("PA sbj") + theme_bw()

    g2 = ggplot(data=data.frame(PA_ov),aes(PA_ov)) + geom_histogram(bins = 12,color="black", fill="white") + geom_vline(xintercept = mean(PA_ov),col="red") + theme_bw() +
      ylab("") + xlab("PA overall")

    data_plot=data.frame(apply(dtws,1,mean));names(data_plot)="dtw"
    g3 = ggplot(data=data_plot,aes(dtw)) + geom_histogram(bins = 12,color="black", fill="white") + geom_vline(xintercept = mean(data_plot$dtw),col="red") + theme_bw() +
      ylab("") + xlab("DTW index")

    print(plot_grid(g1,plot_grid(g2,g3,ncol=1),labels=c("A","B")))
  }
  
  return(dataout)
}
