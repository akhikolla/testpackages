



# not exported
reconstruct.vector<-function(x,index,n)
{
  if(length(x) == n+1) { return(x) }
  y<-rep(NA,n+1)
  y[1]=x[1]
  y[index+1]=x[-1]

  ##deal with start
  i=1
  while(is.na(y[i+1])) i=i+1
  # if(i<n) y[2:(i+1)]=0 
  if(i<n) y[2:(i+1)]=x[2] 
  
  for(i in 2:n){
      if(is.na(y[i+1])) y[i+1]=y[i] 
  }  
  return(y)

}

# not exported
reconstruct.pelt.result<-function(pelt.out,index,n)
{
   if(length(pelt.out$lastchangelike) == n+1) { return(pelt.out) }
   #update the output of F
    F<-reconstruct.vector(pelt.out$lastchangelike,index,n)
   cpts=pelt.out$cpts
    if(length(cpts)>1){
        cpts=index[cpts]
    }
   cpts[length(cpts)]<-n # last changepoint is always the end of the data   
   return(list("lastchangelike"=F,"cpts"=cpts))    
}


#' Most Recent Changepoints. 
#'
#' Detects the Most Recent Changepoints (mrc) for panel data consisting of many related univariate timeseries (Bardwell, Eckley, Fearnhead and Smith, 2016).
#' The method first determines the most recent univariate changepoints using PELT (Killick, Fearnhead and Eckley 2012) and then pools information across the
#' time-series by solving the K-median problem using \link[tbart]{tb.raw} (Teitz and Bart, 1968).  
#'
#' @param data An \eqn{n \times N} matrix or data frame representing a length \eqn{n} time series containing observations of \eqn{N} variables. The data can contain missing values
#' which are indicated by NA.
#' @param cost A string indicating which cost function to use. Possible choices are "mean" (change in mean) or "var" (change in variance).
#' The default value is \code{cost="mean"}.
#' @param alpha The variable-specific penalty used to penalise the addition of a given changepoint into a given variable. This can be specified as a real positive value
#' or as a function of \eqn{n}. The latter form is used when the data contains missing values which leads to time series of different lengths.
#'
#' Default value \code{alpha=function(n) 1.5*log(n)}.
#' @param indexed Boolean value indicating that the first column of \code{data} is an index variable. If \code{indexed=FALSE} an index variable will automatically be generated.
#' Default value is \code{indexed=FALSE}.
#' @param pmax Maximum number of most recent changepoints to search for. Default value \code{pmax=5}.
#' @param mad Boolean value indicating if the variates should be scaled by an estimate of the their standard deviation obtained using mean absolute deviation (Zhang, Nancy, Siegmund and David 2007).
#' This is useful in conjunction with \code{cost="mean"} for which unit variance is assumed. Default value is \code{mad=FALSE}.
#' @param phi Lag 1 autocorrelation to model the temporal dependence in the residuals of the time series assuming a MA(1) process. Default \code{phi=0}. 
#' 
#' @references \insertRef{doi:10.1080/00401706.2018.1438926}{changepoint.mv}
#' @references \insertRef{OR:TEITZBART}{changepoint.mv}
#' @references \insertRef{doi:10.1080/01621459.2012.737745}{changepoint.mv} 
#' 
#' @examples
#' library(changepoint.mv)
#' data(mrcexample)
#' res<-mrc(mrcexample,pmax=2)
#' MDL(res)     # MDL == pmax (possibly under-estimating MDL, retry)
#' res<-mrc(mrcexample,pmax=6)
#' MDL(res)     # MDL = 5 (< pmax)
#' # view the most recent changepoints (corresponding to pmax = 5)
#' unique(cpts.mr(res,p=5)[,1])
#' summary(res) # summary of result
#'
#' @export
mrc <- function(data,cost="mean",alpha=function(n) 1.5*log(n),pmax=5,indexed=FALSE,mad=FALSE,phi=0.0)
{
    # no of series
    n <- nrow(data)    
    usr.data<-process.data(data,indexed)
    
    penalty_function<-NULL
    # set the penalty function
    if(is_function(alpha))
    {
        penalty_function<-alpha 
    }
    else if(is_numeric(alpha))
    {
        penalty_function<-function(n) {alpha}
    }
    else 
    {
        stop("alpha should be scalar or a function") 
    }

    

    # phi correction - only moving average for now
    correction<-NULL
    if(is.na(phi))
        {
            correction<-1
        }
    else if(is_scalar(phi) && 0 <= phi && phi <= 1)
        {
            correction<-1+2*phi
        }
    else
        {
            stop("phi must be a scalar with 0 <= phi <= 1")
        }
    
    if(!(cost == "mean" || cost == "var"))
    {
        stop("mrc only supports cost = mean or var in this version.")
    }

    min.dist<-0
    if(cost=="var")
    {
        min.dist<-1
    }

    
    # apply mad if necessary
    if(mad)
    {
        data<-cbind(usr.data[,1],data.frame(Map(function(i) usr.data[,i]*sqrt(2)/mad(diff(usr.data[,i][!is.na(usr.data[,i])])),2:ncol(usr.data))))
        names(data)<-names(usr.data)
    }
    else
    {
        data<-usr.data
    }
    data<-as.matrix(data[-1])


    not.missing.indices<-Map(function(col) (1:nrow(data))[!is.na(data[,col])],1:ncol(data))
    
    pelt.results<-tryCatch(
        {
            Map(function(X,index)
            {
                pelt.result<-rcppeigen_peltuv(X,cost,correction*penalty_function(length(X)),min.dist)
                pelt.result$lastchangelike[1:length(X)]<-pelt.result$lastchangelike[1:length(X)]+rcppeigen_tail_costs(X,cost,min.dist)
                pelt.result$cpts<-pelt.result$cpts[-1]
                pelt.result<-reconstruct.pelt.result(pelt.result,index,nrow(data))
                return(list(pelt.result$lastchangelike[-(nrow(data)+1)],pelt.result$cpts))
            },
            Map(function(u) data[,u][!is.na(data[,u])],1:ncol(data)),
            not.missing.indices		 				 
            )
        },
        error = function(e) {e$message<-"mrc stopped because of user interrupt";stop();}
        )
     
    G <- matrix(unlist(Map(function(X) X[[1]],pelt.results)),ncol(data),nrow(data),byrow=TRUE)
    if(!is.na(phi))
    {
        G <- G/correction 
    }
    uv.cpts <- Map(function(X) X[[2]] ,pelt.results)

    N <- dim(G)[1]
    location.vec <- 0:(n-1)
    
    penalised.cost <- numeric( pmax )
    mmrc <- vector( "list" , pmax )
    affected <- vector( "list" , pmax )
    # p=1 separate as simpler
    mmrc[[1]] <- tb.raw( G , c(1) ) 
    # say all series are affected by this 1 change
    index <- rep(1,times=N)
    # find which affected (more evidence above threshold)
    affected[[1]] <- index
    # objective cost
    penalised.cost[1] <-  sum( G[ , mmrc[[1]] ] )
    
    if (pmax>1){
        
        for (p in 2:pmax){
            
            # mmrc[[i]] gives locations of the i best locations
            mmrc[[p]] <- tb.raw( G , c(1:p) )
                         0               # affected[[i]], gives each dimension a label from 1:i
                                        # depending on which change it is associated with
            affected[[p]] <- apply( G[ , mmrc[[p]] ] , 1 , which.min )
            # penalised.cost[[i]] gives the objective cost for solving with i different changes/sets
            csum <- 0
            for (i in 1:p){
                wa <- which( affected[[p]] == i)
                csum <- csum + sum( G[ wa , mmrc[[p]][i] ] )  
            }
            penalised.cost[p] <- csum
            
        }
        
    }
    
    locations <- vector("list",pmax)
    for ( i in 1:pmax ){
        locations[[i]] <- location.vec[ mmrc[[i]] ]
    }
    MDL <- which.min(penalised.cost + (N*log(1:pmax)+ (1:pmax)*log(n))/log(2)) ##MDL
    cpts.uv<-Map(function(cpts) c(0,cpts),uv.cpts)
    cpts.mv<-list()
    for ( i in 1:pmax )
        {
            cpts.mr<-locations[[i]][affected[[i]]]
            cpts.mv[[i]]<-Map(function(uv,mr) if(length(uv) > 2) {uv[length(uv)-1]<-mr;return(uv);} else {return(uv);},cpts.uv,cpts.mr)
        }
    return(changepoint.mv.mrc.class(data=usr.data,pmax=pmax,cpts.uv=cpts.uv,cost=cost,mad=mad,cpts.mv=cpts.mv,alpha=penalty_function(n),
                                    penalised.costs=penalised.cost,locations=locations,affected=affected))
        
}

