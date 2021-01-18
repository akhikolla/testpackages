#' Prepare mouse-tracking trajectories for state-space modeling via Stan
#'
#' @details The function prepares the mouse-tracking trajectories to be modeled for the state-space analysis. It automatically processes trajectories according to time-normalization, translation, and atan2 conversion. 
#' Users can skip pre-processing by setting \code{preprocess=FALSE}.
#' 
#' The input dataframe \code{X} needs to be organized using the long format with information being organized as nested. In particular, \code{X} must contains the following variables:
#' \describe{
#' \item{sbj}{The ID number of participants}
#' \item{trial}{The ID number of trials}
#' \item{factors}{1,...,Q factors for the categorical variables of the design. They may have different levels.  }
#' \item{timestep}{The ID number of the recorded x-y trajectories}
#' \item{x}{The recorded x-trajectories associated to trials and experimental levels}
#' \item{y}{The recorded y-trajectories associated to trials and experimental levels}
#' }
#' 
#' See \code{\link{language}} and \code{\link{congruency}} as examples of datasets format required by \pkg{ssMousetrack} package.
#' 
#' @param X (dataframe) a data frame of x-y trajectories and experimental design (see \code{Details})
#' @param preprocess (boolean) indicates whether x-y trajectories should be pre-processed (default \code{preprocess=TRUE})
#' @param N (integer) number of timesteps for trajectory normalization (default \code{N=61})
#' @param Z.formula (character) a formula of the contrasts for the model matrix Z (see \code{\link{model.matrix}})
#' @param Z.contrast (character) type of contrasts (default: treatment) for the model matrix Z (see \code{\link{model.matrix}})
#' @param yT (numeric) position in angles of the target. The default option yT="AUTO" will automatically determine the target position from the observed data
#' @param yD (numeric) position in angles of the distractor. The default option yD="AUTO" will automatically determine the target position from the observed data
#' @return a list containing (i) the new dataframe of the pre-processed dataset (\code{X_processed}) and (ii) the needed data for \code{\link{run_ssm}}
#' @export
#' @examples
#' 
#' data(congruency)
#' dataout <- prepare_data(X = congruency,preprocess = TRUE,Z.formula = "~congruency*plausibility")
#' str(dataout)
#' 

prepare_data <- function(X=NULL,preprocess=TRUE,N=61,Z.formula=NULL,Z.contrast="treatment",yT="AUTO",yD="AUTO"){
  if(N<1)
    stop("Positive integers should be provided for N")
  if(is.null(X))
    stop("The input dataset X must be provided")
  if(!is.logical(preprocess))
    stop("preprocess can be either TRUE or FALSE")
  if(!is.data.frame(X))
    stop("X must be a dataframe")
  if(is.null(Z.formula))
    stop("The model.matrix formula for Z must be provided")
  if((yT!="AUTO") & (yD!="AUTO")){
    if(yT<=0 | yT > pi)
      stop("The parameter yT must lie in the range (0,pi]")
    if(yD<=0 | yD > pi)
      stop("The parameter yD must lie in the range (0,pi]")
    if(yT >= yD)
      stop("The parameter yD must be greater yT")
  }

  
  # Get information on factors
  fcts <- names(Filter(f = is.factor,x = (X))) #get all the factors from X
  iid <- which(names(X) %in% names(Filter(f = is.factor,x = (X)))) #position of the factors in X  
  sbjs <- unique(X$sbj)
  
  # Preprocessing
  if(preprocess){
    X_normalized <- data.frame()
    for(i in 1:length(sbjs)){ #loop inside sbj
      Xi <- X[X$sbj==sbjs[i],]
      lvls <- unique(Xi[,iid]);names(lvls)=rep("",length(iid))
      lvls <- as.matrix(lvls)
      for(j in 1:dim(lvls)[1]){ #loop inside factors
        iidx <- which(apply(mapply(function(h){Xi[,iid[h]]==lvls[j,h]},1:length(iid)),1,sum)==length(iid))
        Xj <- Xi[iidx,]
        trjs = unique(Xj$trial)
        for(k in 1:length(trjs)){ #loop inside trials
          Xk <- Xj[Xj$trial==trjs[k],]
          # normalize x-y trajectory
          x <- stats::approx(x = 1:length(Xk$x),y = Xk$x,n = N)$y #normalize x
          y <- stats::approx(x = 1:length(Xk$y),y = Xk$y,n = N)$y #normalize y
          # reflection of x
          b <- stats::lm(formula = x~seq(1,N))$coef[2]
          if(b<0){x <- x*-1}
          # reflection of y
          b <- stats::lm(formula = y~seq(1,N))$coef[2]
          if(b<0){y <- y*-1}
          # traslation of x
          x <- x-x[1]
          # traslation of x-y into [-1 1] x [0 1]
          y <- (0.1) + ((y-min(y)) / (max(y)-min(y))) * (1-(0.1))
          x <- x/max(x) 
          # populate with the new trajectory
          X_normalized <- rbind(X_normalized,data.frame(Xk[1,1:(dim(Xk)[2]-3)],timestep=1:N,x,y,row.names = NULL))
        }
      }
    }
    X <- X_normalized
  }
  
  
  # Compute radians and adding small white noise
  X$yrad <- atan2(X$y,X$x) + rnorm(dim(X)[1],0,0.01)    
  
  # Compute other quantities
  I <- length(unique(X$sbj))
  N <- max(X$timestep)
  J <- length(unique(X$trial[X$sbj==sbjs[1]])) # sbjs must have same number of trials 
  
  if(yT=="AUTO" & yD=="AUTO"){
    yT <- mean(X$x[X$timestep==N])
    yD <- abs(pi-yT)
  }
    
  
  Y <- matrix(data = X$yrad,nrow = N,ncol = J*I)
  Z <- stats::model.matrix(formula(Z.formula),data = X[X$timestep==N,],contrasts.arg = Z.contrast)
  D <- compute_D(Y=Y,y_T=yT,y_D=yD)
  
  
  
  dataout <- list(X_processed = X,
                  I = I,
                  N = N,
                  J = J,
                  Y = Y,
                  Z = Z,
                  D = D,
                  yT = yT,
                  yD = yD)
  
  return(dataout)
}
