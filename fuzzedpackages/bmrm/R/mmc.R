

#' Loss function for max-margin clustering
#' 
#' @export
#' @param x numeric matrix representing the dataset (one sample per row)
#' @param k an integer specifying number of clusters to find
#' @param minClusterSize an integer vector specifying the minimum number of sample per cluster. 
#'        Given values are reclycled if necessary to have one value per cluster.
#' @param groups a logical matrix for instance grouping (groups[i,j] TRUE when 
#'        sample i belong to group j).
#' @param minGroupOverlap an integer matrix specifyng the minimum number of instance
#'        per cluster for each group.
#' @param weight a weight vector for each instance
#' @return the loss function to optimize for max margin clustering of the given dataset
#' @import lpSolve
mmcLoss <- function(x, k=3L,minClusterSize=1L,groups=matrix(logical(0),nrow(x),0), 
                    minGroupOverlap=matrix(integer(0),k,ncol(groups)), 
                    weight=1/nrow(x)) {
  if (!is.matrix(x)) stop("x must be a numeric matrix")
  if (nrow(groups)!=nrow(x)) stop("nrow(groups)!=nrow(x)")
  if (nrow(minGroupOverlap)!=k) stop("nrow(minGroupOverlap)!=k")
  if (ncol(minGroupOverlap)!=ncol(groups)) stop("ncol(minGroupOverlap)!=ncol(groups)")
  if (is.null(rownames(x))) rownames(x) <- seq_len(nrow(x))
  weight <- rep_len(weight,nrow(x))
  minClusterSize <- rep_len(minClusterSize,k)
  
  cst <- merge(by="col",
    which(as.matrix(groups),arr.ind=TRUE),
    which(as.matrix(minGroupOverlap)>0,arr.ind=TRUE)
  )
  cst$cstId <- as.integer(factor((cst$col-1)*k + cst$row.y))
  cst$varId <- nrow(x)*(cst$row.y-1L) + cst$row.x
  
  lp.constraints <- list(
    dense.const = rbind(
        # constrain the instances to belong to one and only one cluster
        eq = cbind(seq_len(nrow(x)),seq_len(nrow(x)*k),1),
        # constrain on the minimum size of the clusters
        gt = cbind(as.vector(col(matrix(NA,nrow(x),k)))+nrow(x),seq_len(nrow(x)*k),1),
        # groups constraints
        grp = cbind(nrow(x)+k+cst$cstId,cst$varId,rep(1,length(cst$cstId)))
    ),
    const.dir = c(eq=rep_len("==",nrow(x)),gt=rep_len(">=",k),gt=rep_len(">=",sum(minGroupOverlap>0))),
    const.rhs = c(eq=rep_len(1,nrow(x)),gt=minClusterSize,grp=minGroupOverlap[minGroupOverlap>0])
  )
  
  function(w) {
    W <- matrix(w, ncol(x),k)
    F <- x %*% W
    
    # compute loss incured for predicting a given instance in a given cluster
    R <- weight*array(rowSums(pmax(1 + (rep(1,ncol(F)) %x% F) - as.vector(F),0))-1,dim(F))
    
    # find cluster assignment minimizing the loss and satifying the constraints
    Y <- local({
      opt <- lp("min",
              objective.in = as.vector(R),
              dense.const = lp.constraints$dense.const,
              const.dir = lp.constraints$const.dir,
              const.rhs = lp.constraints$const.rhs,
              binary.vec = seq_along(R)
      )
      if (opt$status!=0) stop("LP problem not feasible")
      matrix(opt$solution,nrow(R),ncol(R))
    })
    
    G <- 1-Y+F-rowSums(F*Y)
    G <- ifelse(G>0,1,0)
    G <- Y*rowSums(G) - G
    G <- G*weight
    
    w <- as.vector(W)
    lvalue(w) <- sum(R*Y)
    gradient(w) <- as.vector(crossprod(x,-G))
    attr(w,"Y") <- Y
    class(w) <- "mmcLoss"
    w
  }
}



#' Convenient wrapper function to solve max-margin clustering problem on a dataset
#' 
#' Solve max-margin clustering problem with multiple random starting points to avoid being trap by local minima.
#' The random starting points are determined by randomly assigning N0 samples to each cluster and solving for multi-class SVM
#' 
#' @export
#' @import parallel
#' @param x numeric matrix representing the dataset (one sample per row)
#' @param k an integer specifying number of clusters to find
#' @param N0 number of instance to randomly assign per cluster when determining a random starting point.
#'        The classification dataset it defines is used to train a multi-class SVM whose solution is used 
#'        as the starting point of current MMC iteration.
#' @param LAMBDA the complexity parameter for nrbm()
#' @param seeds the random seeds to use
#' @param nrbmArgsSvm arguments to nrbm() when solving for multi-class SVM problem
#' @param nrbmArgsMmc arguments to nrbm() when solving for max-margin clustering problem
#' @param mc.cores number of core to use when running the random iterations in parallel
#' @param ... additional arguments are passed to mmcLoss()
#' @return the MMC model matrix
#' @examples
#'    # -- Prepare a 2D dataset to cluster with an intercept
#'    x <- cbind(intercept=100,scale(data.matrix(iris[c(1,3)]),center=TRUE,scale=FALSE))
#' 
#'    # -- Find max-margin clusters
#'    y <- mmc(x,k=3,LAMBDA=0.001,minClusterSize=10,seeds=5)
#'    table(y,iris$Species)
#'    
#'    # -- Plot the dataset and the MMC decision boundaries
#'    gx <- seq(min(x[,2]),max(x[,2]),length=100)
#'    gy <- seq(min(x[,3]),max(x[,3]),length=100)
#'    Y <- outer(gx,gy,function(a,b){predict(y,cbind(100,a,b))})
#'    image(gx,gy,Y,asp=1,main="MMC clustering",xlab=colnames(x)[1],ylab=colnames(x)[2])
#'    points(x[,-1],pch=19+y)
mmc <- function(x,k=2L,N0=2L,LAMBDA=1,seeds=1:50,
                nrbmArgsSvm=list(maxCP=10L,MAX_ITER=100L),
                nrbmArgsMmc=list(maxCP=20L,MAX_ITER=300L),
                mc.cores=getOption("mc.cores",1L),...) {  
  nrbmArgsMmc$riskFun <- mmcLoss(x,k=k,...)
  nrbmArgsSvm$LAMBDA <- nrbmArgsMmc$LAMBDA <- k*LAMBDA
  
  mmcFit <- function(seed) {
    tryCatch({
      # select a starting point w0 by randomly selecting N0 samples in each cluster 
      # and train a multi-class SVM on them
      set.seed(seed)
      n0 <- min(nrow(x),k*N0)
      i0 <- sample(seq_len(nrow(x)),n0)
      y0 <- as.factor(sample(rep_len(seq_len(k),n0)))
      
      nrbmArgsSvm$riskFun <- softMarginVectorLoss(x[i0,],y0)
      w0 <- do.call(nrbm,nrbmArgsSvm)
      
      # run MMC solver starting at w0
      nrbmArgsMmc$w0 <- w0
      do.call(nrbm,nrbmArgsMmc)
    },error = function(e) {
      print(e)
      w <- "error";attr(w,"obj") <- NA
      w
    })
  }
  err <- mclapply(mc.cores=mc.cores,seeds,function(i) attr(mmcFit(i),"obj"))
  err <- simplify2array(err)
  names(err) <- seeds

  # find the minimum of all runs
  W <- mmcFit(seeds[which.min(err)])
  W <- matrix(W,ncol(x),k)
  rownames(W) <- colnames(x)
  
  L <- nrbmArgsMmc$riskFun(W)
  y <- max.col(attr(L,"Y"))
  attr(y,"W") <- W
  attr(y,"errors") <- err
  attr(y,"best.seed") <- seeds[which.min(err)]
  class(y) <- "mmc"
  return(y)
}



#' Predict class of new instances according to a mmc model
#' 
#' @export
#' @param object a mmc object
#' @param x a matrix similar to the dataset used, i.e. where rows are instances for
#'          which class must be predicted
#' @param ... unused, present to satisfy the generic predict() prototype 
#' @return a integer vector whose length match nrow(x) and containing the predicted 
#'         class for each of the given instances.
predict.mmc <- function(object,x,...) {
  max.col(x %*% attr(object,"W"))
}




