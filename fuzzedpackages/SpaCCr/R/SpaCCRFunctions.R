#' Solve Spatial Convex Clustering problem for missing data
#'
#' @param X A subject (n) by variable (p) matrix; the data
#' @param w A vector of length p-1; weights for clustering
#' @param gamma A positive scalar; regularization parameter
#' @param nu A positive scalar; augmented Lagrangian paramter
#' @param verbose Logical; should messages be printed?
#' @param tol.base A small positive scalar; convergence tolerance for base SpaCC problem.
#' @param tol.miss A small positive scalar; convergence tolerance for missing data problem.
#' @param max.iter.base A positive integer; maximum number of iterations for base SpaCC problem
#' @param max.iter.miss A positive integer; maximum number of iterations for missing data problem
#' @param Uinit An n by p matrix; initial value for U
#' @param Vinit An n by p-1 matrix; initial value for V
#' @param Laminit An n by p-1 matrix; initial value for Lam
#' @return A list with elements U,V, and Lam
#' @export
#' @examples
#' library(dplyr)
#' library(tidyr)
#' data("methy")
#' methy <- methy[1:20,1:10]
#' Coordinates <- methy$Genomic_Coordinate
#' methy %>%
#'   tbl_df() %>%
#'   select(-Chromosome,-Genomic_Coordinate) %>%
#'   gather(Subject,Value,-ProbeID) %>%
#'   spread(ProbeID,Value) -> X
#' SubjectLabels <- X$Subject
#' X <- X[,-1] %>% as.matrix()
#' X[1:5,1:5]
#' nsubj <- nrow(X)
#' nprobes <- ncol(X)
#' nweights <- choose(nprobes,2)
#' diff.vals <- diff(Coordinates)
#' too.far <- diff.vals > 20000
#' sig = 1/5e3
#' w.values <- exp(-sig*diff.vals)
#' w.values[too.far] = 0
#'
#' verbose=TRUE
#' tol.base = 1e-4
#' tol.miss = 1e-4
#' max.iter.base=5000
#' max.iter.miss=500
#' bo <-t(scale(t(X),center=TRUE,scale=FALSE))
#' bo[is.na(bo)] <- mean(bo,na.rm=TRUE)
#' best.gam = 1
#' Sol <- SpaCC_Missing(t(scale(t(X),center=TRUE,scale=FALSE)),
#'                          w.values,
#'                          gamma = best.gam,
#'                          nu=1/nsubj,
#'                          verbose=TRUE,
#'                          tol.base=tol.base,
#'                          tol.miss=tol.miss,
#'                          max.iter.base=max.iter.base,
#'                          max.iter.miss=max.iter.miss,
#'                          bo,
#'                          t(diff(t(bo))),
#'                          t(diff(t(bo))))
SpaCC_Missing <- function(X,w,gamma,nu=1/nrow(X),verbose=FALSE,tol.base=1e-4,tol.miss=1e-4,max.iter.base=5000,max.iter.miss=500,Uinit,Vinit,Laminit) {
  n <- nrow(X)
  p <- ncol(X)
  miss.inds <- is.na(X)
  Unew <- Uinit
  Vnew <- Vinit
  Lamnew <- Laminit

  err = 1
  iter = 0
  while(err > tol.miss & (iter < max.iter.miss)) {
    iter = iter +1
    Uold <- Unew
    Vold <- Vnew
    Lamold <- Lamnew

    Tnew <- X
    Tnew[miss.inds] <- Unew[miss.inds]
    sol <- SpaCC(Tnew,w,gamma,nu,verbose,tol.base,max.iter.base,Uold,Vold,Lamold)
    Unew <- sol[[1]]
    Vnew <- sol[[2]]
    Lamnew <- sol[[3]]

    err <- norm(Unew - Uold,type = 'f')/norm(Uold,type='f')

  }
  if(verbose){
    print(paste("Miss Iter",iter))
    print(paste("Missing Error:",err))
  }
  list(U=Unew,V=Vnew,Lam=Lamnew)

}

#' Solve Spatial Convex Clustering problem for path of regularization parameters
#'
#' @param X A subject (n) by variable (p) matrix; the data
#' @param w A vector of length p-1; weights for clustering
#' @param gamma.seq A vector of positive scalars; regularization parameter sequence
#' @param nu A positive scalar; augmented Lagrangian paramter
#' @param verbose Logical; should messages be printed?
#' @param tol.base A small positive scalar; convergence tolerance for base SpaCC problem.
#' @param tol.miss A small positive scalar; convergence tolerance for missing data problem.
#' @param max.iter.base A positive integer; maximum number of iterations for base SpaCC problem
#' @param max.iter.miss A positive integer; maximum number of iterations for missing data problem
#' @return A list with elements UPath, VPath, LamPath, and gamma.seq
#' @export
#' @examples
#' NULL
SpaCC_Path <- function(X,w,gamma.seq,nu=1/nrow(X),verbose=FALSE,tol.base=1e-4,tol.miss=1e-4,max.iter.base=5000,max.iter.miss=500){
  n <- nrow(X)
  p <- ncol(X)
  gamma.seq <- sort(gamma.seq,decreasing = TRUE)
  ngam <- length(gamma.seq)

  bo <- X
  bo[is.na(bo)] = mean(bo,na.rm=TRUE)
  Uinit <- bo
  Vinit <- t(diff(t(bo)))
  Laminit <- Vinit

  UPath <- array(0,dim=c(n,p,ngam))
  VPath <- array(0,dim=c(n,p-1,ngam))
  LamPath <- array(0,dim=c(n,p-1,ngam))

  iter = 0
  for(gam in gamma.seq){
    iter = iter +1
    print(paste('gamma iter is:',iter))
    sol <- SpaCC_Missing(X,w,gam,nu,verbose,tol.base,tol.miss,max.iter.base,max.iter.miss,
                         Uinit,Vinit,Laminit)
    UPath[,,iter] <- sol$U
    VPath[,,iter] <- sol$V
    LamPath[,,iter] <- sol$Lam

  }
  list(UPath=UPath,VPath=VPath,LamPath=LamPath,gamma.seq=sort(gamma.seq,decreasing = TRUE))
}



#' Solve Spatial Convex Clustering problem for path of regularization parameters in parallel
#'
#' @param X A subject (n) by variable (p) matrix; the data
#' @param w A vector of length p-1; weights for clustering
#' @param gamma.seq A vector of positive scalars; regularization parameter sequence
#' @param nu A positive scalar; augmented Lagrangian paramter
#' @param verbose Logical; should messages be printed?
#' @param tol.base A small positive scalar; convergence tolerance for base SpaCC problem.
#' @param tol.miss A small positive scalar; convergence tolerance for missing data problem.
#' @param max.iter.base A positive integer; maximum number of iterations for base SpaCC problem
#' @param max.iter.miss A positive integer; maximum number of iterations for missing data problem
#' @param ncores A positive integer; number of cores to use
#' @return A list with elements UPath, VPath, LamPath, and gamma.seq
#' @export
#' @examples
#' NULL
SpaCC_Path_Parallel <- function(X,w,gamma.seq,nu=1/nrow(X),verbose=FALSE,tol.base=1e-4,tol.miss=1e-4,max.iter.base=5000,max.iter.miss=500,ncores=2){
  n <- nrow(X)
  p <- nrow(X)
  gamma.seq <- sort(gamma.seq,decreasing=TRUE)

  ChromBlocks <- GetParBlocks(X,w)
  BlockList <- ChromBlocks$Blocks
  wList <- ChromBlocks$weights
  wIndex <- ChromBlocks$weight.index
  nblocks <- length(BlockList)
  print(paste("Number of Blocks is", nblocks))

  ngam <- length(gamma.seq)
  if(nblocks == 1) {
    ret.me <- SpaCC_Path(X,wList[[1]], gamma.seq, nu, verbose=verbose,tol.base=tol.base,tol.miss=tol.miss,max.iter.base=max.iter.base,max.iter.miss=max.iter.miss)
    UPath <- ret.me$UPath
    VPath <- ret.me$VPath
    LamPath <- ret.me$VPath
  } else{
    weight.sum <- unlist(lapply(wList, sum))
    zero.logic <- weight.sum == 0
    isolated.probes <- lapply(wIndex[zero.logic], function(x){
      x[1]
    })
    full.path <- mcmapply(SpaCC_Path,
                          BlockList[!zero.logic],
                          wList[!zero.logic],
                          rep(list(gamma.seq),times=sum(!zero.logic)),
                          nu,
                          verbose,
                          tol.base,
                          tol.miss,
                          max.iter.base,
                          max.iter.miss,
                          mc.cores = ncores,
                          SIMPLIFY = FALSE,
                          mc.preschedule=FALSE)
    UList <- list()
    UList[!zero.logic] <- lapply(full.path, function(SpaCCPathResult){
      SpaCCPathResult$UPath
    })
    UList[zero.logic] <- lapply(isolated.probes,function(prb){
      na.inds <- is.na(X[,prb])
      mean.val <- mean(X[,prb],na.rm=TRUE)
      iso.prb <- X[,prb]
      iso.prb[na.inds] = mean.val
      ret.me <- array(iso.prb,dim=c(n,1,ngam))
      ret.me
    })
    UPath <- abind(UList,along=2)

    VList <- list()
    non.iso.ind = 1
    next.add = 1
    for(vlogic in seq_along(zero.logic)) {
      if(!zero.logic[vlogic]){
        VList[[next.add]] <- full.path[[non.iso.ind]][[2]]
        next.add = next.add +1
        non.iso.ind = non.iso.ind + 1
      }
      if(vlogic != length(zero.logic)){
        VList[[next.add]] <- array(rep(1,times=n),dim=c(n,1,ngam))
        next.add = next.add +1
      }
    }
    VPath <- abind(VList,along=2)

    LamList <- list()
    non.iso.ind = 1
    next.add = 1
    for(vlogic in seq_along(zero.logic)) {
      if(!zero.logic[vlogic]){
        LamList[[next.add]] <- full.path[[non.iso.ind]][[3]]
        next.add = next.add +1
        non.iso.ind = non.iso.ind + 1
      }
      if(vlogic != length(zero.logic)){
        LamList[[next.add]] <- array(rep(1,times=n),dim=c(n,1,ngam))
        next.add = next.add +1
      }
    }
    LamPath <- abind(LamList,along=2)
    rm(full.path)
    gc()


  }
  list(UPath = UPath,VPath = VPath,LamPath = LamPath,gamma.seq=gamma.seq)
}




#' Perform Cross Validation to select gamma/sparsity level
#'
#' @param X A subject (n) by variable (p) matrix; the data
#' @param w A vector of length p-1; weights for clustering
#' @param gamma.seq A vector of positive scalars; regularization parameter sequence
#' @param nfolds A positive scalar; number of cross validation folds
#' @param nu A positive scalar; augmented Lagrangian paramter
#' @param verbose Logical; should messages be printed?
#' @param tol.base A small positive scalar; convergence tolerance for base SpaCC problem.
#' @param tol.miss A small positive scalar; convergence tolerance for missing data problem.
#' @param max.iter.base A positive integer; maximum number of iterations for base SpaCC problem
#' @param max.iter.miss A positive integer; maximum number of iterations for missing data problem
#' @param parallel A logical; should CV paths be done in parallel?
#' @param frac A positive scalar between 0 and 1; fraction of hold out set to utilize
#' @return A list with elements: ErrMat - a length(gamma.seq) by nfold matrix containing
#' error on out of fold data; SpMat - a length(gamma.seq) by nfold matrix containing sparsity
#' levels; gamma.seq - original gamma.seq sorted largest to smallest
#' @export
#' @examples
#'library(dplyr)
#'library(tidyr)
#'data("methy")
#'methy <- methy[1:20,1:10]
#'Coordinates <- methy$Genomic_Coordinate
#'methy %>%
#'  tbl_df() %>%
#'  select(-Chromosome,-Genomic_Coordinate) %>%
#'  gather(Subject,Value,-ProbeID) %>%
#'  spread(ProbeID,Value) -> X
#'SubjectLabels <- X$Subject
#'X <- X[,-1] %>% as.matrix()
#'nsubj <- nrow(X)
#'nprobes <- ncol(X)
#'nweights <- choose(nprobes,2)
#'diff.vals <- diff(Coordinates)
#'too.far <- diff.vals > 20000
#'sig = 1/5e3
#'w.values <- exp(-sig*diff.vals)
#'w.values[too.far] = 0
#'
#'verbose=TRUE
#'tol.base = 1e-4
#'tol.miss = 1e-4
#'max.iter.base=5000
#'max.iter.miss=500
#'ngam = 20
#'gamma.seq <- exp(seq(log(1e-1),log(1e1),length.out=ngam))
#'CVRes <- SpaCC_CV(X=t(scale(t(X),center=TRUE,scale=FALSE)),
#'                  w=w.values,
#'                  gamma.seq=gamma.seq,
#'                  nfolds=5,
#'                  nu=1/nsubj,
#'                  verbose=TRUE,
#'                  tol.base=tol.base,
#'                  tol.miss=tol.miss,
#'                  max.iter.base=max.iter.base,
#'                  max.iter.miss=max.iter.miss,
#'                  parallel=FALSE,frac = 1)
SpaCC_CV <- function(X,w,gamma.seq,nfolds=5,nu=1/nrow(X),verbose=FALSE,tol.base=1e-4,tol.miss=1e-4,max.iter.base=5000,max.iter.miss=500,parallel=FALSE,frac=1) {
  Subject <- Value <- Cluster <- Probe <- MeanValue <- NULL
  if(parallel){
    Path_Function <- SpaCC_Path_Parallel
  } else {
    Path_Function <- SpaCC_Path
  }
  n <- nrow(X)
  p <- ncol(X)
  gamma.seq <- sort(gamma.seq,decreasing = TRUE)
  ngam <- length(gamma.seq)
  miss.inds <- is.na(X)
  non.miss.inds <- which(!miss.inds)
  n.non.miss <- length(non.miss.inds)
  folds <- sample(1:nfolds,size=n.non.miss,replace = TRUE)

  ErrMat <- matrix(0,nrow=ngam,ncol=nfolds)
  SpMat <- matrix(0,nrow=ngam,ncol=nfolds)
  for(fold in 1:nfolds){
    print(paste("Fold #",fold,"of",nfolds))
    X.fold <- X
    out.inds <- non.miss.inds[folds==fold]
    out.inds <- sample(out.inds,size=floor(frac*length(out.inds)),replace=FALSE)
    X.fold[out.inds] <- NA
    print(paste("Solving Path Problem",length(out.inds)))
    Path <- Path_Function(X=X.fold,w=w,gamma.seq=gamma.seq,nu=nu,verbose=verbose,tol.base=tol.base,tol.miss=tol.miss,max.iter.base=max.iter.base,max.iter.miss=max.iter.miss)
    print("Computing Error")
    for(gam.iter in seq_along(gamma.seq)){
      VThreshed <- Path$VPath[,,gam.iter]
      VThreshed[abs(VThreshed)<sqrt(log(p)/n)*sd(X.fold,na.rm=TRUE)] = 0
      clustsThreshed <- GetClusters(VThreshed)
      NEstRegion <- length(unique(clustsThreshed$cluster))
      NEstRegion
      SpMat[gam.iter,fold] = NEstRegion
      X.tmp <- X.fold
      X.tmp[is.na(X.fold)] <- Path$UPath[,,gam.iter][is.na(X.fold)]
      #t(X.fold) %>%
      t(X.tmp) %>%
        as.data.frame() %>%
        tbl_df() %>%
        mutate(
          Cluster = clustsThreshed$cluster) %>%
        gather(Subject,Value,-Cluster) %>%
        group_by(Subject,Cluster) %>%
        summarise(
          MeanValue = mean(Value,na.rm=TRUE)
        ) -> tmp1
      #t(X.fold) %>%
      t(X.tmp) %>%
        as.data.frame() %>%
        tbl_df() %>%
        mutate(
          Cluster = clustsThreshed$cluster,
          Probe = 1:ncol(X.fold)) %>%
        gather(Subject,Value,-Cluster,-Probe) -> tmp2
      tmp2 %>%
        left_join(tmp1,by=c('Subject','Cluster')) %>%
        select(-Value,-Cluster) %>%
        spread(Subject,MeanValue) %>%
        select(-Probe) %>%
        as.matrix() %>% t() -> tmp

      ErrMat[gam.iter,fold] = mean( (tmp[out.inds] - X[out.inds]) ^2,na.rm=TRUE)
    }
    rm(Path)
    gc()
  }
  list(ErrMat=ErrMat,gamma.seq=gamma.seq,SpMat=SpMat)

}

#' Function to compute blocks for parallization; should not be called directly
#'
#' @param X An n by p data matrix
#' @param w A vector of positive scalars of length p-1
#' @export
#' @examples
#' NULL
GetParBlocks <- function(X,w) {
  X.n <- nrow(X)
  X.p <- ncol(X)
  BlockList <- list()
  wList <- list()
  wIndex <- list()
  w.zinds <- which(w == 0)
  if(sum(w.zinds) == 0) {
    ret.me <- list(Blocks=list(X),
                   weights=list(w),
                   weight.index=list(1:length(w)))

  } else {
    start.ind = 1
    for(blk.ind in seq_along(w.zinds)) {
      z.ind <- w.zinds[blk.ind]

      BlockList[[blk.ind]] <- as.matrix(X[,start.ind:z.ind])
      wList[[blk.ind]] <- w[start.ind:(z.ind - 1)]
      wIndex[[blk.ind]] <- start.ind:(z.ind - 1)

      start.ind = z.ind + 1
    }
    BlockList[[blk.ind + 1]] = as.matrix(X[,start.ind:X.p])
    wList[[blk.ind + 1]] = w[start.ind:length(w)]
    wIndex[[blk.ind + 1]] = start.ind:length(w)

    wList[[blk.ind+1]][is.na(wList[[blk.ind+1]])] = 0
    ret.me <- list(Blocks=BlockList,
                   weights=wList,
                   weight.index =wIndex)

  }

  return(ret.me)
}

#' Plot subjects' copy number data with cluster means overlayed for a single chromosome
#'
#' @param Location A vector of length p with chromosomal locations
#' @param X A variable (p) by subject (n) data matrix
#' @param Cluster A vector of length p with cluster labels
#' @param NSubj A positve integer; number of randomly selected subjects to plot.
#' @param lowery,uppery Scalars; limits for y-axis
#' @importFrom stats sd
#' @importFrom utils head
#' @import parallel
#' @import abind
#' @import tidyr
#' @import dplyr
#' @import ggplot2
#' @export
#' @examples
#' NULL
CNVPlotSeriesMeans <- function(Location, X, Cluster,NSubj=3,lowery=-1,uppery=1) {
  Subject <- PrbMeasure <- MeanMeasure <- NULL
  X <- X[,sample(1:(ncol(X)),size=NSubj)]
  SimMatClust <- cbind.data.frame(Location, X, Cluster)
  NSubj <- ncol(X)
  names(SimMatClust) <-c(
    'Location',
    paste('Subject',1:NSubj, sep=''),
    'Cluster')

  rhs <- SimMatClust %>%
    gather(Subject,
           PrbMeasure,
           -Location,
           -Cluster) %>%
    group_by(Cluster, Subject) %>%
    summarise(
      MeanMeasure = mean(PrbMeasure,na.rm=TRUE),
      nProbe = n()
    )
  lhs <- SimMatClust %>%
    gather(Subject,
           PrbMeasure,
           -Location,
           -Cluster)
  inner_join(lhs,rhs,by=c('Subject','Cluster')) %>%
    ggplot(aes(x=Location,y=PrbMeasure)) +
    geom_point() +
    geom_line(aes(x=Location,y=MeanMeasure,
                  color=as.factor(Cluster)),size=I(2)) +
    facet_wrap(~Subject,nrow=NSubj) +
    scale_y_continuous(limits=c(lowery,uppery)) +
    guides(color=FALSE) +
    ylab('Copy Number Variation: Log2Ratio') +
    xlab('Genomic Location')
}

#' Compute Clusters from fusions
#'
#' @param V An n by p-1 data matrix
#' @export
#' @examples
#' NULL
GetClusters <- function(V){
  nprobes <- ncol(V) + 1
  clust.label <- 1
  clusters <- c()
  clusters[1] <- clust.label
  for(probe.iter  in 2:nprobes){
    if(sum(V[,probe.iter-1])!=0){
      clust.label = clust.label + 1
    }
    clusters[probe.iter] <- clust.label
  }
  clust.size <- unlist(lapply(unique(clusters),function(x){sum(clusters==x)}))
  list(cluster=clusters,size=clust.size)
}

#' A function for plotting cross validation errors
#'
#' @param ErrMat A matrix of error outputted by \code{SpaCC_CV}
#' @param rule An interger indicating which CV rule to choose
#' @param gamma.seq The sequence of regularization parameters
#' @export
#' @examples
#'library(dplyr)
#'library(tidyr)
#'data("methy")
#'methy <- methy[1:20,1:10]
#'Coordinates <- methy$Genomic_Coordinate
#'methy %>%
#'  tbl_df() %>%
#'  select(-Chromosome,-Genomic_Coordinate) %>%
#'  gather(Subject,Value,-ProbeID) %>%
#'  spread(ProbeID,Value) -> X
#'SubjectLabels <- X$Subject
#'X <- X[,-1] %>% as.matrix()
#'nsubj <- nrow(X)
#'nprobes <- ncol(X)
#'nweights <- choose(nprobes,2)
#'diff.vals <- diff(Coordinates)
#'too.far <- diff.vals > 20000
#'sig = 1/5e3
#'w.values <- exp(-sig*diff.vals)
#'w.values[too.far] = 0
#'
#'verbose=TRUE
#'tol.base = 1e-4
#'tol.miss = 1e-4
#'max.iter.base=5000
#'max.iter.miss=500
#'ngam = 20
#'gamma.seq <- exp(seq(log(1e-1),log(1e1),length.out=ngam))
#'CVRes <- SpaCC_CV(X=t(scale(t(X),center=TRUE,scale=FALSE)),
#'                  w=w.values,
#'                  gamma.seq=gamma.seq,
#'                  nfolds=5,
#'                  nu=1/nsubj,
#'                  verbose=TRUE,
#'                  tol.base=tol.base,
#'                  tol.miss=tol.miss,
#'                  max.iter.base=max.iter.base,
#'                  max.iter.miss=max.iter.miss,
#'                  parallel=FALSE,frac = 1)
#'PlotCV(CVRes$ErrMat,gamma.seq = CVRes$gamma.seq,rule = 1)
PlotCV <- function(ErrMat,rule=2,gamma.seq) {
  resp <- se <- NULL
  AvgErr <- apply(ErrMat,1,mean)
  SdErr <- apply(ErrMat,1,sd)
  min.err.ind <- which.min(AvgErr)
  if(rule == 1) {
    best.gam <- gamma.seq[min.err.ind]
  } else if(rule == 2){
    best.gam <- gamma.seq[min(which(AvgErr <= AvgErr[which.min(AvgErr)] + SdErr[which.min(AvgErr)]))]
  } else if (rule == 3){
    best.gam <- gamma.seq[min(which(AvgErr <= AvgErr[which.min(AvgErr)] + SdErr[which.min(AvgErr)])) - 1]
  }
  else{
    stop('Unknown Rule')
  }

  ErrDf <- data.frame(
    resp = AvgErr,
    se = SdErr
  )
  limits <- aes(ymax = resp + se, ymin = resp - se)
  qplot(gamma.seq,AvgErr,data=ErrDf) + geom_errorbar(limits,width=0.2) + scale_x_log10() + geom_vline(xintercept=best.gam)

}

#' Get optimal cross validated gamma values by various rules
#' @param ErrMat Matrix of cross validated errors outputted by \code{GetCVErrMat}
#' @param rule A number indicating how optimal gamma should be chosen. 1 for minimum cv error,
#' 2 for 1 standard error rule
#' @param gamma.seq sequence of regularization parameters used for cross validation.
#' @return A scalar. Optimal gamma selected by CV rule.
#' @export
#' @examples
#' NULL
GetGammaCV <- function(ErrMat,rule=1,gamma.seq) {
  gamma.seq <- sort(gamma.seq,decreasing=TRUE)
  AvgErr <- apply(ErrMat,1,mean)
  SdErr <- apply(ErrMat,1,sd)
  min.err.ind <- which.min(AvgErr)
  if(rule == 1) {
    best.gam <- gamma.seq[min.err.ind]
  } else if(rule ==2) {
    best.gam <- gamma.seq[min(which(AvgErr <= AvgErr[which.min(AvgErr)] + SdErr[which.min(AvgErr)]))]
  } else if (rule ==3){
    best.gam <- gamma.seq[min(which(AvgErr <= AvgErr[which.min(AvgErr)] + SdErr[which.min(AvgErr)])) - 1]
  }else{
    stop('Unknown Rule')
  }
  return(best.gam)
}


#' Threshold differences
#' @param V an n x p-1 matrix of differences
#' @param X an n x p matrix
#' @param mult scalar to multiply standard deviation
#' @param thresh.value optional user specified threshold value.
#' @return VThreshed an n x p-1 matrix of thresholded differences
#' @export
#' @examples
#' NULL
ThreshV <- function(V,X,mult=1,thresh.value=NULL){
  if(is.null(thresh.value)){
    p <- ncol(V)+1
    n <- nrow(V)
    thresh.value <- mult*sd(X,na.rm=TRUE)*sqrt(log(p)/n)
  }
  ret <- V
  ret[abs(V)<=thresh.value]=0
  ret
}

#' Plots methylation data by Genomic Coordinates for a given chromosomal region with cluster means overlayed for each subject.
#'
#' @param X A Subject by Probe data matrix for a single chromosome of CNV data
#' @param Coord A vector of Genomic Coordinates for a single chromosome
#' @param Cluster Cluster labels for each probe
#' @param SubjInd A vector of numeric indicies corresponding to the Subjects to be plotted.
#' @param Start Genomic Coordinate minimum
#' @param End Genomic Coordinate maximum
#' @export
#' @examples
#'library(dplyr)
#'library(tidyr)
#'data("methy")
#'methy <- methy[1:20,1:10]
#'Coordinates <- methy$Genomic_Coordinate
#'methy %>%
#'  tbl_df() %>%
#'  select(-Chromosome,-Genomic_Coordinate) %>%
#'  gather(Subject,Value,-ProbeID) %>%
#'  spread(ProbeID,Value) -> X
#'SubjectLabels <- X$Subject
#'X <- X[,-1] %>% as.matrix()
#'nsubj <- nrow(X)
#'nprobes <- ncol(X)
#'nweights <- choose(nprobes,2)
#'diff.vals <- diff(Coordinates)
#'too.far <- diff.vals > 20000
#'sig = 1/5e3
#'w.values <- exp(-sig*diff.vals)
#'w.values[too.far] = 0
#'
#'verbose=TRUE
#'tol.base = 1e-4
#'tol.miss = 1e-4
#'max.iter.base=5000
#'max.iter.miss=500
#'ngam = 20
#'gamma.seq <- exp(seq(log(1e-1),log(1e1),length.out=ngam))
#'CVRes <- SpaCC_CV(X=t(scale(t(X),center=TRUE,scale=FALSE)),
#'                  w=w.values,
#'                  gamma.seq=gamma.seq,
#'                  nfolds=5,
#'                  nu=1/nsubj,
#'                  verbose=TRUE,
#'                  tol.base=tol.base,
#'                  tol.miss=tol.miss,
#'                  max.iter.base=max.iter.base,
#'                  max.iter.miss=max.iter.miss,
#'                  parallel=FALSE,frac = .1)
#'PlotCV(CVRes$ErrMat,gamma.seq = CVRes$gamma.seq,rule = 1)
#'best.gam <- GetGammaCV(CVRes$ErrMat,rule = 1,gamma.seq = CVRes$gamma.seq)
#'bo <-t(scale(t(X),center=TRUE,scale=FALSE))
#'bo[is.na(bo)] <- mean(bo,na.rm=TRUE)
#'Sol <- SpaCC_Missing(t(scale(t(X),center=TRUE,scale=FALSE)),
#'                         w.values,
#'                         gamma = best.gam,
#'                         nu=1/nsubj,
#'                         verbose=TRUE,
#'                         tol.base=tol.base,
#'                         tol.miss=tol.miss,
#'                         max.iter.base=max.iter.base,
#'                         max.iter.miss=max.iter.miss,
#'                         bo,
#'                         t(diff(t(bo))),
#'                         t(diff(t(bo))))
#'VThreshed <- Sol$V
#'clustsThreshed <- GetClusters(VThreshed)
#'NEstRegion <- length(unique(clustsThreshed$cluster))
#'NEstRegion
#'VThreshed <- ThreshV(Sol$V,X,mult = 1)
#'clustsThreshed <- GetClusters(VThreshed)
#'NEstRegion <- length(unique(clustsThreshed$cluster))
#'NEstRegion
#'start.coord <- 2e5
#'end.coord <- 4e5
#'MethyRegionPlot(X,Coordinates,clustsThreshed$cluster,SubjInd = 1:3,Start=start.coord,End=end.coord)
MethyRegionPlot <- function(X,Coord,Cluster,SubjInd=1:3,Start,End){
  Subject <- PrbMeasure <- Location <- MeanMeasure <- NULL
  X <- t(X[SubjInd,])
  NSubj <- length(SubjInd)
  SimMatClust <- cbind.data.frame(Coord, X, Cluster)
  names(SimMatClust) <-c(
    'Location',
    paste('Subject',1:NSubj, sep=''),
    'Cluster')

  rhs <- SimMatClust %>%
    gather(Subject,
           PrbMeasure,
           -Location,
           -Cluster) %>%
    group_by(Cluster, Subject) %>%
    summarise(
      MeanMeasure = mean(PrbMeasure),
      nProbe = n()
    )
  lhs <- SimMatClust %>%
    gather(Subject,
           PrbMeasure,
           -Location,
           -Cluster)
  inner_join(lhs,rhs,by=c('Subject','Cluster')) %>%
    ggplot(aes(x=Location,y=PrbMeasure)) +
    geom_point() +
    geom_line(aes(x=Location,y=MeanMeasure,
                  color=as.factor(Cluster))) +
    facet_wrap(~Subject,nrow=NSubj) +
    scale_y_continuous(limits=c(-1.1,1.1)) +
    scale_x_continuous(limits=c(Start,End))+
    guides(color=FALSE)
}


#' Performs Spatial Convex Clustering for methylation data
#' @param X A subject (n) by variable (p) matrix; the data
#' @param Coordinates a vector listing genomic coordinates
#' @param gamma.seq a vector of regularization parameters
#' @param dist.cutoff maximum distance at which probes should be regularized
#' @param sig positive scalar controling spatial weight decay
#' @param weights a vector of spatial weights
#' @param center should data be centered
#' @param scale should data be scaled
#' @param nfolds number of folds for cross validation
#' @param nu parameter for augmented lagrangian
#' @param tol.base tolerance level for base function
#' @param tol.miss tolerance for missing function
#' @param max.iter.base maximum number of iterations for base function
#' @param max.iter.miss maximum number of iterations for missing function
#' @param frac fration of fold to use for cross validation
#' @param parallel should algorithm be run in parallel
#' @param gam.rule cross validation rule
#' @param thresh.mult multiplier for threshold value
#' @param thresh.value value of threshold
#' @return Labels a vector of cluster labels
#' @export
#' @examples
#'data("methy")
#'methy <- methy[1:20,1:10]
#'library(dplyr)
#'library(tidyr)
#'Coordinates <- methy$Genomic_Coordinate
#'methy %>%
#'  tbl_df() %>%
#'  select(-Chromosome,-Genomic_Coordinate) %>%
#'  gather(Subject,Value,-ProbeID) %>%
#'  spread(ProbeID,Value) -> X
#'SubjectLabels <- X$Subject
#'X <- X[,-1] %>% as.matrix()
#'verbose=TRUE
#'tol.base = 1e-4
#'tol.miss = 1e-4
#'max.iter.base=5000
#'max.iter.miss=500
#'ngam = 20
#'gamma.seq <- exp(seq(log(1e-1),log(1e1),length.out=ngam))
#'ClusterLabels <- SpaCC_Methy(X = X,Coordinates = Coordinates,gamma.seq = gamma.seq)
SpaCC_Methy <- function(X,
                        Coordinates,
                        gamma.seq,
                        dist.cutoff = 20000,
                        sig=1/5e3,
                        weights=NULL,
                        center=TRUE,
                        scale=FALSE,
                        nfolds=5,
                        nu=NULL,
                        tol.base=1e-4,
                        tol.miss=1e-4,
                        max.iter.base=5000,
                        max.iter.miss=500,
                        frac=.1,
                        parallel=FALSE,
                        gam.rule=2,
                        thresh.mult=1,
                        thresh.value=NULL){
  nsubj <- nrow(X)
  nprobes <- ncol(X)
  ngam <- length(gamma.seq)
  if(is.null(weights)){
    nweights <- choose(nprobes,2)
    diff.vals <- diff(Coordinates)
    too.far <- diff.vals > dist.cutoff
    w.values <- exp(-sig*diff.vals)
    w.values[too.far] = 0
  }
  if(is.null(nu)){
    nu <- 1/nsubj
  }
  CVRes <- SpaCC_CV(X=t(scale(t(X),center=center,scale=scale)),
                    w=w.values,
                    gamma.seq=gamma.seq,
                    nfolds=nfolds,
                    nu=nu,
                    verbose=TRUE,
                    tol.base=tol.base,
                    tol.miss=tol.miss,
                    max.iter.base=max.iter.base,
                    max.iter.miss=max.iter.miss,
                    parallel=parallel,
                    frac=frac)
  best.gam <- GetGammaCV(CVRes$ErrMat,rule = gam.rule,gamma.seq = CVRes$gamma.seq)
  bo <-t(scale(t(X),center=TRUE,scale=FALSE))
  bo[is.na(bo)] <- mean(bo,na.rm=TRUE)
  Sol <- SpaCC_Missing(t(scale(t(X),center=TRUE,scale=FALSE)),
                           w.values,
                           gamma = best.gam,
                           nu=1/nsubj,
                           verbose=TRUE,
                           tol.base=tol.base,
                           tol.miss=tol.miss,
                           max.iter.base=max.iter.base,
                           max.iter.miss=max.iter.miss,
                           bo,
                           t(diff(t(bo))),
                           t(diff(t(bo))))
  if(is.null(thresh.value)){
    VThreshed <- ThreshV(Sol$V,X,mult = thresh.mult)
  } else{
    VThreshed <- ThreshV(Sol$V,X,mult = thresh.mult,thresh.value = thresh.value)
  }
  clustsThreshed <- GetClusters(VThreshed)
  NEstRegion <- length(unique(clustsThreshed$cluster))
  clustsThreshed$cluster
}
