#' Partitioning algorithm of a set of variables (associated with mode 2) oh a  three-way array
#'
#' Each group of variables is associated with a one-rank PARAFAC model (comp x loading x weight).
#' Moreover, a Non Negativity (NN) constraint may be added to the model, so that the loading coefficients have positive values.
#' Return an object of class clv3w.
#'
#' @usage CLV3W_kmeans(X,K,mode.scale=0,NN=FALSE,init=10,cp.rand=5)
#'
#' @param X      : a three way array - variables of mode 2 will be clustered
#' @param K      : number of clusters
#' @param mode.scale : scaling parameter applied to X, by default  centering of X (for mode 2 x mode 3) is done. By default no scaling (mode.scale=0) \cr
#'             0 : no scaling only centering - the default \cr
#'             1 : scaling with standard deviation of  (mode 2 x mode 3) elements \cr
#'             2 : global scaling (each block i.e. each mode 2 slice will have the same inertia ) \cr
#'             3 : global scaling (each block i.e. each mode 3 slice will have the same inertia )
#' @param     NN : non Negativity constraint to be added on the loading coefficients. By default no constraint (NN=FALSE)   \cr
#'          TRUE : a non negativity constrained is applied on the loading coefficients to set them as positive values \cr
#'         FALSE : loading coefficients may be either positive or negative
#' @param   init : either the number of random starts i.e. partitions generated for the initialisation (By default init=10)
#' @param cp.rand : number of random starts associated with the one rank Candecomp/Parafac model (By default cp.rand=10)
#' @return \item{}{
#'          \itemize{
#'                \item {clusters}{ :  in line 1, the groups membership in the initial partition; in line 2 the final groups membership }
#'                \item {comp}{ : the latent components of the clusters associated with the first mode }
#'                \item {loading}{ : the vector of loadings  associated with the second mode by cluster }
#'                \item {weigth}{ : the vector of weights  associated with the third mode by cluster }
#'                \item {criterion}{ : vector of loss giving for each cluster the residual amount between the sub-array and its reconstitution associated with the cluster one rank PARAFAC model}
#'                \item {niter}{ : number of iterations of the partitioning alorithm}
#'          }}
#'  @return \item{param}{ contains the clustering parameters
#'          \itemize{
#'                \item {X}{ : the scaled three-way array}
#'          }}
#' @return call : call of the method

#' @seealso summary.clv3W, print.clv3W
#'
#' @author  Veronique Cariou, \email{veronique.cariou@oniris-nantes.fr}
#' @references Wilderjans, T. F., & Cariou, V. (2016). CLV3W: A clustering around latent variables approach to detect panel disagreement in three-way conventional sensory profiling data. Food quality and preference, 47, 45-53.
#' @references Cariou, V., & Wilderjans, T. F. (2018). Consumer segmentation in multi-attribute product evaluation by means of non-negatively constrained CLV3W. Food Quality and Preference, 67, 18-26.
#'
#'
#' @examples data(coffee)
#' ## Cluster Analysis of coffee sensory descriptors with block scaling
#' ## to set the assessors to the same footing
#' res.coffee <- CLV3W_kmeans(coffee,K=2,NN=TRUE,mode.scale=3,init=1,cp.rand=1)
#' summary(res.coffee)
#' get_partition(res.coffee)
#'
#' @export
#'

CLV3W_kmeans<-function(X,K,mode.scale=0,NN=FALSE,init=10,cp.rand=5)
{
  if (!is.array(X))
    stop("Non convenient X array")
  if (length(dim(X))!=3)
    stop("Non convenient X array : must be a three-way array")
  if (is.null(K))
    stop("The number of clusters K must be defined as a parameter")
  p <- dim(X)[2] # number of elements  (mode 2)
  n <- dim(X)[1] # number of elements (mode 1) # in sensory ususally products
  q <- dim(X)[3] # number of elements (mode 3) # in QDA assessors

  meanX<-plyr::aaply(.data=X,.margins=c(2,3),.fun=mean) # centering each descriptor for each assessor
    #new version
  X<-plyr::aaply(X,1,function(s) {s-meanX})
  if (is.null(dimnames(X))) {
    dimnames(X)[[1]] <- 1:n
    dimnames(X)[[2]] <- 1:p
    dimnames(X)[[3]] <- 1:q
  }
  else {
    if (is.null(dimnames(X)[[1]])) {
      dimnames(X)[[1]] <- 1:n
    }
    if (is.null(dimnames(X)[[2]])) {
      dimnames(X)[[1]] <- 1:p
    }
    if (is.null(dimnames(X)[[3]])) {
      dimnames(X)[[1]] <- 1:q
    }
  }
  if (mode.scale==1) {
    # each column mode2 x mode 3 will be standradized
    sdX<-plyr::aaply(.data=X,.margins=c(2,3),.fun=sd)
    X<-plyr::aaply(X,1,function(s) {s/sdX})
  }
  else {

    if (mode.scale==2) {
      # tables associated with mode 2 are set at the same footing that is to say the same inertia
      sst<-function(x) {return(x^2)}
      sstX<-plyr::aaply(.data=X,.margins=c(2,3),.fun=sst)
      localInertia<-apply(sstX,1,sum) #same inertia for all consumers
      globalInertia<-sum(localInertia)

      for (ip in 1:p) {
        X[,ip,]<-X[,ip,]*(sqrt(globalInertia)/sqrt(p*localInertia[ip]))
      }
    }
    else if (mode.scale==3) {
      # tables associated with mode 3 are set at the same footing that is to say the same inertia
      sst<-function(x) {return(x^2)}
      sstX<-plyr::aaply(.data=X,.margins=c(2,3),.fun=sst)
      localInertia<-apply(sstX,2,sum)
      globalInertia<-sum(localInertia)

      for (iq in 1:q) {
        X[,,iq]<-X[,,iq]*globalInertia/(q*localInertia[iq])
      }

    }
  }

  ## test init
  groupes_init <- NULL
  maxiter=100
  if (length(init)==1) {
    nrand=init
    groupes_init <- vector("list",nrand)
    for (nr in 1:nrand) {
      groupes_init[[nr]] <- sample.int(K, size = p, replace = TRUE)
      while (length(unique(groupes_init[[nr]]))!=K) {
        groupes_init[[nr]] <- sample.int(K, size = p, replace = TRUE)
      }
    }
  }
  else {
    if (is.vector(init) && (length(init)==p)) {
      nrand=1
      groupes_init <- vector("list",nrand)
      groupes_init[[1]] <- init
    }
    else
      stop("Wrong initialisation : either the number of random starts or an initial partition of mode2 items")
  }

  lcritere <- vector("numeric",length=nrand)
  clusters <- matrix(0,nrow=2,ncol=p)

    ######################################################################################
    #  Kmeans step
    ######################################################################################
    if (K > 1) {
      min_critere<-0
      for (nr in 1:nrand) {

        groupes_tmp <- groupes_init[[nr]]
        oldcrit=100
        for (i in 1:maxiter) {              # criterion, component and (possibly) loadings, in each cluster
          critere <-rep(0,K)
          comp<-matrix(0,nrow=n,ncol=K)
          weights<-matrix(0,nrow=q,ncol=K)
          loadings<-matrix(0,nrow=p,ncol=K)
          for (k in 1:K) {
            ind <- which(groupes_tmp == k)
            if (length(ind) > 0) {
              if (length(ind)>1) {
                Xnew<-X[,ind,]
                res<-CP1_MS(Xnew,NN)
                comp[,k]<- as.vector(res$u)
                weights[,k] <- as.vector(res$w)
                loadings[ind,k]<-as.vector(res$v)
                critere[k]<-  res$aloss
              }
              else {
                Xnew<-as.matrix(X[,ind,])
                res<-svd(Xnew)
                Xhat<-res$d[1]*res$u[,1]%o%res$v[,1]
                critere[k]<-sum(as.vector((Xnew-Xhat)^2)) #check
                comp[,k]<- as.vector(res$u[,1])
                weights[,k]<- as.vector(res$v[,1])
                loadings[ind,k]<-res$d[1]
              }
            }
            else {
              print('groupe vide')
            }
          }    # end of the loop k

          old_group<-groupes_tmp
          # reassign the variables to the groups
          L<-matrix(0,nrow=p,ncol=K)
          for (j in 1:p) {
            y<-as.vector(X[,j,])
            for (k in 1:K) {
              x<-as.vector(comp[,k]%o%weights[,k])
              beta<-cov(x,y)/var(x)
              if (NN & beta<0)
                beta <- 0
              L[j,k]<-sum((y-beta*x)^2)
            }
            minj<-which(L[j,]==min(L[j,]))
            if (length(minj)>1) {   ## randomization if same criterion for several clusters
              minj = sample(minj,1,replace=FALSE)
            }
            groupes_tmp[j] = minj
          }
          #checking empty clusters
          gp_idx<-setdiff(1:K,groupes_tmp)
          for (k in gp_idx) {
            set.minj<-apply(L,1,min)
            maxj<-which(set.minj==max(set.minj))
            L[maxj,]<-0
            groupes_tmp[maxj]<-k
          }
          #end checking empty clusters
          if (length(which((old_group == groupes_tmp) == FALSE, arr.ind = TRUE)) == 0)
            break

        }
        critere <-rep(0,K)
        comp<-matrix(0,nrow=n,ncol=K)
        weights<-matrix(0,nrow=q,ncol=K)
        loadings<-matrix(0,nrow=p,ncol=K)

        for (k in 1:K) {
          ind <- which(groupes_tmp == k)
          if (length(ind) > 0) {
            if (length(ind)>1) {
              Xnew<-X[,ind,]
              res<-CP1_MS(Xnew,NN)
              comp[,k]<- as.vector(res$u)
              weights[,k] <- as.vector(res$w)
              loadings[ind,k]<-as.vector(res$v)
              critere[k]<-  res$aloss
            }
            else {
              Xnew<-as.matrix(X[,ind,])
              res<-svd(Xnew)
              Xhat<-res$d[1]*res$u[,1]%o%res$v[,1]
              critere[k]<-sum(as.vector((Xnew-Xhat)^2)) #?ventuellement ? revoir
              comp[,k]<- as.vector(res$u[,1])
              weights[,k]<- as.vector(res$v[,1])
              loadings[ind,k]<-1

            }
          }
        }    # end of the loop k
        lcritere[nr]<-sum(critere)
        if ((nr==1) || (sum(min_critere)>sum(critere))) {
          clusters[1,] <- groupes_init[[nr]]
          clusters[2,]<-groupes_tmp
          min_comp<-comp
          min_loadings <- loadings
          min_weights<-weights
          min_critere<-  critere
          niter<-i
        }
      }
    }
    else {
      res<-CP1_MS(X,NN)
      K=1
      clusters[1,] <-rep(1,length=p)
      clusters[2,]<-rep(1,length=p)
      min_comp<- as.matrix(res$u)
      min_weights <- as.matrix(res$w)
      min_loadings <- as.matrix(res$v)
      min_critere <- res$aloss
      niter=1
    }
      #inertie.pct = T/max(T) * 100
      rownames(min_comp) <- dimnames(X)[[1]]
      colnames(min_comp) <- paste("Comp", c(1:K), sep = "")
      rownames(min_loadings) <- dimnames(X)[[2]]
      colnames(min_loadings) <- paste("Comp", c(1:K), sep = "")
      rownames(min_weights) <- dimnames(X)[[3]]
      colnames(min_weights) <- paste("Comp", c(1:K), sep = "")
      critere <-  min_critere
      dimnames(clusters)  <- list(c("init","last"),dimnames(X)[[2]] )
#      listcc = list(groups=clusters,comp=min_comp,loading=loadings, weight=min_weights,critere=min_critere,lcritere=lcritere,niter=niter,appel=match.call())
      listcc = list(clusters=clusters,comp=min_comp,loading=min_loadings, weight=min_weights,criterion=min_critere,niter=niter,call=match.call())
      clv3Wclt<- listcc
      class(clv3Wclt) <- "clv3w"
      clv3Wclt$param <- list(X=X)
      clv3Wclt$call   <- match.call()
      glob.call <- as.list(clv3Wclt$call)
      if (!"mode.scale" %in% names(glob.call))
        glob.call$mode.scale = mode.scale
      if (!"NN" %in% names(glob.call))
        glob.call$NN = NN
      clv3Wclt$call <- as.call(glob.call)
      return(clv3Wclt)
}
