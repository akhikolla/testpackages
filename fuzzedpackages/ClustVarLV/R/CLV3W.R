#' Hierarchical clustering of variables (associated with mode 2 three-way array)  with consolidation
#'
#' Hierarchical Cluster Analysis of a set of variables (mode 2) given a three-way array with a further consolidation step.
#' Each group of variables is associated with a one-rank PARAFAC model (comp x loading x weight).
#' Moreover, a Non Negativity (NN) constraint may be added to the model, so that the loading coefficients have positive values.
#' Return an object of class clv3w.
#'
#' @usage CLV3W(X,mode.scale=0,NN=FALSE,moddendoinertie=TRUE,gmax=20,graph=TRUE,cp.rand=10)
#'
#' @param X      : a three way array - variables of mode 2 will be clustered
#' @param mode.scale : scaling parameter applied to X, by default  centering of X (for mode 2 x mode 3) is done. By default no scaling (mode.scale=0) \cr
#'             0 : no scaling only centering - the default \cr
#'             1 : scaling with standard deviation of  (mode 2 x mode 3) elements \cr
#'             2 : global scaling (each block i.e. each mode 2 slice will have the same inertia ) \cr
#'             3 : global scaling (each block i.e. each mode 3 slice will have the same inertia )
#' @param     NN : non Negativity constraint to be added on the loading coefficients. By default no constraint (NN=FALSE)   \cr
#'          TRUE : a non negativity constrained is applied on the loading coefficients to set them as positive values \cr
#'         FALSE : loading coefficients may be either positive or negative
#' @param   moddendoinertie : dendrogram. By default it is based on the delta clustering criterion (moddendoinertie =TRUE) \cr
#'          TRUE : dendrogram associated with the clustering criterion delta \cr
#'         FALSE : dendrogram associated with the the height (cumulative delta)
#' @param   gmax : maximum number of partitions for which the consolidation will be done (default : gmax=11)
#' @param  graph : boolean, if TRUE, the graphs associated with the dendrogram and the evolution of the aggregation criterion are displayed (default : graph=TRUE)
#' @param cp.rand : number of random starts associated with the one rank Candecomp/Parafac model (By default cp.rand=10)
#' @return \item{tabres}{ Results of the hierarchical clustering algorithm.
#'         In each line you find the results of one specific step of the hierarchical clustering.
#'         \itemize{
#'                \item {Columns 1 and 2}{ : the numbers of the two groups which are merged}
#'                \item {Column 3}{ : name of the new cluster}
#'                \item {Column 4}{ : the value of the aggregation criterion for the Hierarchical Ascendant Clustering (delta) : delta loss}
#'                \item {Column 5}{ : the loss value of the clustering criterion for the HAC}
#'                \item {Column 6}{ : the percentage of explained inertia of the data array X}
#'                \item {Column 7}{ : the loss value of the clustering criterion  after consolidation}
#'                \item {Column 8}{ : the percentage of explained inertia of the data array X after consolidation }
#'                \item {Column 9}{ : number of iterations in the partitioning algorithm. \cr
#'                Remark : A zero in columns 7 to 9 indicates that no consolidation was done}
#'        }}
#' @return \item{hclust}{ contains the results of the HCA }
#' @return \item{partition K}{ contains a list for each number of clusters of the partition, K=1 to gmax with
#'          \itemize{
#'                \item {clusters}{ :  in line 1, the groups membership before consolidation; in line 2 the groups membership after consolidation}
#'                \item {comp}{ : the latent components of the clusters associated with the first mode (after consolidation)}
#'                \item {loading}{ : the vector of loadings  associated with the second mode by cluster (after consolidation)}
#'                \item {weigth}{ : the vector of weights  associated with the third mode by cluster (after consolidation)}
#'                \item {criterion}{ : vector of loss giving for each cluster the residual amount between the sub-array and its reconstitution associated with the cluster one rank PARAFAC model (after consolidation)}
#'          }}
#' @return \item{param}{ contains the clustering parameters
#'          \itemize{
#'                \item {gmax}{ :  maximum number of partitions for which the consolidation has been done}
#'                \item {X}{ : the scaled three-way array}
#'          }}
#' @return call : call of the method
#' @seealso CLV3W_kmeans, get_comp, get_loading, get_partition, plot, plot_var.clv3w,
#'
#' @author  Veronique Cariou, \email{veronique.cariou@oniris-nantes.fr}
#' @references Wilderjans, T. F., & Cariou, V. (2016). CLV3W: A clustering around latent variables approach to detect panel disagreement in three-way conventional sensory profiling data. Food quality and preference, 47, 45-53.
#' @references Cariou, V., & Wilderjans, T. F. (2018). Consumer segmentation in multi-attribute product evaluation by means of non-negatively constrained CLV3W. Food Quality and Preference, 67, 18-26.
#'
#'
#' @examples data(ciders)
#' ## Cluster Analysis of cider sensory descriptors with block scaling
#' ## to set the assessors to the same footing
#' res.cider<-CLV3W(ciders,mode.scale=3,NN=FALSE,moddendoinertie=FALSE,gmax=20,graph=FALSE,cp.rand=5)
#' plot(res.cider,type="delta")
#' plot(res.cider,type="dendrogram")
#' print(res.cider)
#' summary(res.cider,2)
#' get_comp(res.cider,2)
#' get_loading(res.cider,2)
#' get_weight(res.cider,2)
#'
#'
#' @importFrom stats sd cov var as.dendrogram runif
#' @importFrom grDevices dev.new
#' @importFrom graphics plot title barplot
#' @importFrom plyr aaply
#'
#' @export
#'



 CLV3W<-function(X,mode.scale=0,NN=FALSE,moddendoinertie=TRUE,gmax=20,graph=TRUE,cp.rand=10)
{
 p <- dim(X)[2] # number of elements  (mode 2)
 n <- dim(X)[1] # number of elements (mode 1) # in sensory ususally products
 q <- dim(X)[3] # number of elements (mode 3) # in QDA assessors

 #centering mode2 x mode3
 meanX<-plyr::aaply(.data=X,.margins=c(2,3),.fun=mean) # centering each descriptor for each assessor
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

  gmax <- min(p,gmax)
  gmax <- max(2,gmax)
  globalInertia <- sum(X^2)
  stotal<-globalInertia

#  ###############################################################################
#  ## Initialisations
#  ###############################################################################

 #groupes : vector identifying the number of groups at each step
 # ex 5 var and 3 groups       groupes=c(1,2,1,3,3)
 groupes <- 1:p

 #fusions : vector which gives the variable number, (negative), if the variable forms one group
 #          alternatively it indicates the level of agregation
 fusions <- -groupes
 #initialisation of the criterion
 crit<- vector("numeric",length=p)

 for (ij in 1:p) {
   res<-svd(X[,ij,])
   Xhat<-res$d[1]*res$u[,1]%o%res$v[,1]
   crit[ij]<-sum((X[,ij,]-Xhat)^2) #Error criterion  # voir eigen
 }

 #inertie : clustering criterion at each step
 fit <- sum(crit)
 sbegin <- sum(crit)
 inertie <- sum(crit)
 ncluster <- p        # number of clusters


#  storing
  results = matrix(0,p-1,9)
  resultscc <- list()


 # delta : aggregation criterion at each step
 # hmerge : history of the aggregation  (each individual var is associated to its negative number)
 critstep <-  matrix(nrow=p, ncol=p)
 deltamin<-matrix(nrow=p,ncol=p)
 hmerge <-  matrix(0,nrow=p-1, ncol=2)
 delta <- matrix(0,nrow=p-1, ncol=1)
 ordr <- c()

    pg <- max(groupes) # number of groups
    for (i in 1:(pg - 1)) {
        for (j in (i + 1):pg) {
          Xnew<- X[, c(i, j),]
          res<-CP1_MS(Xnew,NN,cp.rand)
          critstep[i,j] = res$aloss
          # variation in the clustering criterion
          deltamin[i,j] <-  -crit[i] - crit[j] + critstep[i,j]
        }   # end  j
    }       # end  i


 ###############################################################################
 # loops : level
 ###############################################################################

 for (level in 1:(p-1)) {


    ##################################################################
    # Hierarchical step
    ##################################################################

    # For which couple of clusters is Delta minimum?
		cmerge<-which(deltamin==min(deltamin,na.rm=TRUE),arr.ind=TRUE)

		# Is the minimum value unique?
		if (nrow(cmerge)>1) cmerge<-cmerge[1,]

		ind1 <- which(groupes == cmerge[1])
    ind2 <- which(groupes == cmerge[2])
    ind<-c(fusions[ind1[1]],fusions[ind2[1]])
    hmerge[level,]<-ind
    listind<-c(ind1,ind2)

    # for the dendrogram : order of the var
    ordr <- c(ordr,listind);
    dupli<-duplicated(ordr,fromLast=FALSE)
    duplibis<-duplicated(ordr,fromLast=TRUE)
    if (sum(dupli)>0) {
      idem<-which(dupli)
      dupli[(length(ordr)-length(listind)+1):(length(ordr))]<-as.logical(1-dupli[(length(ordr)-length(listind)+1):(length(ordr))])
      dif<-which(dupli)
      if(length(dif)>0) {
        num<-ordr[dif]
        ou<-which(duplibis)
        tempo=matrix(0,nrow=1,ncol=length(ordr)-length(idem))
        tempo[1:max(ou)]<-ordr[1:max(ou)]
        tempo[max(ou)+1]<-num
        if (length(tempo)>(max(ou)+1))  tempo[(max(ou)+2):(length(ordr)-length(idem))]<-ordr[(max(ou)+1):(length(ordr)-length(listind))]
        ordr<-tempo
      }
      else {
        ordr<-ordr[-((length(ordr)-length(listind)+1):length(ordr))]
      }
    }

    ######## renewal of the parameters #########################################
    fusions[listind]<- level
    crit[listind[1]]<- critstep[cmerge[1],cmerge[2]];
    autre <- listind[2:length(listind)]
    crit[autre]<-NA

    inertie <- rbind(inertie, (inertie[level] + deltamin[cmerge[1],cmerge[2]]))     #inertie=sum(crit, na.rm = TRUE) also
	  delta[level]<-deltamin[cmerge[1],cmerge[2]]
		inertie_inv<- -inertie+sbegin
	  ncluster <- ncluster-1

    results[level,1:6]<- c(ind,level,delta[level],inertie[level+1],100*(1-inertie[level+1]/stotal))

	  # All variables of the cluster cmerge[2] join the cluster cmerge[1]:
    groupes[groupes==cmerge[2]]<-cmerge[1]
    # The cluster numbers above cmerge[2] are renamed:
    groupes[groupes>cmerge[2]]<-groupes[groupes>cmerge[2]]-1
    # Update the matrix deltamin and matrix critstep:
		# Remove row and column of deltamin and critstep corresponding to the second merged cluster:
    deltamin<-deltamin[-cmerge[2],-cmerge[2]]
    critstep<-critstep[-cmerge[2],-cmerge[2]]


    # Update the values concerning the new cluster:
    gr2<-which(groupes==cmerge[1])
    if (cmerge[1]>1)	{
			for (iter in 1:(cmerge[1]-1)){
				gr1<-which(groupes==iter)
        Xnew<-X[,c(gr1,gr2),]
  			res<-CP1_MS(Xnew,NN,cp.rand)
  			critstep[iter,cmerge[1]] = res$aloss
       deltamin[iter,cmerge[1]]<--crit[gr1[1]]-crit[gr2[1]]+critstep[iter,cmerge[1]]
			}
    }
    if (cmerge[1]<max(groupes))	 {
			for (iter in (cmerge[1]+1):max(groupes))	{
				gr1<-which(groupes==iter)
				Xnew<-X[,c(gr2,gr1),]
				res<-CP1_MS(Xnew,NN,cp.rand)
				critstep[cmerge[1],iter] = res$aloss
		    deltamin[cmerge[1],iter]<--crit[gr2[1]]-crit[gr1[1]]+critstep[cmerge[1],iter]
			}
    }


    ######################################################################################
    #  Consolidation, if the number of clusters has attained gmax
    ######################################################################################
    if (ncluster <= gmax & ncluster > 1 ) {
    cc_consol <- t(t(groupes))
    K <- ncluster
    T<-c()
    maxiter=100
    oldcrit=100

    for (i in 1:maxiter) {              # criterion, component and (possibly) loadings, in each cluster
        critere <-rep(0,K)
        comp<-matrix(0,nrow=n,ncol=K)
        weights<-matrix(0,nrow=q,ncol=K)
        groupes_tmp <- cc_consol[,i]
        for (k in 1:K) {
          ind <- which(groupes_tmp == k)
          if (length(ind) > 0) {
            if (length(ind)>1) {
              Xnew<-X[,ind,]
              res<-CP1_MS(Xnew,NN)
              comp[,k]<- as.vector(res$u)
              weights[,k] <- as.vector(res$w)
              critere[k]<-  res$aloss
            }
            else {
              Xnew<-as.matrix(X[,ind,])
              res<-svd(Xnew)
              Xhat<-res$d[1]*res$u[,1]%o%res$v[,1]
              critere[k]<-sum((Xnew-Xhat)^2) #criterion update
              comp[,k]<- as.vector(res$u[,1])
              weights[,k]<- as.vector(res$v[,1])
            }
          }
          else {
            print(c('groupe vide',k))
          }
        } # end of the loop k

        T = cbind(T, sum(critere))
        # re-aassigning the variables to the groups
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
              groupes_tmp[j] = sample(minj,1,replace=FALSE)
          }
          else {
            groupes_tmp[j] = minj
          }
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
        if ((length(which((cc_consol[, i] == groupes_tmp) == FALSE, arr.ind = TRUE)) == 0) | (abs(oldcrit-sum(critere))<1e-06)) {
          break
        }
        cc_consol = cbind(cc_consol, groupes_tmp)
        oldcrit<-sum(critere)
    }
    #a last step to the final computation of the PARAFAC latent components associated to each group
      niter<-ncol(cc_consol)
      critere <-rep(0,K)
      comp<-matrix(0,nrow=n,ncol=K)
      weights<-matrix(0,nrow=q,ncol=K)
      loadings<-matrix(0,nrow=p,ncol=K)
      groupes_tmp <- cc_consol[,niter]
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
            if (all(weights[,k]<0)) {
              weights[,k]<--weights[,k]
              if (NN==TRUE) {
                comp[,k]<- - comp[,k]
              } else {
                loadings[,k]<- - loadings[,k]
              }
            }
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
      }    # end of the loop k
    #inertie.pct = T/max(T) * 100
    rownames(comp) <- dimnames(X)[[1]]
    colnames(comp) <- paste("Comp", c(1:K), sep = "")
    rownames(loadings) <- dimnames(X)[[2]]
    colnames(loadings) <- paste("Comp", c(1:K), sep = "")
    rownames(weights) <- dimnames(X)[[3]]
    colnames(weights) <- paste("Comp", c(1:K), sep = "")
    rownames(cc_consol) <- dimnames(X)[[2]]
    names(cc_consol) = NULL
    initgroupes<-cc_consol[,1]
    lastgroupes<-cc_consol[,ncol(cc_consol)]
    listcc = list(clusters = rbind(initgroupes,lastgroupes),  comp=comp,weight=weights,loading=loadings,criterion=critere)
    results[level,7:9]<- c(sum(critere),100*(1-sum(critere)/stotal), i)
    # composante[[K]]=comp
    # idgr[[K]]= cc_consol
     resultscc[[K]] <- listcc

    }  # end of the loop "if ncluster <= gmax"
 }      # end of loop level

  #the results when the variables form only one cluster
  results[p-1,7:9]<- c(results[p-1,5],results[p-1,6], 0)
  groupes<-as.matrix(rep(1,p))
  rownames(groupes) <- dimnames(X)[[2]]
  res<-CP1_MS(X,NN)
  comp<- as.matrix(res$u)
  weights <- as.matrix(res$w)
  loadings <- as.matrix(res$v)
  rownames(comp) <- dimnames(X)[[1]]
  colnames(comp) <- "Comp1"

  rownames(loadings) <- dimnames(X)[[2]]
  colnames(loadings) <- "Comp1"
  rownames(weights) <- dimnames(X)[[3]]
  colnames(weights) <- "Comp1"
  list = list(clusters = rbind(t(groupes),t(groupes)),  comp=comp, loading=loadings,weight=weights,criterion=res$aloss)

  resultscc[[1]] <- list

 # the results when every variables is a cluster itself
   if (gmax == p) {
        groupes<-1:p
        comp<-matrix(0,nrow=n,ncol=p)
        weights<-matrix(0,nrow=q,ncol=p)
        loadings <- diag(p)
        critere <- vector("numeric",length=p)
        for (k in 1:p) {
            ind <- which(groupes == k)
            res<-svd(X[,ind,])
            Xhat<-res$d[1]*res$u[,1]%o%res$v[,1]
            comp[,k]<-res$u[,1]
            weights[,k]<-res$v[,1]
           critere[k]<-sum(as.vector((X[,ind,]-Xhat)^2)) #check
        }
        rownames(comp) <-dimnames(X)[[1]]
        colnames(comp) <- paste("Comp", c(1:p), sep = "")
        rownames(loadings) <- dimnames(X)[[2]]
        colnames(loadings) <- paste("Comp", c(1:p), sep = "")
        rownames(weights) <- dimnames(X)[[3]]
        colnames(weights) <- paste("Comp", c(1:p), sep = "")

        list = list(clusters = matrix(rep(t(groupes),2),nrow=2,ncol=p,byrow=TRUE),  comp=comp,weight=weights,loading=loadings,criterion=critere)
        resultscc[[gmax]] <- list
   }


 ################################################################################
  # Adding a column for the delta criterion
  cc.consol <- min(which(results[,7]!=0))
  if (length(cc.consol)==0) {
    loss.crit.cc <- cbind(results[1:cc.consol,5],results[(cc.consol+1),7])
  } else {
    loss.crit.cc <- results[,7]
  }
  colnames(results)= c("merg1","merg2","new.clust","delta.crit.hac","loss.crit.hac","%S0expl.hac","loss.crit.cc","%S0expl.cc","iter")
  names(resultscc) = paste("partition",1:gmax,sep="")
  clv3Wclt=resultscc
  clv3Wclt$tabres<- results
#  clv3Wclt$globalInertia<-globalInertia
#  clv3Wclt$sbegin<-sbegin


 ###############################################################################


      if (moddendoinertie==FALSE){ ## le faire en delta cumsum
        resultscah=list(labels=dimnames(X)[[2]], inertie=inertie, height=cumsum(results[,4]), merge=hmerge, order=ordr) }
      else if (moddendoinertie==TRUE) {
          resultscah=list(labels=dimnames(X)[[2]],inertie=inertie, height=delta, merge=hmerge,order=ordr,delta=delta )}
      class(resultscah)="hclust"
      mydendC=as.dendrogram(resultscah)
      #reorder
      resultscah$order <- sapply(1:p,function(j){which(dimnames(X)[[2]]==labels(mydendC)[j])})
      clv3Wclt$hclust <- resultscah
      clv3Wclt$param <- list(gmax=gmax,X=X)


  if (graph) {
        # graph of Delta
    delta <- results[,4]
    delta[(length(delta)-gmax+3):length(delta)] <- diff(results[(length(delta)-gmax+2):length(delta),7])

        dev.new() ##
    if (moddendoinertie==FALSE){ ##
      graphics::barplot(cumsum(delta)[(length(delta)-gmax+2):length(delta)],col=4,xlab="Nb clusters after aggregation", ylab="cum.delta", main="Evolution of the aggregation criterium",axisnames=TRUE,names.arg=paste(gmax:2,"->",(gmax-1):1),cex.names=0.5) }
        else if (moddendoinertie==TRUE) {
        graphics::barplot(delta[(length(delta)-gmax+2):length(delta)],col=4,xlab="Nb clusters after aggregation", ylab="delta", main="Evolution of the aggregation criterium",axisnames=TRUE,names.arg=paste(gmax:2,"->",(gmax-1):1),cex.names=0.5) }
        # Dendrogram #
        dev.new()
        plot(mydendC, type ="rectangle", xlab="", main="CLV3W Dendrogram")
  }
  clv3Wclt$call   <- match.call() ## changer le match.call
  glob.call <- as.list(clv3Wclt$call)
  if (!"mode.scale" %in% names(glob.call))
    glob.call$mode.scale = mode.scale
  if (!"NN" %in% names(glob.call))
    glob.call$NN = NN
  if (!"moddendoinertie" %in% names(glob.call))
    glob.call$moddendoinertie = moddendoinertie
  if (!"mode.gmax" %in% names(glob.call))
    glob.call$gmax = gmax
  clv3Wclt$call <- as.call(glob.call)

  class(clv3Wclt) <- c("clv3w","clv3wHCA")
  return(clv3Wclt)
}
