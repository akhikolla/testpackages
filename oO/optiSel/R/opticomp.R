
"opticomp"<-function(f, phen, obj.fun="NGD", lb=NULL, ub=NULL, ...){
  phenAsDataTable <- "data.table" %in% class(phen)
  phen <- as.data.frame(phen)
  if(phenAsDataTable){setDF(phen)}
  if(!("data.frame" %in% class(phen))){
    stop("Argument 'phen' is not a data frame.\n")
  }
  if(!("Indiv" %in% colnames(phen))){
    stop("Column 'Indiv' with IDs of individuals is missing in data frame 'phen'.\n")
  }
  phen$Indiv <- as.character(phen$Indiv)
  if(any(is.na(phen$Indiv))){
    stop("Some Individual IDs are NA in data frame phen.\n")
  }
  if(any(duplicated(phen$Indiv))){
    stop("Some Individuals appear twice in data frame phen.\n")
  }
  rownames(phen)<-phen$Indiv
  if(!("Breed" %in% colnames(phen))){
    stop("Column 'Breed' is missing in data frame 'phen'.\n")
  }
  if(!is.matrix(f)){stop("f is not a matrix.\n")}
  if(!all(phen$Indiv %in% colnames(f))){
    stop("Some individuals from phen are missing in f.\n")
  }
  
  f <- f[phen$Indiv, phen$Indiv]
  X <- aggregate(f, list(phen$Breed), mean)
  rownames(X)<-X[,"Group.1"]
  X <- t(as.matrix(X[,-1]))
  X <- aggregate(X, list(phen$Breed), mean)
  rownames(X)<-X[,"Group.1"]
  f <- as.matrix(X[,-1])
  
  l <- rep(0,ncol(f))
  u <- rep(1,ncol(f))
  names(l)<-colnames(f)
  names(u)<-colnames(f)
  if(!is.null(lb)){
    if(!all(names(lb)%in%phen$Breed)){stop("Some breeds in lb are missing in phen.\n")}
    if(any(lb<0)){stop("Some components of lb are negative.\n")}
    if(sum(lb)>1){stop("sum(lb)>1.\n")}
    l[names(lb)]<-lb
    }
  if(!is.null(ub)){
    if(!all(names(ub)%in%phen$Breed)){stop("Some breeds in ub are missing in phen.\n")}
    idb <- intersect(names(lb), names(ub))
    if(any(lb[idb]>ub[idb])){stop("Some components of ub are smaller than lb\n")}
    u[names(ub)]<-ub
    }
  H <- matrix(c(f),ncol=ncol(f),nrow=nrow(f))
  gc()
  # maximize neutral gene diversity
  if(obj.fun=="NGD"){
    Res   <- solve.QP(Dmat=2*H,dvec=rep(0,nrow(H)),Amat=cbind(diag(nrow(H)),-diag(nrow(H)),1),bvec=c(l,-u,1), ...)
    bc    <- setNames(Res$solution, colnames(f))
    value <- (1-Res$value)
  }
  #maximize neutral trait diversity
  if(obj.fun=="NTD"){
    F    <- diag(f)
    eins <- rep(1,nrow(f))
    Dmat <- 2*(eins%*%t(eins)-(F%*%t(eins)-2*f+eins%*%t(F)))
    dvec <- -F
    Res  <- solve.QP(Dmat=Dmat,dvec=dvec,Amat=cbind(diag(nrow(H)),-diag(nrow(H)),1),bvec=c(l,-u,1), ...)
    bc   <- setNames(Res$solution, colnames(f))
    value<- t(bc)%*%(eins-F)+t(bc)%*%(F%*%t(eins)-2*f+eins%*%t(F))%*%bc
  }
   
  #Compute genetic distance between phen$Breeds
  A <- matrix(diag(f),ncol=ncol(f),nrow=nrow(f),byrow=TRUE)
  Dist <- sqrt((A+t(A))/2-f)
  
  list(bc=bc, value=value, f=f,  Dist=Dist)
}

