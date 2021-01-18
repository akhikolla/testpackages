Bagg_Surv <-
function(xdata, Y.names, P.names, T.names, method = 'R2', args.rpart, args.parallel = list(numWorkers = 1), Bag = 100) 
{
  #set.seed(seed) ; 
  methodSurvR2 = list(eval = .surveval, split = .survsplitR2, init = .survinit)
  methodSurvLR = list(eval = .surveval, split = .survsplitLR, init = .survinit)
  if(method == 'R2') methodsurv = methodSurvR2 ;
  if(method == 'LR') methodsurv = methodSurvLR ;
  
  xdata <- xdata[order(xdata[, Y.names[1]]), ] ;
  time1 <- Sys.time() ; n <- nrow(xdata) ; Ind_data <- 1:n ;
  row.names(xdata) <- Ind_data
  IND_SAMP <- lapply(1:Bag, function(u) sample.int(n, size = n, replace = TRUE)) ;
  IND_OOB <- lapply(IND_SAMP, function(v) setdiff(Ind_data, v)) ;
  
  wrapper <- function(xdata_bag) {
##    version a consolider pour tenir compte des variables agissant sur le plateau
#     beta <- coxph(Surv(xdata_bag[, Y.names[1]], xdata_bag[, Y.names[2]]) ~ xdata_bag[, P.names])$coefficients;
#     coveffect <- exp(colSums(beta * t(xdata_bag[, P.names]))) ;
#     YSTRATE <- coveffect ;
    
    # YSTRATE <- xdata_bag[, P.names] + 1 ;
    YSTRATE <- exp(xdata_bag[, P.names]) ;
    tree_surv <- rpart(as.formula(paste(paste('Surv(', paste(Y.names, collapse = ','),')', sep = ''), '~', paste(T.names, collapse = '+'), sep = '')), 
                       data = xdata_bag, method = methodsurv,  control = args.rpart, weights = YSTRATE, y = FALSE ) ;
#     tree_surv$functions <- NULL
    tree_surv$terms  <- NULL
    return(tree_surv)
  } 
  
  numWorkers = args.parallel$numWorkers ;
  cat("\n ncores = ", numWorkers, " for bagging trees !\n")
  List_xdatas <- lapply(IND_SAMP, function(w) return(xdata[w, ])) ;
  MaxTreeList <- mclapply(List_xdatas, wrapper, mc.cores = getOption("mc.cores", 
                                                                     numWorkers), mc.preschedule = TRUE, mc.silent = TRUE)
  
  tree_ind <- sapply(MaxTreeList, function(treee) inherits(treee, "rpart")) ;
  MaxTreeList <- MaxTreeList[tree_ind] ;
  IIND_SAMP <- IND_SAMP[tree_ind] ;
  IND_OOB <- IND_OOB[tree_ind] ;
  N_leafs <- sapply(MaxTreeList, function(xxtree) {
    return(sum(xxtree$frame$var == "<leaf>"))
  });
  MaxTreeList[which(N_leafs == 1)] <- NULL ;
  IIND_SAMP[which(N_leafs == 1)] <- NULL ;
  IND_OOB[which(N_leafs == 1)] <- NULL ;
  
  depthvar <- function(treee){
    varrs <- as.character(treee$frame$var) ;
    nodes <- as.numeric(row.names(treee$frame)) ;
    varrs <- varrs[order(nodes)] ;
    nodes <- sort(nodes) ;
    depth <- rep(0, max(nodes)) ;
    for(k in nodes[-1]){
      depth[k] <- depth[k/2] +1 ;
    }
    depth <- c(0, depth[depth!=0]) ;
    leafs <- which(varrs == '<leaf>') ;
    varrs <- varrs[-leafs] ;
    depth <- depth[-leafs] ;
    names(depth) <- varrs ;
    depth2 <- sort(tapply(depth, names(depth), min)) ;
    return(list(depth = depth, depth2 = depth2)) ;  
  } 
  
  chenIMP <- function(treee){
    improvee <- treee$splits[, 'improve']*treee$frame[treee$frame[, 'var']!="<leaf>", ][, 'yval'] ;
    DEPTH1 <- depthvar(treee)$depth ;
    if(is.null(names(improvee))) names(improvee) <- names(DEPTH1) ;  
    DEPTH1 <- DEPTH1[names(improvee)] ;
    Chen_IMP <- 2^(-DEPTH1)*improvee ;
    DEV_IMP <- improvee ;
    return(list(Chen_IMP = Chen_IMP, DEV_IMP = DEV_IMP))
  }
  
  DIS <- lapply(MaxTreeList, function(ww) chenIMP(ww)$DEV_IMP) ;
  DIS <- sort(tapply(unlist(DIS), names(unlist(DIS)), sum), decreasing = T) ;
  DIS <- DIS/sum(DIS)*100 ;
  
  DDIS <- lapply(MaxTreeList, function(ww) chenIMP(ww)$Chen_IMP) ;
  DDIS <- sort(tapply(unlist(DDIS), names(unlist(DDIS)), sum), decreasing = T) ;
  DDIS <- DDIS/sum(DDIS)*100 ;
  
  depthBAG <- lapply(MaxTreeList, function(treee) depthvar(treee)$depth2) ;
  depth <- unlist(depthBAG) ;
  DEPTH <- sort(tapply(depth, names(depth), mean)) ;
  
  NEWBAG <- length(MaxTreeList) ;
  
#   LIST_VAR <- lapply(1 : NEWBAG, function(zzz){
#     return(unique(as.character(MaxTreeList[[zzz]]$frame$var[MaxTreeList[[zzz]]$frame$var!= '<leaf>'])))
#   }) ;
#   LIST_VARR <- names(DIS) ;
#   OCCUR <- t(sapply(LIST_VAR, function(uu){
#     return(LIST_VARR %in% uu)
#   }))
#   colnames(OCCUR) <- LIST_VARR ; 
  
  time2 <- Sys.time() ;
  Timediff <- difftime(time2, time1) ;
  
  return(list(MaxTreeList = MaxTreeList, IIS = DIS, DIIS = DDIS, DEPTH = DEPTH, IND_OOB = IND_OOB,
              IIND_SAMP = IIND_SAMP, IND_SAMP = IND_SAMP, Bag = NEWBAG, indrpart = tree_ind, Timediff = Timediff)) ;
}
