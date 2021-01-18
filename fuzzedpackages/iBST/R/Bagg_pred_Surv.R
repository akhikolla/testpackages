Bagg_pred_Surv <-
function(xdata, Y.names, P.names, resBag, args.parallel = list(numWorkers = 1), new_data = data.frame(), OOB = FALSE)
{ 
  time1 <- Sys.time() ;
  numWorkers <- args.parallel$numWorkers ;
  xdata <- xdata[order(xdata[, Y.names[1]]), ] ;
  n <- dim(xdata)[1] ;
  row.names( xdata ) <- 1:n ;
  time_vec <- unique(xdata[, Y.names[1]]) ;
  # statuss <- xdata[, Y.names[2]];
  List_xTrees_xDatas <- lapply(1 : resBag$Bag, function(j) {
    return(list(resBag$MaxTreeList[[j]], xdata[resBag$IIND_SAMP[[j]], ]))
  }) ;
  
  ## Prediction du risque cumule dans un arbre via Nelson Aalen et Breslow
  Mat <- matrix(nrow = length(time_vec), ncol = 2) ;
  Mat[, 1] <- time_vec ; Mat[, 2] <- 0 ;
  
  ## cumulative risk estimated via Aalen and Breslow for a given tree
  wrapper1 <- function(tree_data){
    tempdata <- tree_data[[2]] ;
    temptree <- tree_data[[1]] ;
    tempind <- tree2indicators(temptree);
    
    ## the plateau effect 
    beta <- coxph(Surv(tempdata[, Y.names[1]], tempdata[, Y.names[2]]) ~ strata(temptree$where) + tempdata[, P.names])$coefficients;
    
    ## list of items for each leaf node
    list_leaf <- lapply(sort(unique(temptree$where)), function(zz) {
      floor(as.numeric(names(temptree$where[temptree$where == zz])))
    })
    names( list_leaf ) <- unlist( tempind ) ;
    haz_list <- lapply(list_leaf, function(izi){
      leafdat <- xdata[izi, ] ; 
      leafdat <- leafdat[order(leafdat[, Y.names[1]]), ] ;
      
      ## Pour Nelson AAlen
      IRisk <- sapply(leafdat[, Y.names[1]], function(uu) as.integer(leafdat[, Y.names[1]] >= uu)) ;
      NRISK <- colSums(IRisk) ;
      MARGRISK1 <- leafdat[, Y.names[2]]/NRISK ; 
      MARGRISK1 <- tapply(MARGRISK1, leafdat[, Y.names[1]], sum) ;
      Mat_NAa <- Mat ; Mat_NAa[time_vec %in% (unique(leafdat[, Y.names[1]])), 2] <- MARGRISK1 ;
      Mat_NAa[, 2] <- cumsum(Mat_NAa[, 2]) ;
      
      ## Pour Breslow
      coveffect <- exp(colSums(beta * t(leafdat[, P.names]))) ;
      RISKCOV <- t(IRisk) %*% coveffect ;
      MARGRISK2 <- (leafdat[, Y.names[2]]/RISKCOV)*coveffect ;
      MARGRISK2 <- tapply(MARGRISK2, leafdat[, Y.names[1]], sum) ;
      Mat_BRe <- Mat ; Mat_BRe[time_vec %in% unique(leafdat[, Y.names[1]]), 2] <- MARGRISK2 ;
      Mat_BRe[, 2] <- cumsum(Mat_BRe[, 2]) ;
      return(list(Mat_NAa, Mat_BRe)) ;
    })
    
    names( haz_list ) = names( list_leaf ) ;
    return( haz_list )
  }
  
  ## list with a tree prediction at first level, the leaf prediction at second level
  ## and the NAa or Bre estimator at third level
  cat("\n ncores = ", numWorkers, " for prediction !\n")
  
  List_Haz_Bagg <- mclapply(List_xTrees_xDatas, wrapper1, mc.cores = getOption("mc.cores", 
                                                                               numWorkers), mc.preschedule = TRUE, mc.silent = TRUE)
  ## NAa estimate for each tree and each node of the tree
  List_Haz_BaggNAa <-  lapply(List_Haz_Bagg, function(zzz){ lapply(zzz, function(uuu) uuu[[1]])})
  
  mathazNa = lapply(List_Haz_BaggNAa, function(wz){ sapply(wz, function(uu) uu[, 2]) })
  
  tabhazNAa = lapply(mathazNa, function(cyp) {matt = cyp;
                                              colnames(matt) = paste('leaf', 1:dim(cyp)[2], sep = '');
                                              return(matt)})
  
  List_Haz_BaggBRe <-  lapply(List_Haz_Bagg, function(zzz){ lapply(zzz, function(uuu) uuu[[2]])})
  
  mathazBr = lapply(List_Haz_BaggBRe, function(wz){ sapply(wz, function(uu) uu[, 2]) })
  
  tabhazBRe = lapply(mathazBr, function(cyp) {matt = cyp;
                                              colnames(matt) = paste('leaf', 1:dim(cyp)[2], sep = '') ;
                                              return(matt)})
  
  rm(mathazNa, mathazBr, List_Haz_Bagg, List_Haz_BaggNAa,  List_Haz_BaggBRe, List_xTrees_xDatas) ;
  
  ## predicted cumulative risk for each individual of the learning sample per column, 
  ## list for all the Aalen  predictors
  INDICATORS <- lapply(resBag$MaxTreeList, function(ztz){
    return(tree2indicators(ztz)) ;
  }) ;
  

  predNAa <- lapply(1 : resBag$Bag, function(uuu){
    tempind <- INDICATORS[[uuu]] ;
    vecc <- sapply(tempind, function(vv) with(xdata, eval(parse(text = vv)))) ;
    vectt = apply(vecc, 1, function(u) which(u == TRUE)) ;
    hazz <-  (tabhazNAa[[uuu]])[, vectt] ;
    return(hazz) ;
  })

  ## predicted cumulative risk for each individual of the learning sample per column,
  ## list for all the Breslow predictors

  predBR <- lapply(1 : resBag$Bag, function(uuu){
    tempind <- INDICATORS[[uuu]] ;
    vecc <- sapply(tempind, function(vv) with(xdata, eval(parse(text = vv)))) ;
    vectt = apply(vecc, 1, function(u) which(u == TRUE)) ;
    hazz <-  (tabhazBRe[[uuu]])[, vectt] ;
    return(hazz)
  })


  ## individual prediction of cumulative risk via the ensemble method
  ## using the learning sample
  ## *********** Fantastic  *********** ##

  PREDNA <- matrix(nrow = length(time_vec), ncol = dim(xdata)[1]) ;
  PREDBRE <- matrix(nrow = length(time_vec), ncol = dim(xdata)[1]) ;

  for(i in 1:dim(xdata)[1]){
    predna1 <- sapply(predNAa, function(vvv) vvv[, i])
    PREDNA[, i] <- rowMeans(predna1)
    predbr1 <- sapply(predBR, function(vvv) vvv[, i])
    PREDBRE[, i] <- rowMeans(predbr1)
  }

  rm(predNAa, predBR, predna1, predbr1 ) ;
  
  ## OOB prediction with OOB estimation for each permuted variable
  
  if (OOB) {
    UNIQUE_IND_OOB <- sort(unique(unlist(resBag$IND_OOB)))
    
    if(length(UNIQUE_IND_OOB) != n) stop('Please add the number of trees in the bagging sequence!!')
    
    LIST_VAR <- lapply(1:resBag$Bag, function(zzz){
      return(as.character(resBag$MaxTreeList[[zzz]]$frame$var[resBag$MaxTreeList[[zzz]]$frame$var!= '<leaf>']))
    }) ;
    
    LIST_VAR = unique(unlist(LIST_VAR)) ;
    LIST_VAR_BRE <- paste(LIST_VAR, 'BRE', sep = '')
    
#     LIST_DATA_VAR <- lapply(LIST_VAR, function(ww){
#       temp = xdata ;
#       temp[, ww] <- xdata[, ww][sample(1:n)]
#       return(temp)
#     }) ;
    
    wrapper2 <- function(uvw){
      indexoob <- sapply(resBag$IND_OOB, function(vvv) is.element(uvw,vvv)) ;
      treelistoob <- resBag$MaxTreeList[indexoob] ;
      indicatorsoob <- INDICATORS[indexoob]
      tabhazNAaoob <- tabhazNAa[indexoob]
      tabhazBReoob <- tabhazBRe[indexoob]
      
      newdataoob <- xdata[uvw, ] ;
      
      predNAaoob <- sapply(1 : length(treelistoob), function(uuu){
        tempind <- indicatorsoob[[uuu]] ;
        vecc <- sapply(tempind, function(vv) with(newdataoob, eval(parse(text = vv)))) ;
        vectt = which(vecc == TRUE) ;
        hazz <-  (tabhazNAaoob[[uuu]])[, vectt] ;
        return(hazz) ;
      })
      
      predBRoob <- sapply(1 : length(treelistoob), function(uuu){
        tempind <- indicatorsoob[[uuu]] ;
        vecc <- sapply(tempind, function(vv) with(newdataoob, eval(parse(text = vv)))) ;
        vectt = which(vecc == TRUE) ;
        hazz <-  (tabhazBReoob[[uuu]])[, vectt] ;
        return(hazz) ;
      })
      ## prediction OOB pour un idividu via la methode ensembliste
      predNAoob <- rowMeans(predNAaoob) ;
      predBREoob <- rowMeans(predBRoob) ; 
      
#       PREDVARNAOOB <- sapply(1:length(LIST_VAR), function(var){
#         predvaoob <- sapply(1 : length(treelistoob), function(uuu){
#           tempind <- indicatorsoob[[uuu]] ;
#           vecc <- sapply(tempind, function(vv) with(LIST_DATA_VAR[[var]][uvw, ], eval(parse(text = vv)))) ;
#           vectt = which(vecc == TRUE) ;
#           hazz <-  (tabhazNAaoob[[uuu]])[, vectt] ;
#           return(hazz) ;
#         }) ;
#         PREDVAROOB <- rowMeans(predvaoob) ;
#         return(PREDVAROOB) ;
#       }) ;
      
#       PREDVARBREOOB <- sapply(1:length(LIST_VAR), function(var){
#         predvaoob <- sapply(1 : length(treelistoob), function(uuu){
#           tempind <- indicatorsoob[[uuu]] ;
#           vecc <- sapply(tempind, function(vv) with(LIST_DATA_VAR[[var]][uvw, ], eval(parse(text = vv)))) ;
#           vectt = which(vecc == TRUE) ;
#           hazz <-  (tabhazBReoob[[uuu]])[, vectt] ;
#           return(hazz) ;
#         }) ;
#         PREDVAROOB <- rowMeans(predvaoob)
#         return(PREDVAROOB)
#       })
      
#       PREDMATHAZZ <- cbind(predNAoob, PREDVARNAOOB, predBREoob, PREDVARBREOOB) ;
#       colnames(PREDMATHAZZ) <- c('oobNA', LIST_VAR, 'oobBR', LIST_VAR_BRE) ;
      PREDMATHAZZ <- cbind(predNAoob, predBREoob)
      colnames(PREDMATHAZZ) <- c('oobNA', 'oobBR')
      
      return(PREDMATHAZZ) ;
        
    } ;
    
    cat("\n ncores = ", numWorkers, " for oob (VIMP + OOB) !\n")
    
    CUMHAZOOB <- mclapply(UNIQUE_IND_OOB, wrapper2, mc.cores = getOption("mc.cores", numWorkers), 
                          mc.preschedule = TRUE, mc.silent = TRUE) ;
    
    HAZOOBNA <- sapply(CUMHAZOOB, function(predmat){return(predmat[, 'oobNA'])}) ;
    HAZOOBBRE <- sapply(CUMHAZOOB, function(predmat){return(predmat[, 'oobBR'])}) ;
    
#     HAZVARNA <- list()
#     HAZVARBRE <- list()
#     for(j in 1:length(LIST_VAR)){
#       HAZVARNA[[j]] <- sapply(CUMHAZOOB, function(predmat){return(predmat[, LIST_VAR[j]])}) ;
#       HAZVARBRE[[j]] <- sapply(CUMHAZOOB, function(predmat){return(predmat[, LIST_VAR_BRE[j]])}) ;
#     }
    
    ## Compute the OOB predictor by predictor
    wrapper3 <- function(cyp)
      {
      newdatatemp <- xdata[resBag$IND_OOB[[cyp]], ] ;
      newdatatemp <- newdatatemp[order(newdatatemp[, Y.names[1]]), ] ;
      nn <- dim(newdatatemp)[1] ;
      time_veccc <- unique(newdatatemp[, Y.names[1]]);
      tempind <- INDICATORS[[cyp]] ;
      vecc <- sapply(tempind, function(vv) with(newdatatemp, eval(parse(text = vv)))) ;
      vectt = apply(vecc, 1, function(u) which(u == TRUE)) ;
      hazznaa <-  (tabhazNAa[[cyp]])[, vectt] ;
      hazzbre <-  (tabhazBRe[[cyp]])[, vectt] ;
      survna <- exp(-hazznaa) ; survbre <- exp(-hazzbre) ; 
      
      ## position of time_veccc within the initial time-vec
      IND <- sapply(time_veccc, function(vv) sum(time_vec <= vv))
      survna <- survna[IND, ]
      survbre <- survbre[IND, ]
      
      Riskk <- t(sapply(time_veccc, function(uu)  as.integer(newdatatemp[, Y.names[1]] >= uu))) ;
     
      ## helpfull to check the dimensional equation 
      rsk <- sapply(time_veccc, function(uu)  as.integer(newdatatemp[, Y.names[1]] <= uu)) ;
      RRISKK <- t(rsk*newdatatemp[, Y.names[2]]) ;
      
      NRISKK <- colSums(t(Riskk)) ;
      cens <- tapply( (1- newdatatemp[, Y.names[2]]),  newdatatemp[, Y.names[1]], sum)
      # ndeath <- tapply( ( newdatatemp[, Y.names[2]]),  newdatatemp[, Y.names[1]], sum)
      MARGRISK1 <- cens/NRISKK ;
      # MARGRISKM <- ndeath/NRISKK ;
      
      KMSURVC <- exp(-cumsum(MARGRISK1));
      # SURVKMnew <- exp(-cumsum(MARGRISKM));
      
      ## faire que KMSURV ait la taille des donnees test de base
      vectrep <- colSums(sapply(time_veccc, function(uu) uu == newdatatemp[, Y.names[1]]))
      KMSURVCI <- rep(KMSURVC, vectrep)
      ## KMSURVCnew <- rep(KMSURVnew, vectrep)
      ## Score de Brier en considerant la survie censuree estimee via KM
      
      matbstnaoob <- t(t(survna**2 * RRISKK) * (1/KMSURVCI)) + (1 - survna)**2 * t(1 - rsk) * (1/KMSURVC) ;
      
      bstnaoob <- rowMeans(matbstnaoob);
      
      matbstbreoob <- t(t(survbre**2 * RRISKK) * (1/KMSURVCI)) + (1 - survbre)**2 * t(1 - rsk) * (1/KMSURVC) ;
      
      bstbreoob  <- rowMeans(matbstbreoob);
      
      ibsnaoob = sum(diff(time_veccc) * (bstnaoob[-1] + bstnaoob[-length(bstnaoob)]) / 2)/max(time_veccc);
      
      ibsbreoob = sum(diff(time_veccc) * (bstbreoob[-1] + bstbreoob[-length(bstbreoob)]) / 2)/max(time_veccc);
      
      ## compute the IBS on OOB data after permuting each variable
      ibsoobvarpbp <- sapply(1:length(LIST_VAR), function(var){
        tempdat <- newdatatemp ;
        tempdat[, LIST_VAR[var]] <- newdatatemp[, LIST_VAR[var]][sample(1:nn)] ;
        
        vecc <- sapply(tempind, function(vv) with(tempdat, eval(parse(text = vv)))) ;
        vectt = apply(vecc, 1, function(u) which(u == TRUE)) ;
        hazznaa <-  (tabhazNAa[[cyp]])[, vectt] ;
        hazzbre <-  (tabhazBRe[[cyp]])[, vectt] ;
        survna <- exp(-hazznaa) ; survbre <- exp(-hazzbre) ; 
        survna <- survna[IND, ] ; survbre <- survbre[IND, ] ;
        
        matbstnaoob <- t(t(survna**2 * RRISKK) * (1/KMSURVCI)) + (1 - survna)**2 * t(1 - rsk) * (1/KMSURVC) ;
        bstnaoob <- rowMeans(matbstnaoob);
        
        matbstbreoob <- t(t(survbre**2 * RRISKK) * (1/KMSURVCI)) + (1 - survbre)**2 * t(1 - rsk) * (1/KMSURVC) ;
        bstbreoob  <- rowMeans(matbstbreoob);
        
        ibsnaoobv = sum(diff(time_veccc) * (bstnaoob[-1] + bstnaoob[-length(bstnaoob)]) / 2)/max(time_veccc);
        
        ibsbreoobv = sum(diff(time_veccc) * (bstbreoob[-1] + bstbreoob[-length(bstbreoob)]) / 2)/max(time_veccc);
        
        return(c(ibsnaoobv, ibsbreoobv))
      })
      
      ibstree <- cbind(c(ibsnaoob, ibsbreoob), ibsoobvarpbp) ;
      
      return(ibstree) ;
      }
    
 IBSOOBPBP <- mcmapply(wrapper3, 1:length(resBag$IND_OOB), mc.cores = getOption("mc.cores", numWorkers), mc.silent = TRUE) ;
    
    ## the IBS corresponding to Breslow are available in even indexes
    indbre <- which((1:nrow(IBSOOBPBP) %% 2) == 0) ;
    IBSOOBPBPBRE <- IBSOOBPBP[ indbre, ] ;
    oobibspbpbre <- IBSOOBPBPBRE[1, ] ;
    ## mean oob error when a new tree is added
    oobibspbpbre <- cumsum(oobibspbpbre)/1:length(oobibspbpbre) ;
    
    IBSOOBPBPBRE <- apply(IBSOOBPBPBRE, 1, mean) ;
    
    IBSOOBPBPNA <- IBSOOBPBP[setdiff(1:nrow(IBSOOBPBP), indbre), ] ;
    oobibspbpna <- IBSOOBPBPNA[1, ] ;
    ## mean oob error when a new tree is added
    oobibspbpna <- cumsum(oobibspbpna)/1:length(oobibspbpna) ;
    IBSOOBPBPNA <- apply(IBSOOBPBPNA, 1, mean) ;
    
    vimpoobpbpna <- IBSOOBPBPNA[-1] - IBSOOBPBPNA[1] ;
    names(vimpoobpbpna) <- LIST_VAR ;
    vimpoobpbpna <- sort(vimpoobpbpna, decreasing = T) ;
    
    vimpoobpbpbre <- IBSOOBPBPBRE[-1] - IBSOOBPBPBRE[1] ;
    names(vimpoobpbpbre) <- LIST_VAR ;
    vimpoobpbpbre <- sort(vimpoobpbpbre, decreasing = TRUE) ;
  }
  
  
  
  ## prediction on a new dataset ( Test sample)
  ## **************  Fantastic  ******************  ##
  
  if (nrow(new_data) != 0) {
    new_data <- new_data[order(new_data[, Y.names[1]]), ]
    time_vecnew <- unique(new_data[, Y.names[1]]);
    
    ## list of predictions for each tree
    
    predNAanew <- lapply(1 : resBag$Bag, function(uuu){
      tempind <- INDICATORS[[uuu]] ;
      vecc <- sapply(tempind, function(vv) with(new_data, eval(parse(text = vv)))) ;
      vectt = apply(vecc, 1, function(u) which(u == TRUE)) ;
      hazz <-  (tabhazNAa[[uuu]])[, vectt] ;
      return(hazz) ;
    })
    
    predBRnew <- lapply(1 : resBag$Bag, function(uuu){
      tempind <- INDICATORS[[uuu]] ;
      vecc <- sapply(tempind, function(vv) with(new_data, eval(parse(text = vv)))) ;
      vectt = apply(vecc, 1, function(u) which(u == TRUE)) ;
      hazz <-  (tabhazBRe[[uuu]])[, vectt] ;
      return(hazz)
    }) 
    
    PREDNAnew <- matrix(nrow = length(time_vec), ncol = dim(new_data)[1]) ;
    PREDBREnew <- matrix(nrow = length(time_vec), ncol = dim(new_data)[1]) ;
    for(j in 1:dim(new_data)[1]){
      predna1 <- sapply(predNAanew, function(vvv) vvv[, j])
      PREDNAnew[, j] <- rowMeans(predna1)
      predbr1 <- sapply(predBRnew, function(vvv) vvv[, j])
      PREDBREnew[, j] <- rowMeans(predbr1)
    }
    
    rm(predNAanew, predBRnew, predna1, predbr1) ;
    
    ## prediction of the cumulative risk for the new observed times
    
    Time_index <- sapply(new_data[, Y.names[1]] , function(wz) sum(time_vec <= wz))
    Time_index[Time_index == 0] <- 1
    SURV_NAnew <- exp(-diag(PREDNAnew[Time_index, ]));
    SURV_BREnew <- exp(-diag(PREDBREnew[Time_index, ]));
    
    ## Integrated Brier score for new data
    
    SURVNAnew <- exp(-PREDNAnew) ; SURVBREnew <- exp(-PREDBREnew) ; 
    IND <- sapply(time_vecnew, function(vv) sum(time_vec <= vv))
    IND[IND == 0] <- 1
    SURVNAnew <- SURVNAnew[IND, ]
    SURVBREnew <- SURVBREnew[IND, ]
    
    Risknew <- t(sapply(time_vecnew, function(uu)  as.integer(new_data[, Y.names[1]] >= uu))) ;
    rsknew <- sapply(time_vecnew, function(uu)  as.integer(new_data[, Y.names[1]] <= uu)) ;
    RRISKnew <- t(rsknew*new_data[, Y.names[2]]) ;
    
    NRISKnew <- colSums(t(Risknew)) ;
    censnew <- tapply( (1- new_data[, Y.names[2]]),  new_data[, Y.names[1]], sum)
    ndeathnew <- tapply( ( new_data[, Y.names[2]]),  new_data[, Y.names[1]], sum)
    MARGRISK1new <- censnew/NRISKnew ;
    MARGRISKMnew <- ndeathnew/NRISKnew ;
    
    KMSURVCnew <- exp(-cumsum(MARGRISK1new));
    
    SURVKMnew <- exp(-cumsum(MARGRISKMnew));
    
    ## faire que KMSURV ait la taille des donnees test de base
    vectrep <- colSums(sapply(time_vecnew, function(uu) uu == new_data[, Y.names[1]]))
    KMSURVCInew <- rep(KMSURVCnew, vectrep)
    ## KMSURVCnew <- rep(KMSURVnew, vectrep)
    ## Score de Brier en considerant la survie censuree estimee via KM
    
    MATBSTNAKMnew <- t(t(SURVNAnew**2 * RRISKnew) * (1/KMSURVCInew)) + (1 - SURVNAnew)**2 * t(1 - rsknew) * (1/KMSURVCnew) ;
    BSTNAKMnew <- rowMeans(MATBSTNAKMnew);
    
    MATBSTBREKMnew <- t(t(SURVBREnew**2 * RRISKnew) * (1/KMSURVCInew)) + (1 - SURVBREnew)**2 * t(1 - rsknew) * (1/KMSURVCnew);
    BSTBREKMnew <- rowMeans(MATBSTBREKMnew);
    
    MATSURVKMnew <- matrix(SURVKMnew, nrow = length(time_vecnew), ncol = length(new_data[, Y.names[1]]), byrow = FALSE) ;
    MATBSTKMnew <- t(t(MATSURVKMnew**2 * RRISKnew) * (1/KMSURVCInew)) + (1 - MATSURVKMnew)**2 * t(1 - rsknew) * (1/KMSURVCnew) ;
    BSTKMnew <- rowMeans(MATBSTKMnew);
    
    
    
    IBSNAKMnew = sum(diff(time_vecnew) * (BSTNAKMnew[-1] + BSTNAKMnew[-length(BSTNAKMnew)]) / 2)/max(time_vecnew);
    
    IBSBRKMnew = sum(diff(time_vecnew) * (BSTBREKMnew[-1] + BSTBREKMnew[-length(BSTBREKMnew)]) / 2)/max(time_vecnew);
    
    IBSKMnew = sum(diff(time_vecnew) * (BSTKMnew[-1] + BSTKMnew[-length(BSTKMnew)]) / 2)/max(time_vecnew);
    
    rm(BSTBREKMnew, BSTNAKMnew, MATBSTNAKMnew, MATBSTBREKMnew, KMSURVCInew, KMSURVCnew,
       Risknew, RRISKnew, MATSURVKMnew)
  }
  
  
  
  ## *** Procedure pour le calcul de la survie pour la  censure (necessaire pour Brier!!!) *** ##
  ## *** il me semble assez naturel d utiliser un estimateur type KM sur l echantillon d apprentissage!? *** ##
  ## *** Trois approches sont evaluees :  KM sur echantillon d apprentissage, Nelson Aalen via la *** ## 
  ## *** methode ensembliste en inversant les roles (observations - censure), et Breslow via la methode ensembliste ... *** ##
  #   
#   SURVNA <- exp(-PREDNA) ; SURVBRE <- exp(-PREDBRE) ; 
  # SURVNAC <- exp(-PREDNAC) ; SURVBREC <- exp(-PREDBREC) ;
  Risk <- t(sapply(time_vec, function(uu)  as.integer(xdata[, Y.names[1]] >= uu))) ;
  ## helpfull to check the dimensional equation 
  RSK <- sapply(time_vec, function(uu)  as.integer(xdata[, Y.names[1]] <= uu)) ;
  RRISK <- t(RSK*xdata[, Y.names[2]]) ;
  
  #   
  ## estimate the survival function of censored observations with the KM estimated
  NRISK <- colSums(t(Risk)) ;
  cens <- tapply(1- xdata[, Y.names[2]],  xdata[, Y.names[1]], sum)
  death <- tapply(xdata[, Y.names[2]],  xdata[, Y.names[1]], sum)
  MARGRISK1 <- cens/NRISK ; 
  MARGRISKKM <- death/NRISK ;
#   KMSURVC <- survfit(Surv(xdata[, Y.names[1]], 1 - xdata[, Y.names[2]]) ~ 1)$surv
  KMSURVC <- exp(-cumsum(MARGRISK1));
  SURVKM <- exp(-cumsum(MARGRISKKM));
  ## faire que KMSURVC ait la taille des donnees de base
  vectrep <- colSums(sapply(time_vec, function(uu) uu == xdata[, Y.names[1]]))
  KMSURVCI <- rep(KMSURVC, vectrep)
  ## Score de Brier en considerant la survie censuree estimee via KM
  
#   MATBSTNAKM <- SURVNA**2 * RRISK * matrix(1/KMSURVCI, nrow = length(time_vec), ncol = length(KMSURVCI), byrow = TRUE)
#   + (1 - SURVNA)**2 * (1 - Risk) * matrix(1/KMSURVC, nrow = length(time_vec), ncol = length(xdata[, Y.names[1]]), byrow = FALSE);
#   BSTNAKM <- rowMeans(MATBSTNAKM);
#   
#   MATBSTBREKM <- SURVBRE**2 * RRISK * matrix(1/KMSURVCI, nrow = length(time_vec), ncol = length(KMSURVCI), byrow = TRUE) 
#   + (1 - SURVBRE)**2 * (1 - Risk) * matrix(1/KMSURVC, nrow = length(time_vec), ncol = length(xdata[, Y.names[1]]), byrow = FALSE);
#   BSTBREKM <- rowMeans(MATBSTBREKM);
#   
  MATSURVKM <- matrix(SURVKM, nrow = length(time_vec), ncol = length(xdata[, Y.names[1]]), byrow = FALSE) ;
  MATBSTKM <- t(t(MATSURVKM**2 * RRISK) * (1/KMSURVCI)) + (1 - MATSURVKM)**2 * t(1 - RSK) * (1/KMSURVC) ;
  BSTKM <- rowMeans(MATBSTKM);
  
  if(OOB){
    SURVNAOOB <- exp(- HAZOOBNA ) ;
    MATBSTNAOOB <- t(t(SURVNAOOB**2 * RRISK) * (1/KMSURVCI)) + (1 - SURVNAOOB)**2 * t(1 - RSK) * (1/KMSURVC) ;
    BSTNAOOB <- rowMeans(MATBSTNAOOB);
    IBSNAOOB <- sum(diff(time_vec) * (BSTNAOOB[-1] + BSTNAOOB[-length(BSTNAOOB)]) / 2)/max(time_vec);
    
    SURVBREOOB <- exp(- HAZOOBBRE ) ;
    MATBSTBREOOB <- t(t(SURVBREOOB**2 * RRISK )* (1/KMSURVCI)) + (1 - SURVBREOOB)**2 * t(1 - RSK) * (1/KMSURVC);
    BSTBREOOB <- rowMeans(MATBSTBREOOB);
    IBSBREOOB <- sum(diff(time_vec) * (BSTBREOOB[-1] + BSTBREOOB[-length(BSTBREOOB)]) / 2)/max(time_vec);
    
#     IBSVAR <- sapply(1:length(LIST_VAR), function(uvu){
#       SURVVARNA <- exp(-HAZVARNA[[uvu]]) ;
#       MATBSTNAVAR <- t(t(SURVVARNA**2 * RRISK) * (1/KMSURVCI)) + (1 - SURVVARNA)**2 * t(1 - RSK) * (1/KMSURVC);
#       BSTNAVAR <- rowMeans(MATBSTNAVAR);
#       IBSNAVAR <- sum(diff(time_vec) * (BSTNAVAR[-1] + BSTNAVAR[-length(BSTNAVAR)]) / 2)/max(time_vec) ;
#       
#       SURVVARBRE <- exp(-HAZVARBRE[[uvu]]) ;
#       MATBSTBREVAR <- t(t(SURVVARBRE**2 * RRISK) * (1/KMSURVCI)) + (1 - SURVVARBRE)**2 * t(1 - RSK) * (1/KMSURVC);
#       BSTBREVAR <- rowMeans(MATBSTNAVAR);
#       IBSBREVAR <- sum(diff(time_vec) * (BSTBREVAR[-1] + BSTBREVAR[-length(BSTBREVAR)]) / 2)/max(time_vec) ;
#       return(c(IBSNAVAR, IBSBREVAR))
#     }) ;
    
#     VIMPIBSNA <- IBSVAR[1, ] - IBSNAOOB ;
#     names(VIMPIBSNA) <- LIST_VAR ;
#     VIMPIBSNA  <- sort(VIMPIBSNA, decreasing = TRUE )
#     
#     VIMPIBSBRE <- IBSVAR[2, ] - IBSBREOOB ;
#     names(VIMPIBSBRE) <- LIST_VAR ;  
#     VIMPIBSBRE <- sort(VIMPIBSBRE, decreasing = TRUE)
  }
  
#   IBSNAKM = sum(diff(time_vec) * (BSTNAKM[-1] + BSTNAKM[-length(BSTNAKM)]) / 2)/max(time_vec);
#   
#   IBSBRKM = sum(diff(time_vec) * (BSTBREKM[-1] + BSTBREKM[-length(BSTBREKM)]) / 2)/max(time_vec);
#   
  IBSKM = sum(diff(time_vec) * (BSTKM[-1] + BSTKM[-length(BSTKM)]) / 2)/max(time_vec);
  
  time2 <- Sys.time() ;
  Timediff <- difftime(time2, time1) ;
  

return(list(PREDNA = PREDNA, PREDBRE = PREDBRE, tabhazNAa = tabhazNAa, tabhazBRe = tabhazBRe, OOB = ( if(OOB) list( IBSKM = IBSKM, IBSNAOOB = IBSNAOOB, IBSBREOOB = IBSBREOOB,  vimpoobpbpna = vimpoobpbpna, vimpoobpbpbre = vimpoobpbpbre, 
                                 oobibspbpna = oobibspbpna, oobibspbpbre = oobibspbpbre, SURVNAOOB = SURVNAOOB,  SURVBREOOB  =  SURVBREOOB, BSTKM = BSTKM,
                                 BSTNAOOB = BSTNAOOB, BSTBREOOB = BSTBREOOB) else NULL), Timediff = Timediff, 
            TEST = (if (length(new_data) != 0) list(IBSNAKMnew = IBSNAKMnew, IBSBRKMnew = IBSBRKMnew, IBSKMnew = IBSKMnew,
                                                    SURVNAnew = SURVNAnew, SURVBREnew = SURVBREnew, SURV_NAnew = SURV_NAnew, SURV_BREnew = SURV_BREnew) else NULL)));

}
