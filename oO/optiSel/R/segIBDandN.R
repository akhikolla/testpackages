globalVariables('thisCore')

"segIBDandN"<-function(files, Native, map, minSNP=20, unitP="Mb", minL=1.0, unitL="Mb", a=0.0, keep=NULL, skip=NA, cskip=NA, cores=1){
  ##################################################
  # Convert data tables to data frames             #
  ##################################################
  mapAsDataTable <- "data.table" %in% class(map)
  map <- as.data.frame(map)
  if(mapAsDataTable){setDF(map)}
  
  ##################################################
  #       Data preparation                         #
  ##################################################
  if(!("Name"  %in% colnames(map))){ stop("Column 'Name' is missing in map.")}
  if(!("Chr"   %in% colnames(map))){ stop("Column 'Chr'  is missing in map.")}  
  if(!is.null(keep)){
    keep <- as.character(keep)
    keep <- setdiff(keep, c(NA))
  }
  checkMap(map, unitL, unitP)
  map$Name      <- as.character(map$Name)
  map$Chr       <- as.character(map$Chr)
  rownames(map) <- map$Name

  names(files)<-str_sub(str_extract(basename(files), "Chr[0-9,a-z]*"),4,-1)
  if(!all(names(files) %in% map$Chr)){
    stop("Some chromosomes are not in the map.")
  }
  if(any(duplicated(names(files)))){
    stop("For some chromosomes different files were provided.")
  }
  ##################################################
  #       Main part                                #
  ##################################################
  
  if(is.vector(Native)){
    names(Native)<-str_sub(str_extract(basename(Native), "Chr[0-9,a-z]*"),4,-1)
    if(!identical(sort(names(files)), sort(names(Native)))){
      stop("Chromosome names do not match for files and Native.\n")
    }
    IndivFileN <- scan(Native[1], nlines=1, what="character",quiet=TRUE)[-1]
  }else{
    IndivFileN <- colnames(Native)
  }
  if(is.na(skip)){  skip <- getskip(files[1])}
  if(is.na(cskip)){cskip <- getcskip(files[1], skip)}
  IndivFileC <- scan(files[1], nlines=1, what="character",quiet=TRUE, skip=skip)
  if(cskip>0){IndivFileC <-IndivFileC[-(1:cskip)]}
  if(length(IndivFileC)==1 || IndivFileC[1]!=IndivFileC[2]){
    IndivFileC <- rep(IndivFileC, each=2)
  }
  Indiv      <- IndivFileC[IndivFileC %in% IndivFileN]
  if(any(Indiv != IndivFileN[IndivFileN %in% IndivFileC])){
    stop("Individuals must be in the same order in all files and matrices.\n")
  }
  if(!is.null(keep)){
    Indiv <- Indiv[Indiv %in% keep]
  }
  indexC     <- which(IndivFileC %in% Indiv)
  indexN     <- which(IndivFileN %in% Indiv)
  NFileC     <- length(IndivFileC)
  NFileN     <- length(IndivFileN)
  NC         <- length(indexC)
  symB       <- scan(files[1], n=cskip+1, skip=1+skip, what="character", quiet=TRUE)[cskip+1]
  
  Chromosomes <- names(files)
  gesL   <- 0
  cMList <- list()
  kbList <- list()
  for(chr in Chromosomes){
    submap <- map[map$Chr==chr, ]
    M <- nrow(submap)
    submap$SNP <- 1:M
    if(unitL %in% colnames(submap)){cM <- submap[, unitL]}else{cM <- 1:M; cat("Using: unitL=Marker number\n")}
    if(unitP %in% colnames(submap)){kb <- submap[, unitP]}else{kb <- 1:M; cat("Using: unitP=Marker number\n")}
    cM   <- (c(0,cM)+c(cM,cM[length(cM)]+cM[1]))/2
    kb   <- (c(0,kb)+c(kb,kb[length(kb)]+kb[1]))/2
    gesL <- gesL + kb[length(kb)] - kb[1]
    cMList[[chr]] <- cM
    kbList[[chr]] <- kb
  }
  
  if(is.na(cores)){
    M <- max(table(map$Chr))
    cores <- detectCores()-1
    if(is.na(cores)){cores <- 1}
    if(is.vector(Native)){
      suppressWarnings(cLim <- floor((memory.limit()-memory.size()-2500)*1000*1000/(10*(3*NC^2/2+1*NC^2+1*M*NC))))
      if(!is.na(cLim)){cores <- min(cores,cLim)}
    }else{
      suppressWarnings(cLim <- floor((memory.limit()-memory.size()-2500)*1000*1000/(10*(2*NC^2/2+ 1*NC^2+2*M*NC))))
      if(!is.na(cLim)){cores <- min(cores,cLim)}
      }
    cores <- max(cores, 1)
    cores <- min(cores, length(Chromosomes))
  }

  gc()
  if(cores==1){
    fROHN <- matrix(0,NC,NC)
    for(chr in Chromosomes){
      cat(paste0("Reading chromosome ", chr, "...  "))
      if(is.vector(Native)){
        fROHN <- fROHN +  rcpp_segIBDandN(as.character(files[chr]), as.character(Native[chr]), as.integer(NFileC), as.integer(NFileN), as.integer(indexC-1), as.integer(indexN-1), as.integer(NC), as.integer(minSNP), as.double(minL), as.double(cMList[[chr]]), as.double(kbList[[chr]]), as.double(a), as.character(symB), as.integer(skip), as.integer(cskip))
      }else{
      fROHN <- fROHN +  rcpp_segIBDandNVersion2(as.character(files[chr]), as.integer(NFileC), as.integer(NC), as.integer(indexC-1), Native[map$Name[map$Chr==chr],indexN], as.integer(minSNP), as.double(minL), as.double(cMList[[chr]]), as.double(kbList[[chr]]), as.double(a), as.character(symB), as.integer(skip), as.integer(cskip))
      }
    }
  }else{
    use_cor <- 1 + ((1:length(Chromosomes)) %% cores)
    Cores   <- unique(use_cor)
    cat(paste0("Using ",cores," cores... "))
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    fROHN <- foreach(thisCore=Cores, .combine='+', .inorder=FALSE) %dopar% {
        x <- matrix(0,NC, NC)
        for(chr in Chromosomes[use_cor==thisCore]){
          if(is.vector(Native)){
            x <- x + rcpp_segIBDandN(as.character(files[chr]), as.character(Native[chr]), as.integer(NFileC), as.integer(NFileN), as.integer(indexC-1), as.integer(indexN-1), as.integer(NC), as.integer(minSNP), as.double(minL), as.double(cMList[[chr]]), as.double(kbList[[chr]]), as.double(a), as.character(symB), as.integer(skip), as.integer(cskip))
          }else{
            x <- x + rcpp_segIBDandNVersion2(as.character(files[chr]), as.integer(NFileC), as.integer(NC), as.integer(indexC-1), Native[map$Name[map$Chr==chr],indexN], as.integer(minSNP), as.double(minL), as.double(cMList[[chr]]), as.double(kbList[[chr]]), as.double(a), as.character(symB), as.integer(skip), as.integer(cskip))
          }
          gc()
        }
        x
    }
    gc()
    stopCluster(cl)
    cat("finished.\n")
  }

  N     <- nrow(fROHN)/2
  fROHN <- (fROHN[2*(1:N)-1,2*(1:N)-1]+ fROHN[2*(1:N)-1,2*(1:N)]+ fROHN[2*(1:N),2*(1:N)-1]+ fROHN[2*(1:N),2*(1:N)])/4
  rownames(fROHN) <- IndivFileC[indexC][2*(1:N)]
  colnames(fROHN) <- IndivFileC[indexC][2*(1:N)]
  fROHN/gesL
}


"checkMap"<-function(map, unitL, unitP){
  if( is.null(map) && (unitP!="SNP" | unitL!="SNP")){
    stop("No marker map was provided, so only SNP can be used to measure coverage and distance.")
  }
  if(!is.null(map)){
    if(!("Name" %in% colnames(map))){
      stop("The marker map must contain column 'Name' with marker names.")
    }
    if(!(unitP %in% c(colnames(map), "SNP"))){
      stop("Variable 'unitP' must be 'SNP' or the name of a column in the marker map.")
    }
    if(!(unitL %in% c(colnames(map), "SNP"))){
      stop("Variable 'unitL' must be 'SNP' or the name of a column in the marker map.")
    }
    if(unitL %in% colnames(map) && !is.numeric(map[,unitL])){
      stop(paste("Column", unitL, "must be numeric.",sep=" "))
    }
    if(unitP %in% colnames(map) && !is.numeric(map[,unitP])){
      stop(paste("Column", unitP, "must be numeric.",sep=" "))
    }
    if(unitL %in% colnames(map) && sum(as.numeric(tapply(map[,unitL],map$Chr,max)))>2^31-1){
      stop(paste("Integer overflow: Column", unitL, "must have smaller values.",sep=" "))
    }
    if(unitP %in% colnames(map) && sum(as.numeric(tapply(map[,unitP],map$Chr,max)))>2^31-1){
      stop(paste("Integer overflow: Column", unitP, "must have smaller values.",sep=" "))
    }
  }  
}