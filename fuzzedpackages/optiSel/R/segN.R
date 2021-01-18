

"segN"<-function(Native, map, unitP="Mb", keep=NULL, cores=1){
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
  checkMap(map, unitP, unitP)
  map$Name <- as.character(map$Name)
  map$Chr       <- as.character(map$Chr)
  rownames(map) <- map$Name
  if(!(unitP %in% c("SNP", colnames(map)))){
    stop("Parameter unitP must be either 'SNP' or the name of a column in the marker map.\n")
  }
  
  if(is.vector(Native)){
    names(Native)<-str_sub(str_extract(basename(Native), "Chr[0-9,a-z]*"),4,-1)
    if(!all(names(Native) %in% map$Chr)){
      stop("Some chromosomes are not in the map.")
    }
    if(any(duplicated(names(Native)))){
      stop("For some chromosomes different files were provided.")
    }
    IndivFileN <- scan(Native[1], nlines=1, what="character",quiet=TRUE)[-1]
  }else{
    IndivFileN <- colnames(Native)
  }
  if(is.null(keep)){
    indexN <- which(IndivFileN %in% IndivFileN)
  }else{
    indexN <- which(IndivFileN %in% keep)
  }
  NFileN <- length(IndivFileN)
  Indiv  <- IndivFileN[indexN]
  NC     <- length(indexN)
  
  if(!is.vector(Native)){
    map <- map[rownames(Native),]
  }
  
  Chromosomes <- unique(map$Chr)
  for(chr in Chromosomes){
    submap <- map[map$Chr==chr, ]
    M      <- nrow(submap)
    submap$SNP <- 1:M
    x    <- submap[, unitP]
    kb   <- (c(0,x)+c(x,x[length(x)]+x[1]))/2
    Nkb  <- kb[2:length(kb)]-kb[1:(length(kb)-1)]
    map$nUnits[map$Chr==chr] <- Nkb
  }
  GenL <- sum(map$nUnits)
  
  if(is.na(cores)){
    cores <- detectCores()-1
    if(is.na(cores)){cores <- 1}
    suppressWarnings(cLim <- floor((memory.limit()-memory.size()-2500)*1000*1000/(10*(2*NC^2))))
    if(!is.na(cLim)){cores <- min(cores,cLim)}
    cores <- max(cores, 1)
    cores <- min(cores, length(Chromosomes))
  }

  gc()
  if(cores==1){
    Res  <- matrix(0,nrow=NC, ncol=NC)
    for(chr in Chromosomes){
      if(is.vector(Native)){
        cat(paste("Reading chromosome ", chr, "...  "))
        Res <- Res + rcpp_segN(as.character(Native[chr]), as.integer(NFileN), as.integer(NC), as.integer(indexN-1), as.double(map$nUnits[map$Chr==chr]))
      }else{
        Units      <- matrix(map$nUnits[map$Chr==chr],ncol=NC, nrow=sum(map$Chr==chr), byrow=FALSE)
        thisNative <- Native[map$Chr==chr, indexN]
        Res        <- Res + t(thisNative)%*%(thisNative*Units)
      }
    }
  }else{
    cat(paste0("Using ",cores," cores... "))
    use_cor <- 1 + ((1:length(Chromosomes)) %% cores)
    Cores   <- unique(use_cor)
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    Res <- foreach(thisCore=Cores, .combine='+', .inorder=FALSE) %dopar% {
      x <- matrix(0,NC, NC)
      for(chr in Chromosomes[use_cor==thisCore]){
        if(is.vector(Native)){
          x <- x + rcpp_segN(as.character(Native[chr]), as.integer(NFileN), as.integer(NC), as.integer(indexN-1), as.double(map$nUnits[map$Chr==chr]))
        }else{
          Units      <- matrix(map$nUnits[map$Chr==chr],ncol=NC, nrow=sum(map$Chr==chr), byrow=FALSE)
          thisNative <- Native[map$Chr==chr, indexN]
          x <- x + t(thisNative)%*%(thisNative*Units)
        }
        gc()
      }
      x
    }
    stopCluster(cl)
    cat("finished.\n")
  }
  
  N   <- ncol(Res)/2
  Res <- ((Res[2*(1:N)-1,2*(1:N)-1] + Res[2*(1:N)-1,2*(1:N)] + Res[2*(1:N),2*(1:N)-1] + Res[2*(1:N),2*(1:N)])/GenL)/4
  rownames(Res) <- Indiv[2*(1:N)]
  colnames(Res) <- Indiv[2*(1:N)]
  Res
}