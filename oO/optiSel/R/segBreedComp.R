
"segBreedComp"<-function(Native, map, unitP="Mb"){
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
  map$Name <- as.character(map$Name)
  map$Chr  <- as.character(map$Chr)  
  checkMap(map, unitP, unitP)
  rownames(map) <- map$Name
 
  if(!(unitP %in% c("SNP", colnames(map)))){
    stop("Parameter unitP must be either 'SNP' or the name of a column in the marker map.\n")
  }
  
  if(unitP=="SNP"){
    map$SNP <- NA
    for(chr in unique(map$Chr)){
      map[map$Chr==chr, "SNP"] <- (1:sum(map$Chr==chr))
    }
  }
  
  map$nUnits<-NA
  for(chr in unique(map$Chr)){
    x  <- map[map$Chr==chr, unitP]
    kb <- (c(0,x)+c(x,x[length(x)]+x[1]))/2
    map[map$Chr==chr, "nUnits"] <- kb[2:length(kb)]-kb[1:(length(kb)-1)]
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
  indexN <- which(IndivFileN %in% IndivFileN)
  NFileN <- length(IndivFileN)
  Indiv  <- IndivFileN[indexN]
  NC     <- length(indexN)
  NCont  <- rep(0, NC)
  BCont <- data.frame(native=rep(0,length(Indiv)))
  
  if(is.vector(Native)){
    mymap<-NULL
    M <- rep(NA, length(names(Native)))
    for(i in 1:length(Native)){
      chr   <- names(Native)[i]
      mymap <- rbind(mymap, map[map$Chr==chr, ])
      M[i]  <- sum(map$Chr==chr)
    }
    BCont <- rcpp_segBreedComp(as.character(Native), as.integer(NFileN), as.integer(NC), as.integer(indexN-1), as.integer(M), as.double(mymap$nUnits))
    sym   <- c("native", sort(setdiff(colnames(BCont), "native")))
    BCont <- BCont[ ,sym]
  }else{
    map    <- map[rownames(Native),]
    Units  <- matrix(map$nUnits,ncol=NC, nrow=nrow(Native), byrow=FALSE)
    Native <- Native[, indexN]
    if(is.character(Native)){
      sym     <- sort(setdiff(c(Native),"1"))
      BSymbol <- c("1",  sym)
      BName   <- c("native", sym)
      for(b in 1:length(BSymbol)){
        BCont[[BName[b]]] <- apply((Native==BSymbol[[b]])*Units,2,sum)
      }
    }else{
      BCont$native  <- apply(Native*Units, 2, sum)
    }
  }
  N     <- nrow(BCont)/2
  BCont <- (BCont[2*(1:N),,drop=FALSE] + BCont[2*(1:N)-1,,drop=FALSE])/(2*sum(map$nUnits))
  BCont <- data.frame(Indiv=Indiv[2*(1:N)], BCont, stringsAsFactors = FALSE)
  rownames(BCont)<- Indiv[2*(1:N)]
  BCont
}