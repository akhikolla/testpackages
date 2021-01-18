
"haplofreq"<-function(files, phen, map, thisBreed, refBreeds="others", minSNP=20, minL=1.0, unitL="Mb", ubFreq=0.01, keep=NULL, skip=NA, cskip=NA, w.dir=NA, what=c("freq","match"), cores=1){
  ##################################################
  # Convert data tables to data frames             #
  ##################################################
  mapAsDataTable <- "data.table" %in% class(map)
  map <- as.data.frame(map)
  if(mapAsDataTable){setDF(map)}
  phenAsDataTable <- "data.table" %in% class(phen)
  phen <- as.data.frame(phen)
  if(phenAsDataTable){setDF(phen)}
  ##################################################
  #       Data preparation                         #
  ##################################################
  if(!("Indiv" %in% colnames(phen))){stop("Column 'Indiv' is missing in phen.")}
  if(!("Breed" %in% colnames(phen))){stop("Column 'Breed' is missing in phen.")}
  if(!("Name"  %in% colnames(map))){ stop("Column 'Name' is missing in map.")}
  if(!("Chr"   %in% colnames(map))){ stop("Column 'Chr'  is missing in map.")}
  phen$Indiv <- as.character(phen$Indiv)
  phen$Breed <- as.character(phen$Breed)
  phen <- phen[!duplicated(phen$Indiv),]
  phen[is.na(phen$Breed),"Breed"] <- "notSpecified"
  map$Name   <- as.character(map$Name)
  map$Chr <- as.character(map$Chr)
  if(unitL=="SNP"){
    map$SNP    <- NA
    for(chr in unique(map$Chr)){
      map[map$Chr==chr, "SNP"] <- (1:sum(map$Chr==chr))
    }
  }
  if(!(unitL   %in% colnames(map))){ stop(paste("Column '", unitL, "'  is missing in map."))}
  rownames(phen)<- phen$Indiv
  rownames(map) <- map$Name
  checkMap(map, unitL, unitL)
  
  if(is.logical(keep)){
    keep <- phen$Indiv[!is.na(keep) & keep]
    }
  if(!is.null(keep)){
    keep <- setdiff(as.character(keep),c(NA))
    if(!all(keep %in% phen$Indiv)){
      cat("The following Individuals are removed because they do not appear in phen:")
      cat(setdiff(keep, phen$Indiv),"\n")
      keep <- keep[keep %in% phen$Indiv]
    }
  }else{
    keep <- phen$Indiv
  }
  
  ##################################################
  # Specify files for reading and writing          #
  ##################################################
  if(!is.list(files)){files<-list(hap.thisBreed=files)}
  if(!("hap.thisBreed" %in% names(files))){stop("List of files must have component 'hap.thisBreed'\n")}
  if(!("hap.refBreeds" %in% names(files))){files$hap.refBreeds<-files$hap.thisBreed}
  names(files$hap.thisBreed) <- str_sub(str_extract(basename(files$hap.thisBreed), "Chr[0-9,a-z]*"),4,-1)
  names(files$hap.refBreeds) <- str_sub(str_extract(basename(files$hap.refBreeds), "Chr[0-9,a-z]*"),4,-1)
  if(is.na(w.dir)){
    files$seg.origin    <- rep("",length(files$hap.thisBreed))
    files$seg.frequency <- rep("",length(files$hap.thisBreed))
  }else{
    if("match" %in% what){
      files$seg.origin    <- paste("match",thisBreed,paste("Chr", names(files$hap.thisBreed),sep=""),"in",paste(refBreeds, collapse="."),"txt",sep=".")
      files$seg.origin    <- file.path(w.dir, files$seg.origin)
    }else{
      files$seg.origin    <- rep("",length(files$hap.thisBreed))
    }
    if("freq" %in% what){
      files$seg.frequency <- paste("freq", thisBreed,paste("Chr", names(files$hap.thisBreed),sep=""),"in",paste(refBreeds, collapse="."),"txt",sep=".")
      files$seg.frequency <- file.path(w.dir, files$seg.frequency)
    }else{
      files$seg.frequency    <- rep("",length(files$hap.thisBreed))
    }
  }
  names(files$seg.origin)    <- names(files$hap.thisBreed)
  names(files$seg.frequency) <- names(files$hap.thisBreed)
  save.frequency   <- !is.na(w.dir) && ("freq" %in% what)
  save.origin      <- !is.na(w.dir) && ("match" %in% what)
  return.result    <- !save.frequency && !save.origin
  return.frequency <- is.na(w.dir) && ("freq" %in% what)
  return.origin    <- is.na(w.dir) && ("match" %in% what)
  
  if(!is.na(w.dir)){
      dir.create(w.dir, showWarnings=FALSE, recursive=TRUE)
  }
  if(!identical(sort(names(files$hap.thisBreed)), sort(names(files$hap.refBreeds)))){
    stop("Chromosome names do not match for hap.thisBreed and hap.refBreeds.\n")
  }
  if(!all(names(files$hap.thisBreed) %in% map$Chr)){
    stop("Some chromosomes are not in the map.")
  }
  if(any(duplicated(names(files$hap.thisBreed)))){
    stop("For some chromosomes different files were provided.")
  }
  if(is.na(skip)){  skip <- getskip(files$hap.thisBreed[1])}
  if(is.na(cskip)){cskip <- getcskip(files$hap.thisBreed[1], skip)}
  ##################################################
  #      Extract information for this breed        #
  ##################################################
  IndivInFileC <- scan(files$hap.thisBreed[1], nlines=1, what="character",quiet=TRUE, skip=skip)
  if(cskip>0){IndivInFileC <-IndivInFileC[-(1:cskip)]}
  if(length(IndivInFileC)==1 || IndivInFileC[1]!=IndivInFileC[2]){
    IndivInFileC <- rep(IndivInFileC, each=2)
  }
  IndexCand    <- which((IndivInFileC %in% keep) & (phen[IndivInFileC,"Breed"]==thisBreed))
  Candidate    <- IndivInFileC[IndexCand]
  NC           <- length(Candidate)
  NFileC       <- length(IndivInFileC)
  if(NC==0){stop("Individuals from this breed are not included in phen or keep.\n")}

  remainingIndiv <- setdiff(keep, Candidate)
  if(any(phen$Indiv %in% remainingIndiv & phen$Breed==thisBreed)){
    cat("The following individuals are removed from the analysis because they are not genotyped:\n")
    print(phen[phen$Indiv %in% remainingIndiv & phen$Breed==thisBreed,])
  }
  ##################################################
  #     Extract information for reference breeds   #
  ##################################################
  IndivInFileR <- scan(files$hap.refBreeds[1], nlines=1, what="character",quiet=TRUE, skip=skip)
  if(cskip>0){IndivInFileR <-IndivInFileR[-(1:cskip)]}
  if(length(IndivInFileR)==1 || IndivInFileR[1]!=IndivInFileR[2]){
    IndivInFileR <- rep(IndivInFileR, each=2)
  }
  NFileR       <- length(IndivInFileR)
  BreedR       <- phen[IndivInFileR,"Breed"]
  keptBreedR   <- BreedR[IndivInFileR %in% keep]
  if(identical(refBreeds, "all")){refBreeds <- unique(keptBreedR)}
  if(identical(refBreeds, "others")){refBreeds <- setdiff(unique(keptBreedR),thisBreed)}
  BreedSymbol <- paste(substring(refBreeds, 1, 1), collapse="")

  NR <- table(keptBreedR)[refBreeds]
  if(length(NR)==0||any(is.na(NR))||max(NR)==0){stop(paste("Individuals from reference breeds are not included in phen or keep.","\n"))}
  
  IndexRef <- matrix(NA, nrow=max(NR), ncol=length(refBreeds), dimnames=list(NULL,refBreeds))
  RefIndiv <- NULL
  for(b in refBreeds){
    if(NR[b]==0){stop(paste("Individuals from reference breed",b,"are not included in phen or keep.","\n"))}
    IndexRef[1:NR[b],b] <- which((IndivInFileR %in% keep) & (phen[IndivInFileR,"Breed"] %in% b))-1
    RefIndiv <- c(RefIndiv, IndivInFileR[IndexRef[1:NR[b],b]])
    }
  storage.mode(IndexRef)<-'integer'

  remainingIndiv <-setdiff(keep, RefIndiv)
  if(any(phen$Indiv %in% remainingIndiv & phen$Breed %in% refBreeds)){
    cat("The following reference individuals are removed from the analysis because they are not genotyped:\n")
    print(phen[phen$Indiv %in% remainingIndiv & phen$Breed %in% refBreeds,])
  }  
  ##################################################
  #    Compute segment frequencies and origin      #
  ##################################################
  symB <- scan(files$hap.thisBreed[1], n=cskip+1, skip=1+skip, what="character", quiet=TRUE)[cskip+1]
  
  Chromosomes <- unique(map$Chr)
  if(save.frequency){
    for(chr in Chromosomes){
      writeLines(paste(c("Name", Candidate), collapse=" "), files$seg.frequency[chr])
    }
  }
  if(save.origin){
    for(chr in Chromosomes){
      writeLines(paste(c("Name",Candidate), collapse=" "), files$seg.origin[chr])
    }
  }

  cMList <- list()
  MNamesList <- list()
  for(chr in Chromosomes){
    submap <- map[map$Chr==chr, ]
    MNamesList[[chr]] <- str_pad(submap$Name,max(nchar(submap$Name))+1, side="right")
    cM <- submap[, unitL]
    cMList[[chr]] <- (c(0,cM)+c(cM,cM[length(cM)]+cM[1]))/2
  }

  if(is.na(cores)){
    cores <- detectCores()-1
    if(is.na(cores)){cores <- 1}
    cores <- max(cores, 1)
    cores <- min(cores, length(Chromosomes))
  }
  
  if(cores==1){  
    Haplo <- list()
    for(chr in Chromosomes){
      cat(paste("Reading chromosome ", chr, "...  "))
      Haplo[[chr]] <- rcpp_haplofreq(as.character(files$hap.thisBreed[chr]), as.character(files$hap.refBreeds[chr]), as.character(files$seg.frequency[chr]), as.character(files$seg.origin[chr]), as.character(MNamesList[[chr]]), as.character(BreedSymbol), as.integer(IndexCand-1), IndexRef, as.integer(NFileC), as.integer(NFileR), as.integer(NC), as.integer(NR), as.integer(minSNP), as.double(minL), as.double(ubFreq), as.double(cMList[[chr]]), as.character(symB), as.integer(skip), as.integer(cskip), as.integer("freq" %in% what), as.integer("match" %in% what))
      cat(paste("M=", length(cMList[[chr]])-1, ", Reference breeds: ", paste(refBreeds, collapse=" "), "\n",sep=""))
      }
  }else{
    cat(paste0("Using ",cores," cores... "))
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    Haplo <- foreach(chr=Chromosomes, .inorder=TRUE) %dopar% {
      rcpp_haplofreq(as.character(files$hap.thisBreed[chr]), as.character(files$hap.refBreeds[chr]), as.character(files$seg.frequency[chr]), as.character(files$seg.origin[chr]), as.character(MNamesList[[chr]]), as.character(BreedSymbol), as.integer(IndexCand-1), IndexRef, as.integer(NFileC), as.integer(NFileR), as.integer(NC), as.integer(NR), as.integer(minSNP), as.double(minL), as.double(ubFreq), as.double(cMList[[chr]]), as.character(symB), as.integer(skip), as.integer(cskip), as.integer("freq" %in% what), as.integer("match" %in% what))
      }
    names(Haplo) <- Chromosomes
    stopCluster(cl)
    cat("finished.\n")
  }
  

  if(return.result){
    Res <- list()
    if(return.frequency){
      Freq  <- NULL
      for(chr in Chromosomes){
        rownames(Haplo[[chr]]$freq) <- map$Name[map$Chr==chr]
        Freq <- rbind(Freq, Haplo[[chr]]$freq)
      }
      colnames(Freq) <- Candidate
      Res$freq <- Freq
      attributes(Res)$map <- map[rownames(Freq),]
    }
    if(return.origin){
      Match <- NULL
      for(chr in Chromosomes){
        rownames(Haplo[[chr]]$match) <- map$Name[map$Chr==chr]
        Match <- rbind(Match, Haplo[[chr]]$match)
      }
      colnames(Match) <- Candidate
      Res$match <- Match
      attributes(Res)$map <- map[rownames(Match),]
    }
    
    attributes(Res)$Indiv     <- Candidate
    attributes(Res)$thisBreed <- thisBreed
    attributes(Res)$refBreeds <- refBreeds
    attributes(Res)$minSNP    <- minSNP
    attributes(Res)$minL      <- minL
    attributes(Res)$unitL     <- unitL
    class(Res)<-"HaploFreq"
    return(Res)
  }else{
    Res <- data.frame(Chr=names(files$hap.thisBreed), stringsAsFactors = FALSE)
    if("freq" %in% what){Res$freq   <- files$seg.frequency}
    if("match" %in% what){Res$match <- files$seg.origin}
    return(Res)
  }
  
}
