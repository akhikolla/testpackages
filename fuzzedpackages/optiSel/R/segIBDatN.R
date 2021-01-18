
"segIBDatN" <- function(files, phen, map, thisBreed, refBreeds="others", ubFreq=0.01, minSNP=20, unitP="Mb", minL=1.0, unitL="Mb", a=0.0, keep=NULL, lowMem=TRUE, skip=NA, cskip=NA, cores=1){
  ##################################################
  # Convert data tables to data frames             #
  ##################################################
  mapAsDataTable <- "data.table" %in% class(map)
  map <- as.data.frame(map)
  if(mapAsDataTable){setDF(map)}
  phenAsDataTable <- "data.table" %in% class(phen)
  phen <- as.data.frame(phen)
  if(mapAsDataTable){setDF(phen)}
  
  ##################################################
  #       Data preparation                         #
  ##################################################
  if(!("Indiv" %in% colnames(phen))){stop("Column 'Indiv' is missing in phen.")}
  if(!("Breed" %in% colnames(phen))){stop("Column 'Breed' is missing in phen.")}
  if(!("Name"  %in% colnames(map))){ stop("Column 'Name' is missing in map.")}
  if(!("Chr"   %in% colnames(map))){ stop("Column 'Chr'  is missing in map.")}
  if(("Sex" %in% colnames(phen)) && all(!is.na(phen$Sex))){
    if(all(phen$Sex=="male")){stop("Females are missing. Sexes can be ignored by removing column 'Sex'.\n")}
    if(all(phen$Sex=="female")){stop("Males are missing. Sexes can be ignored by removing column 'Sex'.\n")}
    if(!all(phen$Sex %in% c("female","male"))){stop("Sexes are not coded as 'female' and 'male'. They can be ignored by removing column 'Sex'.\n")}
  }
  phen$Indiv <- as.character(phen$Indiv)
  phen$Breed <- as.character(phen$Breed)
  if(is.logical(keep)){keep <- phen[keep,"Indiv"]}
  if(is.null(keep)){keep <- phen$Indiv}
  if(!is.null(keep)){
    keep <- as.character(keep)
    keep <- setdiff(keep, c(NA))
  }
  phen <- phen[!duplicated(phen$Indiv),]
  phen[is.na(phen$Breed),"Breed"] <- "notSpecified"
  rownames(phen)<- phen$Indiv
  
  map$Name      <- as.character(map$Name)
  map$Chr       <- as.character(map$Chr)
  rownames(map) <- map$Name
  if(is.character(files)){
    files <- list(hap.thisBreed=files)
  }
  if(is.na(skip)){  skip <- getskip(files$hap.thisBreed[1])}
  if(is.na(cskip)){cskip <- getcskip(files$hap.thisBreed[1], skip)}
  
  ##################################################
  #       Main part                                #
  ##################################################
  cat("Identifying native alleles...\n")
  if("match" %in% names(files)){
    Native <- files$match
  }else{
    if(lowMem){
      wdir <- file.path(tempdir(), "tempHapGFFdBvcw")
      dir.create(wdir, showWarnings=FALSE)
      Res <- haplofreq(files=files, phen=phen, map=map, thisBreed=thisBreed, refBreeds=refBreeds, minSNP=minSNP, minL=minL, unitL=unitL, ubFreq=ubFreq, keep=keep, skip=skip, cskip=cskip, what="match", w.dir=wdir, cores=cores)
      Native <- Res$match
    }else{
      Native <- haplofreq(files=files, phen=phen, map=map, thisBreed=thisBreed, refBreeds=refBreeds, minSNP=minSNP, minL=minL, unitL=unitL, keep=keep, skip=skip, cskip=cskip, what="freq", cores=cores)$freq < ubFreq
    }
  }
  
  cat("Computing probabilities for segments to be shared and native...\n")
  segIBDandN <- segIBDandN(files=files$hap.thisBreed, Native=Native, map=map, minSNP=minSNP, unitP=unitP, minL=minL, unitL=unitL, a=a, keep=keep, skip=skip, cskip=cskip, cores=cores)
  cat("Computing probabilities for segments to be native...\n")
  segN       <- segN(Native=Native, map=map, unitP=unitP, keep=keep, cores=cores)
  segN[segN==0] <- 1e-14
  
  
  
  Res <- optiSolve::ratiofun(Q1=segIBDandN, d1=0, Q2=segN, d2=0, id=rownames(segIBDandN))
  
  Res$mean <- mean(segIBDandN)/mean(segN)
  NC       <- segBreedComp(Native, map, unitP=unitP)$native
  Res$NC   <- setNames(NC, rownames(segIBDandN))
  
  #cat("Kinship at native alleles:", round(Res$mean,     4), "\n")
  #cat("Native Contribution:      ", round(mean(Res$NC), 4), "\n")
  
  if(lowMem & !("match" %in% names(files))){
    file.remove(Native)
    unlink(wdir, recursive=TRUE)
  }
  return(Res)
}
