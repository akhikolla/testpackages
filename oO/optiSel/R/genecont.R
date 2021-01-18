

"genecont" <- function(Pedig, from=NULL, to=NULL){
  PedigAsDataTable <- "data.table" %in% class(Pedig)
 
  if(PedigAsDataTable){
    Pedig <- as.data.frame(Pedig)
    setDF(Pedig)
    }

  colnames(Pedig)[1:3] <- c("Indiv", "Sire", "Dam")
  if(is.null(from)){
    hasOffspring <- (Pedig$Indiv %in% Pedig$Sire)|(Pedig$Indiv %in% Pedig$Dam)
    from <- Pedig$Indiv[hasOffspring]
  }else{
    from <- setdiff(as.character(from), c(NA, "", " ", "0"))
  }
  
  if(is.null(to)){
    to <- Pedig$Indiv
  }else{
    to <- setdiff(as.character(to), c(NA, "", " ", "0"))
  }  
  
  Pedig <- prePed(Pedig, addNum=TRUE, keep=c(from, to))
  Pedig$nOff <- 0
  x <- table(Pedig$Sire)
  Pedig[names(x),"nOff"] <- x
  x <- table(Pedig$Dam)
  Pedig[names(x),"nOff"] <- x
  
  from  <- Pedig$Indiv[Pedig$Indiv %in% from] 
  to    <- Pedig$Indiv[Pedig$Indiv %in% to]
  numAnc  <- match(from, Pedig$Indiv)
  numKeep <- match(to,   Pedig$Indiv)
  rNames  <- as.character(Pedig$Indiv[Pedig$Indiv %in% to])
  cNames  <- from

  GeneCont <- rcpp_genecont(as.integer(Pedig$numSire), as.integer(Pedig$numDam), as.integer(numAnc-1), as.integer(numKeep-1), as.integer(Pedig$Indiv %in% to), rNames, cNames, as.integer(Pedig$nOff))
  if(identical(rNames, to)){return(GeneCont)}
  GeneCont[to, from, drop=FALSE]
}





