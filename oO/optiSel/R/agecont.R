
### ageM=0 if the individual and its sire are born in the same year
### In this case, the individual is born when its sire is in age class 1.

"agecont" <- function(Pedig, use=Pedig$Born >= quantile(Pedig$Born, 0.75), maxAge=NA){
  Pedig <- checkphen(Pedig, columns=c("Indiv", "Born", "Sire", "Dam"), quiet=TRUE)

  Pedig$Born <- floor(Pedig$Born)
  
  Sire <- Pedig[use,"Sire"]
  Dam  <- Pedig[use,"Dam"]
  ageM <- Pedig[use,"Born"] - Pedig[Sire,"Born"]
  ageF <- Pedig[use,"Born"] - Pedig[Dam,"Born"]
  ageM <- ageM[!is.na(ageM) & ageM>=0]
  ageF <- ageF[!is.na(ageF) & ageF>=0]
  
  if(!is.na(maxAge)){
    ageM <- ageM[ageM<=maxAge]
    ageF <- ageF[ageF<=maxAge]
  }
  
  if(min(c(ageM, ageF))==0){
    cat("Individuals whose parents are born in the same year are removed.\n")
    ageM <- ageM[ageM>0]
    ageF <- ageF[ageF>0]
  }
  
  q <- max(c(ageM+1,ageF+1), na.rm=TRUE)
  
  #### Proportion of offspring born when the individual ##########
  #### is in the respective age class ############################
  
  pOff <- data.frame(age  = 1:q, 
                   male   = tabulate(ageM+1, nbins=q)/length(ageM), 
                   female = tabulate(ageF+1, nbins=q)/length(ageF))

  #### Generation interval ######################################
  
  Lm <- sum((pOff$age-1)*pOff$male)   #Lf = 1 + sum(v$female)
  Lf <- sum((pOff$age-1)*pOff$female) #Lm = 1 + sum(v$male)
  #L  <- (Lm + Lf)/2
  
  alpham <- ((Lf-1)/Lf)/((Lm-1)/Lm + (Lf-1)/Lf)
  alphaf <- 1-alpham
  
  #### Proportion of offspring not yet born when the individual ####
  #### enters the respective age class ############################x
  
  vm <- rev(cumsum(rev(pOff$male[-1]  )))
  vf <- rev(cumsum(rev(pOff$female[-1])))
  
  v <- data.frame(age    = pOff$age, 
                  male   = c(vm,0), 
                  female = c(vf,0))
  
  data.frame(age = 1:q, male = alpham*v$male/Lm, female = alphaf*v$female/Lf)
}
