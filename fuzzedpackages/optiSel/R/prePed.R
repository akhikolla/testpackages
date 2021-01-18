

"prePed" <- function(Pedig, keep=NULL, thisBreed=NA, lastNative=NA, addNum=FALSE){
    PedigAsDataTable <- "data.table" %in% class(Pedig)
    Pedig <- as.data.frame(Pedig)
    if(PedigAsDataTable){setDF(Pedig)}
    colnames(Pedig)[1:3]<-c("Indiv", "Sire", "Dam")
    Pedig$Indiv <- as.character(Pedig$Indiv)
    Pedig$Sire  <- as.character(Pedig$Sire)
    Pedig$Dam   <- as.character(Pedig$Dam)
    Pedig$Sire[Pedig$Sire %in% c("","0"," ")] <- NA
    Pedig$Dam[Pedig$Dam  %in% c("","0"," ")] <- NA
    
    if(is.logical(keep)){
      keep <- Pedig$Indiv[keep]
      }
    if(!is.null(keep)){
      keep <- as.character(keep)
      keep <- setdiff(keep, c(NA, ""," ", "0"))
      }
    if(anyDuplicated(Pedig$Indiv)){
      cat("Duplicated IDs were removed.\n")
      ord   <- order(is.na(Pedig$Sire)+is.na(Pedig$Dam))
      Pedig <- Pedig[ord,]
      cat("This includes e.g.\n")
      print(head(Pedig[duplicated(Pedig$Indiv),]))      
      Pedig <- Pedig[!duplicated(Pedig$Indiv),]
    }
    rownames(Pedig)<-Pedig$Indiv
    
    if("Sex" %in% colnames(Pedig)){
      Pedig$Sex <- as.character(Pedig$Sex)
      Pedig$Sex[Pedig$Sex %in% c(""," ")] <- NA
    }else{
      Pedig$Sex<-NA
    }
    
    if("Breed" %in% colnames(Pedig)){
      if(!is.character(Pedig$Breed)){
        Pedig$Breed <- as.character(Pedig$Breed)
      }
      Pedig$Breed[Pedig$Breed %in% c(""," ")] <- NA
    }else{
      if(!is.na(thisBreed)){
        Pedig$Breed <- thisBreed
        }
    }
    if("Born" %in% colnames(Pedig) & !is.numeric(Pedig$Born)){
      Pedig$Born <- as.character(Pedig$Born)
      Pedig$Born[Pedig$Born %in% c(""," ")] <- NA
      Pedig$Born <- as.numeric(Pedig$Born)
    }
    withBreed <- ("Breed" %in% colnames(Pedig))
    withBorn  <- ("Born"  %in% colnames(Pedig))
    
    ######### Cut Pedigree loops #########
    suppressWarnings(ord<-pedigree::orderPed(Pedig[,1:3]))
    if(sum(ord==-1)>0){
      cat("Pedigree loops were detected. We recommend to correct them manually before\n")
      cat("using prePed(). The parents of the following individuals are set to unknown\n")
      cat("to remove the loops.\n")
      print(Pedig[ord==-1, 2:3])
      cat("\n")
    
      Pedig$Sire[ord==-1] <- NA
      Pedig$Dam[ord==-1]  <- NA
      Pedig$Breed[ord==-1]<- "Pedigree Error"
    }
 
    ####### Add imaginary ancestors #######
    if(!is.na(lastNative)){
      rownames(Pedig)<-Pedig$Indiv
      ID <-  Pedig$Indiv[is.na(Pedig$Sire) & !is.na(Pedig$Dam)]
      Pedig[ID, "Sire"] <- paste0("S", ID)
      ID <-  Pedig$Indiv[!is.na(Pedig$Sire) & is.na(Pedig$Dam)]
      Pedig[ID, "Dam"]  <- paste0("D", ID)
    }

    #### Add lines for ancestors, sort pedigree ####
    if(!all(is.na(Pedig$Dam) & is.na(Pedig$Sire))){
      Pedig <- nadiv::prepPed(Pedig)
    }
    Pedig$Indiv <- as.character(Pedig$Indiv)
    rownames(Pedig)<- Pedig$Indiv
    
    Mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    
    ### code sexes as male and females ####
    if(sum(!is.na(Pedig$Sex))>0){
      sexes<-names(table(Pedig$Sex))
      if(length(sexes)>2){
        cat("Warning: The sex has more than 2 levels. Please correct it.\n")
      }
      if(!all(sexes %in% c("male","female", NA))){
        if(length(sexes)==1){sexes<-c(sexes, "dummysex")}
        Mval <- Mode(Pedig$Sex[Pedig$Indiv %in% Pedig$Sire])
        Fval <- Mode(Pedig$Sex[Pedig$Indiv %in% Pedig$Dam])
        if(is.na(Mval)){Mval <- setdiff(sexes, Fval)}
        if(is.na(Fval)){Fval <- setdiff(sexes, Mval)}
        if(!is.na(Mval)&!is.na(Fval)&(Mval!=Fval)){
          Pedig$Sex <- mapvalues(Pedig$Sex, from=c(Mval, Fval), to=c("male","female"))
        }else{
        cat("Meaning of sex labels cannot be determined from pedigree structure.\n")
        }
      }
    }

    #### determine sexes from pedigree structure ####
    wrongMale   <- Pedig$Indiv %in% Pedig$Sire & !(Pedig$Sex %in% c("male", NA))
    wrongFemale <- Pedig$Indiv %in% Pedig$Dam  & !(Pedig$Sex %in% c("female", NA))
    if(sum(wrongMale)+sum(wrongFemale)>0){
      cat("The sex of the following animals was not compatible with the pedigree, so\n")
      cat("it was modified:\n")
      print(Pedig[wrongMale|wrongFemale, 2:3])
      cat("\n")
    }
    Pedig$Sex[Pedig$Indiv %in% Pedig$Sire] <- "male"
    Pedig$Sex[Pedig$Indiv %in% Pedig$Dam] <- "female"
    
    asSire  <- Pedig$Indiv %in% Pedig$Sire
    asDam   <- Pedig$Indiv %in% Pedig$Dam
    Dummies <- c(paste0("Founder", 1:20), paste0("Migrant", 1:20))
    if(any(asSire & asDam & !(Pedig$Indiv %in% Dummies))){
      cat("The following individuals appear as sire and as dam:\n")
      print(Pedig[asSire & asDam & !(Pedig$Indiv %in% Dummies), 1:3])
      stop("Please correct this")
    }
    
    ######  prune Pedigree   #######
    if(!is.null(keep)){
      Pedig <- nadiv::prunePed(Pedig, phenotyped=keep)
      Pedig$Indiv<- as.character(Pedig$Indiv)
      Pedig$Sire <- as.character(Pedig$Sire)
      Pedig$Dam  <- as.character(Pedig$Dam)
    }
    
    ####           Correct wrong breed names              ###
    if(withBreed & !is.na(thisBreed)){
      numSire <- match(Pedig$Sire, Pedig$Indiv)
      numDam <- match(Pedig$Dam, Pedig$Indiv)
      
      hasWrongBreed <- (!is.na(Pedig$Sire) & !is.na(Pedig$Dam) & !(Pedig$Breed %in% c(thisBreed, "Pedigree Error")) & (Pedig$Breed[numSire] %in% thisBreed) & (Pedig$Breed[numDam] %in% thisBreed))
      while(any(hasWrongBreed)){
        cat(paste0("The breed name of ",sum(hasWrongBreed)," individuals is changed to ",thisBreed))
        cat(paste0(" because both parents are ", thisBreed,".\nThis includes\n"))
        print(head(Pedig[hasWrongBreed,]))
        Pedig$Breed[hasWrongBreed] <- thisBreed
        hasWrongBreed <- (!is.na(Pedig$Sire) & !is.na(Pedig$Dam) & !(Pedig$Breed %in% c(thisBreed, "Pedigree Error")) & (Pedig$Breed[numSire] %in% thisBreed) & (Pedig$Breed[numDam] %in% thisBreed))
      }
    }
    
    
    ####             Estimate missing breeds              ###
    ###    Animals with missing breeds are assumed to     ###
    ### be from the same breed as most of their offspring ###
    if(withBreed){
      ID  <- Pedig$Indiv[is.na(Pedig$Breed)]
      Tab <- Pedig[(Pedig$Sire%in% ID | Pedig$Dam %in% ID) & !is.na(Pedig$Breed),c("Sire", "Dam", "Breed")] 
      Pedig[ID, "Breed"] <- tapply(c(Tab$Breed, Tab$Breed), list(c(Tab$Sire, Tab$Dam)), Mode)[ID]
    }
    
    ####   convert breed name of founders   ###
    ####  born after lastNative to unknown  ###
    if(withBreed & withBorn & !is.na(lastNative)){
      isFounder <- is.na(Pedig$Sire) & is.na(Pedig$Dam) & (Pedig$Breed==thisBreed | is.na(Pedig$Breed))
      notNative <- isFounder & Pedig$Born>lastNative
      names(notNative)<-Pedig$Indiv
      ID  <- Pedig$Indiv[is.na(notNative)]
      Tab <- Pedig[Pedig$Sire%in% ID | Pedig$Dam %in% ID & !is.na(Pedig$Born), c("Sire", "Dam", "Born")] 
      geb <- tapply(c(Tab$Born, Tab$Born), list(c(Tab$Sire, Tab$Dam)), min) - 0
      notNative[ID] <- geb[ID]>lastNative
      Pedig[is.na(notNative) | notNative, "Breed"] <- "unknown"
    }
    
    if(withBorn){
      ### Add generation intervall ##############################
      numSire <- match(Pedig$Sire, Pedig$Indiv)
      numDam  <- match(Pedig$Dam,  Pedig$Indiv)
      PBorn   <- cbind(Pedig$Born[numSire], Pedig$Born[numDam])
      PBorn   <- rowMeans(PBorn,  na.rm=TRUE)
      Pedig$I <- (Pedig$Born - PBorn)
      Pedig$I[Pedig$I<0] <- NA
    }
    
    Pedig$Offspring <- (Pedig$Indiv %in% Pedig$Sire)|(Pedig$Indiv %in% Pedig$Dam)

    ######    Add numeric IDs    #######
    if(addNum){
      nP<-nadiv::numPed(Pedig[,1:3])
      nP[nP==-998]<-0
      Pedig$numIndiv <- nP[,1]
      Pedig$numSire <- nP[,2]
      Pedig$numDam <- nP[,3]
    }
    
    cols <- c("Indiv", "Sire", "Dam", "Sex")
    if(withBreed){cols <- c(cols, "Breed")}
    if(withBorn){cols <- c(cols, "Born","I")}
    if(addNum){cols <- c(cols, c("numIndiv", "numSire", "numDam"))}
    cols  <- c(cols, setdiff(colnames(Pedig),  cols))
    Pedig <- Pedig[, cols]

    if(PedigAsDataTable){setDT(Pedig)}
    class(Pedig)<-c("Pedig", class(Pedig))
    Pedig
}