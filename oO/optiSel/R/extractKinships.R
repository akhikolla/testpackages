
"extractKinships"<- function(obj, phen, BreedNames){
  ### Check list obj which may include kinships and native kinships ####
  ### Kinships are converted to class quadFun ##########################
  
  for(i in seq_along(obj)){
    if(!("quadFun" %in% class(obj[[i]])) && !("ratioFun" %in% class(obj[[i]])) && (!is.matrix(obj[[i]]))){
      stop("All '...' arguments must have class 'quadFun', 'ratioFun', or 'matrix'.\n")
    }
    
    Name <- names(obj)[i]
    if(is.null(Name)){stop(paste("A name must be specified for all additional parameters.\n",sep=""))}
    
    if(("quadFun" %in% class(obj[[i]])) || ("ratioFun" %in% class(obj[[i]]))){
      obj[[i]]$name <- Name
    }
    if(is.matrix(obj[[i]])){
      if(is.null(rownames(obj[[i]]))){
        stop(paste0("Rownames must be specified for matrix ", Name, ".\n",sep=""))
      }
      obj[[i]]   <- optiSolve::quadfun(Q=obj[[i]], d=0, id=rownames(obj[[i]]), name=Name)
    }
  }
  
  ### From kinship matrices containing kinships for multiple breeds  ###
  ### kinship matrices for single breeds are created and appended ######
  
  
  Seq <- seq_along(obj)
  for(i in Seq){
    if(("quadFun" %in% class(obj[[i]])) && all(phen$Indiv %in% obj[[i]]$id) && (length(BreedNames)>1)){
      Name <- obj[[i]]$name
      bname <- str_extract(Name, paste(BreedNames, collapse="|"))
      if(!is.na(bname)){stop(paste0("Kinship" , Name, " contains kinships from more than one breed, so the parameter name should not contain a breed name.\n"))}
      for(b in BreedNames){
        Kname <- paste(Name, b, sep=".")
        if(!(Kname %in% names(obj))){
          use <- obj[[i]]$id %in% phen$Indiv[phen$Breed==b]
          obj[[Kname]] <- optiSolve::quadfun(Q=obj[[i]]$Q[use, use], a=obj[[i]]$a[use],  d=0, id=obj[[i]]$id[use], name=Kname)
        }
      }
    }
  }
  
  
  ### Check if kinship matrices contain all required individuals ###
  ### Set component obj[[i]]$breed #################################
  for(i in seq_along(obj)){
    if(length(BreedNames)==1){
      ### Check if all individuals from 'phen' are included in ###############
      ### all kinship matrices ###############################################
      if(!all(phen$Indiv %in% obj[[i]]$id)){
        stop(paste0("Data frame 'phen' contains individuals not included in argument ", obj[[i]]$name, ".\n"))
      }
      obj[[i]]$breed <- BreedNames
    }else{
      ### Check if all kinship matrices contain either all individuals or ###
      ### all individuals from exactly one breed ############################
      if(all(phen$Indiv %in% obj[[i]]$id)){
        if("ratioFun" %in% class(obj[[i]])){stop("Kinships at native alleles must include individuals from only one breed.\n")}
        obj[[i]]$breed <- "across breeds"
      }else{
        thisBreed <- unique(phen[intersect(phen$Indiv, obj[[i]]$id),"Breed"])
        if(length(thisBreed)>1){
          stop(paste0("Argument ",obj[[i]]$name, " must include either all individuals from 'phen', or all individuals from exactly one breed.\n"))
        }
        if(!all(phen$Indiv[phen$Breed==thisBreed] %in% obj[[i]]$id)){
          stop(paste0("Some individuals from breed ", thisBreed, " are missing in argument", obj[[i]]$name, ".\n"))
        }
        obj[[i]]$breed <- thisBreed
      }
    }
  }
  
  
  ### Remove individuals not included in data frame 'phen' from kinships ##############
  ### Enlarge kinship matrices for one breed to include individuals from all breeds ###
  ### and sort individuals in kinship matrices according to 'phen' ####################
  
  for(i in seq_along(obj)){
    obj[[i]] <- adjust(obj[[i]], phen$Indiv)
    if("ratioFun" %in% class(obj[[i]])){
      obj[[i]]$NC <- obj[[i]]$NC[obj[[i]]$id]
    }
  }
  
  return(obj)
}