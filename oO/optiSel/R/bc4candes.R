
"bc4candes" <- function(kinship, phen, BreedNames, bc, quiet){
  if(length(BreedNames)==1){
    return(setNames(c(1), BreedNames))
  }
  
  if(!is.null(bc) && !is.character(bc)){
    return(checkbc(BreedNames, bc))
  }
  
  bcKin <- ifelse(is.character(bc), bc, NA)
  
  ### Get kinship bcKin for estimating the optimum breed composition ##

  if(is.character(bcKin)){
    if(!(bcKin %in% names(kinship))){
      stop(paste0("bc=", bcKin, " is not a valid name of a kinship.\n"))
    }
    if(!("quadFun" %in% class(kinship[[bcKin]]))){
      stop(paste0("bc=", bcKin, " is not a kinship.\n"))
    }
    if(kinship[[bcKin]]$breed != "across breeds"){
      stop(paste0("Kinship ",bcKin," must contain kinships from all breeds.\n"))
    }
  }else{
    iKin <- NA
    for(i in rev(seq_along(kinship))){
      if(("quadFun" %in% class(kinship[[i]])) && (kinship[[i]]$breed == "across breeds")){
        iKin <- i
      }
    }
    if(is.na(iKin)){
      stop("No kinship is suitable for estimating the optimum breed composition. Please provide parameter 'bc',\n")
    }
    bcKin <- names(kinship)[iKin]
    if(!quiet){
      cat(paste0("Breed contributions minimize ",bcKin," across breeds.\n"))
    }
  }
  
  lbbc <- setNames(rep(0.10/length(BreedNames),length(BreedNames)), BreedNames)
  bc   <- opticomp(kinship[[bcKin]]$Q, phen, lb=lbbc)$bc
  bc   <- checkbc(BreedNames, bc)
 
  return(bc)
}