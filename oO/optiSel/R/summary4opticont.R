
"summary4opticont" <- function(sy, x, cand, condf){
  BreedNames   <- names(cand$breed)
  KinshipNames <- names(cand$kinship)
  TraitNames   <- gettraits(cand$phen)
  
  sy$Breed <- setNames(cand$current$Breed, cand$current$Var)[sy$Var]
  sy$Breed[is.na(sy$Breed)] <- str_extract(sy$Var, paste(BreedNames, collapse="|"))[is.na(sy$Breed)]
  sy$Breed[is.na(sy$Breed)] <- "" 
  
  sy$Name <- str_replace(sy$Var, "bc.", "breed contribution.")
  sy$Name <- str_replace(sy$Name,"scd.", "sex contrib. diff..")
  sy$Name <- str_replace(sy$Name, paste(paste0("\\.",BreedNames), collapse="|"),"")
  
  ####### Add original threshold values to summary ######## 

  use <- sy$Var %in% condf$var
  sy[use, "Bound"] <- condf[sy$Var[use], "val"]
  
  ####### Compute constraint values for summary ######## 
  
  for(i in 1:nrow(sy)){
    var    <- sy$Var[i]
    b      <- sy$Breed[i]
    thisbc <- ifelse(b == "across breeds", sum(x), sum(x[cand$phen$Breed==b]))
    
    if(var %in% TraitNames){
      sy[i,"Val"] <- sum((cand$phen[[var]])*x, na.rm=TRUE)/thisbc
    }
    
    if(var %in% KinshipNames){
      if(class(cand$kinship[[var]])=="quadFun"){
        sy[i,"Val"] <- c(t(x)%*%(cand$kinship[[var]]$Q)%*%x/(thisbc^2)) 
        sy[i,"Val"] <- sy[i,"Val"] + sum(x*cand$kinship[[var]]$a)/thisbc + cand$kinship[[var]]$d
      }
    }
  }
  
  return(sy)
}