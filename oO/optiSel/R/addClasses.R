
###### This function appends ######################################################
### Column 'Age' with age classes 1,2,..., of all individuals                   ###
### Column 'Class' with Breed X Age or Breed X Age X Sex classes                ###
###################################################################################
### Note:                                                                       ###
###  - Time t denotes the current generation or the current year                ###
###  - Individuals born at time t are in age class 1                            ###
###  - For OCS with non-overlapping generations, column 'Born' must contain the ###
###    generation number, not the year-of-birth.                                ###
###  - A breed has BreedxAge classes if there are NA-values in column 'Sex',    ###
###    and BreedxAgexSex classes otherwise.                                     ###
###  - All individual of the same class have same contributions to the          ###
###    population at all times                                                  ###
###################################################################################

"addClasses" <- function(phen, breed, t){
  phen$Age   <- t+1-floor(phen$Born)
  BreedNames <- names(breed)
  phen$Class <- paste(phen$Breed, phen$Age, sep=".Age")
  for(b in BreedNames){
    use <- phen$Breed==b
    if(length(breed[[b]]$Groups)==2){
      phen$Class[use] <- paste(phen$Class[use], phen$Sex[use], sep=".")
    }else{
      phen$Class[use] <- paste(phen$Class[use], "cohort",      sep=".")
      cat(paste0("Breed ", b, " has NA-sexes, so sexes are ignored for this breed.\n"))
    }
  }
  
  phen
}