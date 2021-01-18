

"characterizeBreeds" <- function(phen){
  BreedNames <- unique(phen$Breed)
  breed <- setNames(vector("list", length(BreedNames)), BreedNames)

  for(b in BreedNames){
    breed[[b]] <- list(Groups        = if(any(is.na(phen$Sex[phen$Breed==b]))){c("cohort")}else{c("male","female")},
                       hasCandidates = any(phen$isCandidate[phen$Breed==b]))
  }
  breed
}