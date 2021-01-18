

########################################################################
### A data table is created with one row for each BreedxAgexGroup-class  
### containing the following columns:
###  - Breed: Breed name
###  - Group: contains either 'male'/'female' or 'cohort', 
###           depending on whether sexes are distinguished or not.
###  - Age:   Age cohort. The youngest individuals are in age cohort 1.
###  - Class: Class name
###  - cont0: The contribution each class has in the 
###           idealized population at time t
###  - cont1: The contribution each class has in the 
###           idealized population at time t+1 
###  - n:     Number of individuals in the data set from each class
###  - rcont0:The contribution each classes is able to have 
###           at time t (0 for classes without records)
###  - rcont1:The contribution each classes is able to have 
###           at time t+1 (0 for classes without records)
###  - propOff: The proportion of offspring individuals from a given 
###           BreedXGroup-class are expected have at a given age.
########################################################################



"characterizeClasses" <- function(cont, breed, phen){
  
  ### Create data table with columns Breed, Group, age and Class ###
  ageClass <- reshape2::melt(map(breed,"Groups"))
  ageClass <- apply(ageClass, 1, paste, collapse=".")
  ageClass <- expand.grid(Breed=ageClass, age=1:nrow(cont),KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  ageClass[, c("Group", "Breed")] <- tstrsplit(ageClass$Breed, ".", fixed=TRUE)
  ageClass$Class <- paste0(ageClass$Breed, ".Age", ageClass$age, ".", ageClass$Group)
  ageClass <- ageClass[order(ageClass$Breed, ageClass$Group, ageClass$age),]
  rownames(ageClass) <- ageClass$Class
  AxG <- paste0(ageClass$age, ".", ageClass$Group)
  
  ### Append column cont0 with theoretical contributions of ###
  ### classes to the population at time t0 ####################
  cont$cohort <- cont$male+cont$female
  
  X <- reshape2::melt(cont, 
            id.vars       = "age", 
            measure.vars  = c("male", "female", "cohort"))

  ageClass$cont0 <- setNames(c(X[["value"]]), paste0(X[["age"]], ".", X[["variable"]]))[AxG]

  ### Append column cont1 with theoretical contributions of ###
  ### classes to the population at time t1 ####################
  cont1       <- rbind(cont,0)[-1,]
  cont1$age   <- 1:nrow(cont1)
  
  X <- reshape2::melt(cont1, 
            id.vars       = "age", 
            measure.vars  = c("male", "female", "cohort"))
 
  ageClass$cont1 <- setNames(c(X[["value"]]), paste0(X[["age"]], ".", X[["variable"]]))[AxG]

  ### Append column n with the numbers of individuals in the data set ###
  ageClass$n <- tapply(phen$Class, phen$Class, length)[ageClass$Class]
  ageClass$n[is.na(ageClass$n)] <- 0

  BreedNames <- names(breed)
  for(b in BreedNames){
    if(min(ageClass$n[ageClass$age==1 & ageClass$Breed==b])==0){
      cat("Number of individuals in age class 1:\n")
      print(table(phen$Breed[phen$Age==1], phen$Sex[phen$Age==1]))
      stop(paste0("There is no ", b, " from a youngest age class in the data set."))
    }
  }
  
  ### Append column rcont0 with real contributions of ######
  ### classes to the population at time t0 ####################
  
  ageClass$rcont0 <- ifelse(ageClass$age==1, ageClass$cont0, 0)
  for(b in BreedNames){
    for(grp in breed[[b]]$Groups){
      use   <- ageClass$Breed==b & ageClass$Group==grp & ageClass$age>1
      delta <- sum(ageClass$cont0[use & ageClass$n==0])/sum(use & ageClass$n>0)
      ageClass$rcont0[use & ageClass$n>0] <- ageClass$cont0[use & ageClass$n>0]+delta
      }
  }

  ### Append column rcont1 with real contributions of ######
  ### classes to the population at time t1 ####################
  
  ageClass$rcont1 <- ifelse(ageClass$age==1, ageClass$cont1, 0)
  for(b in BreedNames){
    for(grp in breed[[b]]$Groups){
      use   <- ageClass$Breed==b & ageClass$Group==grp & ageClass$age>1
      delta <- sum(ageClass$cont1[use & ageClass$n==0])/sum(use & ageClass$n>0)
      ageClass$rcont1[use & ageClass$n>0] <- ageClass$cont1[use & ageClass$n>0]+delta
    }
  }
  
  #################
  
  #for(b in BreedNames){
  #  if(!breed[[b]]$hasCandidates){
  #    ageClass$cont1[ageClass$Breed==b]  <- ageClass$cont0[ageClass$Breed==b]
  #    ageClass$rcont1[ageClass$Breed==b] <- ageClass$rcont0[ageClass$Breed==b]
  #  }
  #}
  
  ### Append column propOff with the proportion of offspring ###
  ### an individual has at this age                          ###
  
  pOff        <- cont - rbind(cont,0)[-1,]
  pOff$age    <- 1:nrow(pOff)
  pOff$male   <- pOff$male/sum(pOff$male)
  pOff$female <- pOff$female/sum(pOff$female)
  pOff$cohort <- pOff$cohort/sum(pOff$cohort)
  X <- reshape2::melt(pOff, 
            id.vars       = "age", 
            measure.vars  = c("male", "female", "cohort"))  
 
  ageClass$propOff <- setNames(c(X[["value"]]), paste0(X[["age"]], ".", X[["variable"]]))[AxG]
  
  return(ageClass)
}



