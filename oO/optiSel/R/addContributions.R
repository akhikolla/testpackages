
###### This function appends ######################################################
### Column 'c0' with contributions of individuals to the population at time t   ###
### Column 'c1' with contributions of individuals to the population at time t+1 ###
###################################################################################
### Note:                                                                       ###
###  - A breed has BreedxAge classes if there are NA-values in column 'Sex',    ###
###    and BreedxAgexSex classes otherwise.                                     ###
###  - All individual of the same class have same contributions to the          ###
###    population at all times                                                  ###
###  - sum(phen$c0)=1    within each breed                                      ###
###  - sum(phen$c1)=1-ra within each breed if there are selection candidates    ###
###    for this breed, and phen$c1=phen$c0 otherwise.                           ###
###################################################################################

"addContributions" <- function(phen, ageClass){
  phen$c0 <- setNames(ageClass$rcont0/ageClass$n, ageClass$Class)[phen$Class]
  phen$c1 <- setNames(ageClass$rcont1/ageClass$n, ageClass$Class)[phen$Class]
  phen$c0[is.na(phen$c0)] <- 0
  phen$c1[is.na(phen$c1)] <- 0
  return(phen)
}