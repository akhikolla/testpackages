## Tennis example on p15 of Hankin 2010


library(hyper2)


## write a throwaway function that accounts for the regular players:


f <- function(H){

  ## 1&3 vs 2&4, scoreline 4-4:
  H[c(1,3)] <- H[c(1,3)] + 4
  H[c(2,4)] <- H[c(2,4)] + 4
  H[1:4] <- H[1:4] - 8
  
  ## 1&4 vs 2&3, scoreline 6-7:
  H[c(1,4)] <- H[c(1,4)] + 6
  H[c(2,3)] <- H[c(2,3)] + 7
  H[1:4] <- H[1:4] - 13

  ## 1 vs 3, scoreline 10-14:
  H[1] <- H[1] + 10
  H[3] <- H[2] + 14
  H[c(1,3)] <- H[c(1,3)] - 24

  ## 2 vs 3, scoreline 12-14:
  H[2] <- H[2] + 12
  H[3] <- H[3] + 14
  H[c(2,3)] <- H[c(2,3)] - 26

  ## 1 vs 4, scoreline 10-14:
  H[1] <- H[1] + 10
  H[4] <- H[4] + 14
  H[c(1,4)] <- H[c(1,4)] - 24

  ## 2 vs 4, scoreline 11-10:
  H[2] <- H[2] + 11
  H[4] <- H[4] + 10
  H[c(2,4)] <- H[c(2,4)] - 21

  ## 3 vs 4, scoreline 13-13:
  H[3] <- H[3] + 13
  H[4] <- H[4] + 13
  H[c(3,4)] <- H[c(3,4)] - 26
  
  return(H)
  }

## First do the no-ghost case:
H <- hyper2(d=4)
H <- f(H)

## 1&2 vs 3&4, scoreline 9-2:
H[c(1,2)] <- 9
H[c(3,4)] <- 2
H[1:4] <- H[1:4] - 11
pnames(H) <- c(paste("P",1:4,sep=""))

## Now the ghost, which is player number 5:
Hg <- hyper2(d=5)

## 1&2 vs 3&4 (NB: includes ghost!), scoreline 9-2 (again):
Hg[c(1,2,5)] <- 9
Hg[c(3,4)] <- 2
Hg[1:5] <- Hg[1:5] - 11

Hg <- f(Hg)
pnames(Hg) <- c(paste("P",1:4,sep=""),"ghost")
