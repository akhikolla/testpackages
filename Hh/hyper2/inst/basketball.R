library("hyper2")
data(NBA)
a <- NBA
allplayers <- as.matrix(a[,5:23])

H <- hyper2(pnames= c(colnames(allplayers),"C_possession","W_possession"))
## The two final names are 'ghost' players (in the sense of Hankin
## 2010) whose strength is added to that of the real players.


## players 1-9 Cleveland
## players 10-19 (sic) Warriors:
C_onpitch <- allplayers[,1:9]
W_onpitch <- allplayers[,10:19]

## go through each row:
for(i in seq_len(nrow(a))){
  if(a[i,4] == 'W'){  # Warriors score a point
    onpitch_scored  <- c(names(which(W_onpitch[i,])),"C_possession")
    onpitch_noscore <- names(which(C_onpitch[i,]))
  } else { # Cleveland score a point
    onpitch_scored  <- c(names(which(C_onpitch[i,])),"W_possession")
    onpitch_noscore <- names(which(W_onpitch[i,]))
    
    H[onpitch_scored] %<>% `+`(1)
    H[c(onpitch_scored,onpitch_noscore)] %<>% `-`(1)
  }
}

mH <- maxp(H)
dotchart(mH,col=c(rep("black",9),rep("red",10),rep("blue",2)),pch=16)
