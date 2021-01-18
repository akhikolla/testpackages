library("hyper2")

jj <- read.table("curacao1962.txt",header=FALSE)
results <- as.matrix(jj[,-(1:2)])
players <- as.character(jj$V1)
nationality <- as.character(jj$V2)
rownames(results) <- players
colnames(results) <- players

H <- hyper2(pnames=c("draw","collusion",players))
M <- matrix(0,nrow(results),ncol(results))

for(i in seq_len(nrow(results)-1)){  # i=row
  for(j in seq(from=i+1,to=ncol(results))){ # j = col
    M[i,j] <- 1
    r <- strsplit(results[i,j],'')[[1]]
    if((nationality[i]=="USSR") & (nationality[j]=="USSR")){
      monster <- "collusion"
    } else {
      monster <- "draw" 
    }
    if(any(r=="W")){  # row player wins
      wins <- sum(r=="W")
      H[  players[i]                    ] %<>% inc(wins)
      H[c(players[i],players[j],monster)] %<>% dec(wins)
    }
    if(any(r=="L")){  # row player loses
      losses <- sum(r=="L")
      H[             players[j]         ] %<>% inc(losses)
      H[c(players[i],players[j],monster)] %<>% dec(losses)
    }
    if(any(r=="W")){  # row player draws
      draws <- sum(r=="D")
      H[                        monster ] %<>% inc(draws)
      H[c(players[i],players[j],monster)] %<>% dec(draws)
    }
  }  # j loop closes
} # i loop closes


## Now some optimization.  First optimize freely:
max_support_free <- maxp(H,give=TRUE)$value
ml_p   <- maxp(H)

## Now optimize but constrained so collusion <= draw:

### First find a consistent start point, s:
s <- indep(equalp(H))
small <- 1e-3
s[1] %<>% `+`(small)
s[2] %<>% `-`(small)

## Perform the constrained optimization:
ml_p_constrained <- maxp(H,fcm=c(1,-1,rep(0,size(H)-3)),fcv=0,startp=s)
max_support_constrained <- loglik(H,indep(ml_p_constrained))

support <- max_support_free - max_support_constrained

print(paste("support = ", support,sep=""))
if(support>2){
  print("two units of support criterion exceeded: strong evidence that the Soviets colluded")
} else {
  print("less than two units of support: no evidence that the Soviets colluded")
}

print(paste("p-value = ",pchisq(2*support,df=1,lower.tail=FALSE)))



