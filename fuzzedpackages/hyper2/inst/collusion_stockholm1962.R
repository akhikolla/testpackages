# Analysis of the 1962 world chess championship, Curacao.  Reference:
# C. C. Moul and J. V. C. Nye 2009. "Did the Soviets collude? A
# statistical analysis of championship chess 1940-1978.  Journal of
# Economic Behaviour and Organization 70(2009) 10-21


library("hyper2")

jj <- read.table("stockholm1962.txt",header=FALSE)
results <- as.matrix(jj[,-(1:2)])
players <- as.character(jj$V1)
nationality <- as.character(jj$V2)
rownames(results) <- players
colnames(results) <- players
points <- rowSums(results,na.rm=TRUE)


## Now some comparison between stockholm1962.txt (results) and
## stockholm1962_matches.txt (results2)
H <- hyper2(pnames=c("white","draw",players))
H_coll <- hyper2(pnames=c("white","draw","collusion",players))

restab <- read.table("stockholm1962_matches.txt",header=FALSE)
stopifnot(all(unique(sort(c(as.character(restab$V1),as.character(restab$V2)))) == sort(players)))
results2 <- matrix(0,length(players),length(players))
rownames(results2) <- players
colnames(results2) <- players
diag(results2) <- NA

for(i in seq_len(nrow(restab))){
    white_player <- as.character(restab[i,1])
    black_player <- as.character(restab[i,2])
    match_result <- as.character(restab[i,3])

    nw <- which(players==white_player)
    nb <- which(players==black_player)

    if(match_result == "1-0"){ # white win
        results2[nw,nb] <- 2
        results2[nb,nw] <- 0

        H[c(white_player             ,"white"       )] %<>% inc()
        H[c(white_player,black_player,"white","draw")] %<>% dec()

        if(nationality[nw]=="USSR" & nationality[nb]=="USSR"){
            H_coll[c(white_player             ,"white"       )] %<>% inc()
            H_coll[c(white_player,black_player,"white","collusion")] %<>% dec()
        } else {
            H_coll[c(white_player             ,"white"       )] %<>% inc()
            H_coll[c(white_player,black_player,"white","draw")] %<>% dec()
        }
    } else if(match_result == "0-1"){ # black wins
        results2[nw,nb] <- 0
        results2[nb,nw] <- 2
        H[c(             black_player               )] %<>% inc()
        H[c(white_player,black_player,"white","draw")] %<>% dec()

        if(nationality[nw]=="USSR" & nationality[nb]=="USSR"){
            H_coll[c(             black_player                    )] %<>% inc()
            H_coll[c(white_player,black_player,"white","collusion")] %<>% dec()
        } else {  # collusion not playing
            H_coll[c(             black_player               )] %<>% inc()
            H_coll[c(white_player,black_player,"white","draw")] %<>% dec()
        }

    } else if (match_result == "1/2-1/2"){ # draw
        results2[nw,nb] <- 1
        results2[nb,nw] <- 1
        H[c(                                  "draw")] %<>% inc()
        H[c(white_player,black_player,"white","draw")] %<>% dec()

        if(nationality[nw]=="USSR" & nationality[nb]=="USSR"){
            H_coll[c(                                  "collusion")] %<>% inc()
            H_coll[c(white_player,black_player,"white","collusion")] %<>% dec()
        } else {  # collusion not playing
            H_coll[c(                                  "draw")] %<>% inc()
            H_coll[c(white_player,black_player,"white","draw")] %<>% dec()
        }



    } else {
        stop("not possible")
    }
}

## results [from the square matrix] and results2 [from the 3-column
## dataframe] should match:
stopifnot(all(results[!is.na(results)] == results2[!is.na(results2)]))
small <- 0.01


maxlike_free <- maxp(H_coll,startp=c(small,small,small/2,rep(small,22)),give=TRUE)
print("l")
dput(maxlike_free$value)

jj <- maxlike_free$par
maxlike_free <- maxp(H_coll,startp=jj,give=TRUE)
dput(maxlike_free$value)

jj <- maxlike_free$par
maxlike_free <- maxp(H_coll,startp=jj,give=TRUE,hessian=TRUE)
dput(maxlike_free$value)


freemp <- maxlike_free$par
freemp[2:3] <- mean(freemp[3:2]) +c(small,-small)*0.7

startp <- freemp

for(i in 1:2){
    maxlike_constrained <-
        maxp(H_coll,startp=startp,fcm=c(0,1,-1,rep(0,22)),fcv=0,give=TRUE)
    startp <- maxlike_constrained$par
    dput(maxlike_constrained$value)
}

maxlike_constrained <-
    maxp(H_coll,startp=startp,fcm=c(0,1,-1,rep(0,22)),fcv=0,give=TRUE)

print(maxlike_free$value - maxlike_constrained$value)
print(paste("pvalue = ", pchisq(2*(maxlike_free$value - maxlike_constrained$value),df=1,lower.tail=FALSE)))
