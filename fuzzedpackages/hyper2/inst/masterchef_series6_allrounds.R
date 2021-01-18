library(hyper2)
rm(list=ls())

Hall <-
    hyper2(pnames=c("Amy", "Ben", "Brendan", "Brent", "Byron",
                    "Cecilia", "Colin", "Deepali", "Emelia", "Emily",
                    "Georgia", "Jamie", "Kira", "Laura", "Nick",
                    "Nicole", "Rachael", "Renae", "Sam", "Sarah",
                    "Scott", "Sean", "Steven", "Tash", "Tracy"))  # NB alphabetical order


# players left in after L7b

H <- hyper2(pnames= c("Amy", "Ben", "Brent", "Colin", "Emelia",
                      "Georgia", "Jamie", "Kira", "Laura", "Renae",
                      "Sarah", "Tash", "Tracy"))

allplayers <-
    c("Brent","Laura","Emelia","Jamie","Tracy","Ben","Amy", "Renae",
      "Sarah","Colin","Kira","Georgia","Tash", "Byron","Sam",
      "Rachael","Steven","Sean","Scott", "Emily","Nick","Nicole",
      "Deepali","Brendan","Cecilia"
      )  # NB game order: Brent won, Laura runner-up, through to
         # Cecilia, who withdrew


 ## variable 'doo' is a Boolean, with entries governing whether a
 ## particular round is included in L or not.


doo <- c(
    L1   = FALSE,
    L2a  = FALSE,
    L2b  = FALSE,
    L2c  = FALSE,
    L3a  = FALSE,
    L3b  = FALSE,
    L3c  = FALSE,
    L4a  = FALSE,
    L4b  = FALSE,
    L4c  = FALSE,
    L5a  = FALSE,
    L5b  = FALSE,
    L6a  = FALSE,
    L6b  = FALSE,
    L7a  = FALSE,
    L7b  = TRUE,
    L8a  = TRUE,
    L8b  = TRUE,
    L8c  = FALSE,  # second chance, can't include this with n=13 players
    L8d  = TRUE,
    L8e  = TRUE,
    L9a  = TRUE,
    L9b  = TRUE,
    L10a = TRUE,
    L10b = TRUE,
    L10c = TRUE,
    L11a = TRUE,
    L11b = TRUE,
    L11c = TRUE,
    L11d = TRUE,
    L11e = TRUE,
    L11f = TRUE,
    L12a = TRUE,
    L12b = TRUE,
    L12c = TRUE,
    L12d = TRUE,
    L12e = TRUE,
    L13a = TRUE,
    L13b = TRUE
    )
    
L <- list()  # overall list

## If a date is given in the comments to a multi-stage elimination
## dataset, then the date refers to the date on which a player was
## actually eliminated.  Note that a likelihood function may also
## include selection of contestants for pressure tests.

if(doo[["L1"]]){
L$week1 <-  # 8 May 2014
ggol(H,
     top3       = c("Laura","Jamie","Sean"),
     btm21      = c("Emelia", "Tracy", "Sarah", "Colin", "Kira",
                    "Georgia", "Tash", "Byron", "Scott", "Emily",
                    "Nick", "Nicole"),
     btm9       = c("Brent", "Ben", "Amy", "Renae", "Sam", "Steven"),
     btm3       = c("Rachael", "Deepali"),
     eliminated = c("Brendan")
     )
}

if(doo[["L2a"]]){
L$week2a <- # elimination 4 chosen 11 May; Deepali eliminated 12 May 2014
ggol(H,
            win = c("Sarah"),
             IN = c("Brent","Laura","Emelia","Tracy","Ben","Amy","Renae",
                    "Colin","Kira","Georgia","Tash","Byron","Sam","Rachael",
                    "Steven","Sean","Emily","Nicole"),
           btm4 = c("Jamie","Scott","Nick"),
     eliminated = c("Deepali")
     )

}

if(doo[["L2b"]]){  
L$week2b <- H  # 13 May 2014
L$week2b[c(    #  winning team (11 players)
    "Emelia","Jamie","Ben","Renae","Sarah","Colin","Kira","Georgia",
    "Sam","Sean","Nick"
)] <- 1

L$week2b[c( # winning team union losing team; 22 players
    "Brent","Laura","Emelia","Jamie","Tracy","Ben","Amy","Renae",
    "Sarah","Colin","Kira","Georgia","Tash","Byron","Sam","Rachael",
    "Steven","Sean","Scott","Emily","Nick","Nicole"
)] <- -1
}

if(doo[["L2c"]]){  # 14 May
L$week2c <- ggol(H,  # elimination, field of losing team from 13 May.
            team_lose  = c("Laura","Amy","Tash","Rachael","Emily"),
            btm6       = c("Brent","Tracy","Byron","Steven","Scott"),
            eliminated = c("Nicole")
            )
}

if(doo[["L3a"]]){
L$week3a <-  # Nick eliminated 19 May, 
    ggol(H,
         top3       = c("Laura","Tracy","Amy"),
         IN         = c("Brent","Emelia","Jamie","Ben","Renae","Colin",
                        "Kira","Georgia","Tash","Byron","Sam","Rachael",
                        "Steven","Scott","Emily"),
         btm3       = c("Sarah","Sean"),
         eliminated = c("Nick")
         )
}

if(doo[["L3b"]]){
L$week3b <- H # team challenge 21 May
L$week3b[c(  # Blue team wins:
    "Laura","Amy","Renae","Sarah","Colin","Kira",
    "Tash","Byron","Rachael","Scott"
    )] <- 1

L$week3b[c(  # winning team union losing team
    "Brent","Laura","Emelia","Jamie","Tracy","Ben", "Amy","Renae",
    "Sarah","Colin","Kira","Georgia", "Tash","Byron","Sam","Rachael",
    "Steven","Sean","Scott","Emily")] <- -1
}

if(doo[["L3c"]]){
L$week3c <- ggol(H, # Emily eliminated 22 May
            team_lose  = c("Emelia","Tracy","Ben","Sam"),
            btm6       = c("Brent","Jamie","Georgia"),
            btm3       = c("Steven","Sean"),
            eliminated = c("Emily")
            )
}

if(doo[["L4a"]]){
L$week4a <- ggol(H, # Scott eliminated 26 May
            top3       = c("Brent","Sarah","Tash"),
            IN         = c("Laura","Emelia","Tracy","Ben","Amy","Renae",
                           "Colin","Kira","Georgia","Byron","Sam",
                           "Rachael","Sean"),
            btm3       = c("Jamie","Steven"),
            eliminated = c("Scott")
            )
}

if(doo[["L4b"]]){
L$week4b <- H  # 
L$week4b[c(  # team challenge 28 May; Blue team wins
    "Laura","Emelia","Jamie","Tracy","Ben","Renae","Sarah",
    "Byron","Steven"
)] <- 1

L$week4b[c(  # winning team union losing team
    "Brent","Laura","Emelia","Jamie","Tracy","Ben","Amy",
    "Renae","Sarah","Colin","Kira","Georgia","Tash",
    "Byron","Sam","Rachael","Steven","Sean"
)] <- -1
}

if(doo[["L4c"]]){ # Sean eliminated 29 May; losing team chosen 28 May
L$week4c <- ggol(H,
             team_lose  = c("Amy","Colin","Georgia",
                           "Tash","Sam","Rachael"),
             btm3       = c("Brent","Kira"),
             eliminated = c("Sean")
            )
}

if(doo[["L5a"]]){  # Steven eliminated 2 Jun; pressure test 1 Jun
L$week5a <- ggol(H,
            top3       = c("Brent","Renae","Rachael"),
            IN         = c("Laura","Emelia","Jamie","Tracy",
                           "Amy","Colin","Kira","Georgia","Tash",
                           "Sam"),
            btm4       = c("Ben","Sarah","Byron"),
            eliminated = c("Steven")
            )
}

if(doo[["L5b"]]){
L$week5b <- H   
L$week5b[c(   # Team challenge  4 June; all players
    "Brent","Laura","Emelia","Jamie","Tracy",
    "Ben","Amy","Renae","Sarah","Colin","Kira",
    "Georgia","Tash","Byron","Sam","Rchael"
   )] <- -1 
    
L$week5b[c(  
    "Tracy","Sarah","Kira","Sam"  # Yellow team wins
    )] <- 1


L$week5b[c(  # green team came second
    "Laura","Amy","Georgia","Byron"
)] <- 1

L$week5b[c(   # all players left after winning team sits out (denominator)
    "Brent","Laura","Emelia","Jamie",
    "Ben","Amy","Renae","Colin",
    "Georgia","Tash","Byron","Rachael"
)] <- -1 

L$week5b[c( # third team
    "Emelia","Jamie","Renae","Colin"
)] <- 1

L$week5b[c(   # all players left after winning and second team sits
              # out (denominator)
    "Brent","Emelia","Jamie", "Ben","Renae","Colin",
    "Tash","Rachael"
   )] <- -1 

L$week5b <- ggol(L$week5b,  # Rachael eliminated 5 Jun
            btm4   = c("Brent","Ben","Tash"),
            eliminated = c("Rachel")
            )

}

if(doo[["L6a"]]){
L$week6a <- H  # Sam eliminated 9th Jun
L$week6a <- ggol(H,
            top3       = c("Amy","Colin","Georgia"),
            IN         = c("Brent","Emelia","Jamie","Tracy",
                           "Ben","Renae","Sarah","Kira","Tash"),
            btm3       = c("Laura","Byron"),
            eliminated = c("Sam")
            )
}

if(doo[["L6b"]]){
L$week6b <- H

L$week6b[c(   # teams allocated randomly; red team wins 11 Jun
    "Laura","Ben","Amy","Renea","Kira","Georgia","Tash"
)] <- 1

L$week6b[c(   # all players (denominator)
    "Brent","Laura","Emelia","Jamie","Tracy","Ben","Amy","Renea","Sarah",
    "Colin", "Kira","Georgia","Tash","Byron"
)] <- 1

L$week6b <- ggol(L$week6b,    # Sarah eliminated 12 June
            team_lose  = c("Brent","Emelia","Jamie","Byron"),
            btm3       = c("Tracy","Colin"),
            eliminated = c("Sarah")
            )
}

if(doo[["L7a"]]){  # Byron eliminated 16 Jun;  
L$week7a <- ggol(H,
            top3       = c("Laura", "Tracy","Ben"),
            IN         = c("Brent","Emelia","Jamie","Renae",
                           "Colin","Kira","Georgia"),
            btm3       = c("Amy","Tash"),
            eliminated = c("Byron")
            )
}

 ## up to this point, following players are eliminated: Cecilia,
 ## Brendan, Deepali, Nicole, Nick, Emily, Scott, Sean, Steven,
 ## Rachael, Sam, Byron.

 ## Players left in are: Tash, Georgia, Kira, Colin, Sarah, Renae,
 ## Amy, Ben, Tracy, Jamie, Emelia, Laura, Brent




if(doo[["L7b"]]){
L$week7b <- H  # Team challenge 18 June; 
L$week7b[c(  # red team wins +=1
    "Jamie","Tracy","Ben","Amy","Renae","Georgia"
)] <- 1

L$week7b[c(  # all players -=1
    "Brent","Laura","Emelia","Jamie","Tracy",
    "Ben","Amy","Renae","Colin","Kira",  # NB no Sarah
    "Georgia","Tash"
    )] <- -1

L$week7b <- ggol(L$week7b,  # Tash eliminated 19 Jun; no Laura (she used her immunity pin)
            team_lose  = c("Brent","Emelia","Colin","Kira"),
            eliminated = c("Tash")
    )
}

if(doo[["L8a"]]){
L$week8a <- ggol(H,  # Georgia eliminated 23 June
            win       = c("Tracy"),
            top3      = c("Brent","Emelia","Jamie"),
            IN        = c("Ben","Colin","Kira"),
            btm3      = c("Laura","Amy","Renae"),  # Renae a special case;
            eliminated= c("Georgia")
            )
}

if(doo[["L8b"]]){

jj <- H  # team challenge 
jj[c(   # winners += 1
    "Brent","Emelia","Jamie","Amy","Renae"
)] <- 1

jj[c(   # all players -= 1
    "Brent","Laura","Emelia","Jamie","Tracy","Ben","Amy","Renae","Sarah","Colin","Kira"
)] <- -1

L$week8b <- list(jj)
}

if(doo[["L8c"]]){
jj <- H    # second chance; Sarah won so she is back in the contest
jj[c("Sarah")] <- 1
jj[c(
    "Sarah", "Georgia","Tash","Byron","Sam","Rachael","Steven",
    "Sean","Scott","Emily","Nick","Nicole","Deepali","Brendan",
    "Cecilia"
)] <- -1

L$week8c <- list(jj)
}

if(doo[["L8d"]]){
jj <- H     # team  challenge 25 Jun
jj[c(  # Blue team wins;  +=1
    "Brent","Emelia","Jamie","Ben","Amy", "Renae"
)] <- +1

jj[c(  # all players -= 1
    "Brent","Laura","Emelia","Jamie","Ben",
    "Amy","Renae","Sarah","Colin","Kira"
)] <- -1

L$week8d <- list(jj)
}

if(doo[["L8e"]]){
jj <- H    # Tracy wins power apron 26 Jun
jj[c("Tracy")] <- 1  # winner
jj[c("Brent","Emelia","Jamie","Tracy","Amy","Renae")] <- -1
L$week8e <- list(jj)
}

if(doo[["L9a"]]){    # Kira eliminated 30 Jun
L$week9a <- ggol(H,  # NB: Tracy excluded due to the power apron
            IN         = c("Laura","Emelia","Amy","Renae","Sarah"),
            btm5       = c("Brent","Ben"),
            btm3       = c("Jamie","Colin"),
            eliminated = c("Kira")
            )
}

if(doo[["L9b"]]){
L$week9b <- H  # team challenge 2 Jul

L$week9b[c(  # all contestants -= 1
    "Brent","Laura","Emelia","Jamie","Tracy","Ben","Amy","Renae","Sarah","Colin"
    )] <- -1

L$week9b[c("Emelia","Jamie","Tracy","Ben","Sarah")] <- -1    # Red team wins; +=1

L$week9b <- ggol(L$week9b,   # Colin eliminated 3 July
            team_lose = c("Brent","Laura","Amy"),
            btm2 = "Renae",
            eliminated = c("Colin")
            )
}

if(doo[["L10a"]]){  # Sarah eliminated 7 July
L$week10a <- ggol(H,
             top3       = c("Laura","Ben","Amy"),
             IN         = c("Brent","Jamie"),
             btm4       = c("Renae"),
             btm3       = c("Emelia","Tracy"),
             eliminated = c("Sarah")
             )
}

if(doo[["L10b"]]){  # team challenge 9 July
jj <- H
jj[c("Laura","Jamie")] <- 1  # blue team wins; +=1
jj[c(                        # all players -=1
    "Brent","Laura","Emelia","Jamie",
    "Tracy","Ben","Amy","Renae"
)] <- -1

jj[c("Emelia","Amy")] <- 1   # yellow team came second
jj[c("Brent","Emelia","Tracy", # all remaining players sans blue
       "Ben", "Amy","Renae"
       )] <- -1 

jj[c("Brent","Tracy")] <- 1 # green team came third
jj[c("Brent","Tracy","Ben","Renae")] <- -1 # remaining players

L$week10b <- list(jj)
}

if(doo[["L10c"]]){  # Renae eliminated 10 Jul
L$week10c <- ggol(H,
             win        = c("Laura"),
             IN         = c("Brent","Tracy","Ben"),
             eliminated = c("Renae")
             )
}

if(doo[["L11a"]]){
L$week11a <- ggol(H,  # Heston week, 13 July part 1
             top3=c("Laura","Jamie","Ben"),
             IN = c("Brent","Emelia","Tracy")
             )
} 

if(doo[["L11b"]]){ # Heston week, 13 July part 2
L$week11b <- ggol(H,
             win  = c("Brent"),
             IN   = c("Laura","Emelia","Jamie","Ben","Amy"),
             lose = c("Tracy")
             )
}

if(doo[["L11c"]]){ # Heston week, 14 July
L$week11c <- ggol(H,
             win = c("Amy"),
             IN  = c("Brent","Laura","Emelia","Jamie"),
             lose = c("Ben")
             )
}

if(doo[["L11d"]]){
L$week11d <- ggol(H,  # Heston week, 15 July
             win  = c("Amy"),
             IN   = c("Laura","Emelia","Jamie"),
             lose = c("Brent")
             )
}

if(doo[["L11e"]]){  # Heston week, 16 July
L$week11e <- ggol(H,
             win  = c("Laura"),
             IN   = c("Emelia","Jamie"),
             lose = c("Amy")
             )
}

if(doo[["L11f"]]){
L$week11f <- ggol(H,  # Amy eliminated 17 July
                    win = c("Laura"),
                   top2 = c("Jamie"),
                   btm4 = c("Brent","Tracy","Ben"),
             eliminated = c("Amy")
             )
}

if(doo[["L12a"]]){
L$week12a <- ggol(H,  # Ben eliminated 20 July
                    win = c("Brent","Laura","Emelia","Tracy"),
                   btm2 = c("Jamie"),
             eliminated = c("Ben")
             )
}

if(doo[["L12b"]]){
L$week12b <- ggol(H, # Tracy eliminated 21 July
                    win = c("Laura"),
                   top3 = c("Emelia","Jamie"),
                   btm2 = c("Brent"),
             eliminated = c("Tracy")
             )
}

if(doo[["L12c"]]){  # three-round duel challenge, 22 July
jj <- H
jj[c("Laura")] <- 1
jj[c("Brent","Laura")] <- -1

jj[c("Emelia")] <- 1
jj[c("Emelia","Jamie")] <- -1

jj[c("Laura")] <- powers(jj[c("Laura")]) + 1
jj[c("Emelia","Laura")] <- -1
L$week12c <- list(jj)
}

if(doo[["L12d"]]){  # 23 July?
L$week12d <- ggol(H,
             win = c("Brent"),
             IN = c("Emelia","Jamie")
             )
}

if(doo[["L12e"]]){ # Jamie eliminated 24 Jul
L$week12e <- ggol(H,
             btm3       = c("Brent","Emelia"),
             eliminated = c("Jamie")
             )
}

if(doo[["L13a"]]){ # Emelia eliminated 27 Jul
L$week13a <- ggol(H,
                    win = c("Laura"),
                   btm2 = c("Brent"),
             eliminated = c("Emelia")
             )
}

if(doo[["L13b"]]){
jj <- H    # Final 28 July; Laura eliminated and Brent wins
jj[c("Brent")] <- 1
jj[c("Brent","Laura")] <- -1
L$week13b <- list(jj)
}


`rprop` <- function(n){
    out <- runif(n)
    out <- out/sum(out)
    out[-n]
}

n <- 13   # 13 players; now specify constraints:
UI <- rbind(diag(n-1),-1)
CI <- c(rep(0,n-1),-1)

startp <- function(n){rep(1/n,n-1)}

jj <- # hot-start
  c(0.108618206752482, 0.0745797046860904, 0.134355320426558,
    0.0281960563127382, 0.116976608268827, 6.85045511136332e-09,
    0.106541176599712, 0.0205579373409824, 0.275062072951276,
    0.0407064266223856, 0.0280389311011074, 1.14209384061803e-09)


`ans_unconstrained` <-   # takes about an hour to run without hotstart
constrOptim(
    theta = jj,
    f = function(p){-like_series(p,L)},  # 'L' created sequentially above
    grad=NULL,
    ui = UI, ci=CI,
    control=list(trace=100,maxit=100000)
)


`ans_unconstrained_made_earlier` <-
structure(list(par = c(0.108618206752482, 0.0745797046860904, 
0.134355320426558, 0.0281960563127382, 0.116976608268827, 6.85045511136332e-09, 
0.106541176599712, 0.0205579373409824, 0.275062072951276, 0.0407064266223856, 
0.0280389311011074, 1.14209384061803e-09), value = 66.1965202262011, 
    counts = 0, convergence = 0L, message = NULL, outer.iterations = 1L, 
    barrier.value = 0.000214030311354918), .Names = c("par", 
"value", "counts", "convergence", "message", "outer.iterations", 
"barrier.value"))


## first implement the restriction that Brent >= Laura:
UI = rbind(UI,c(0,0,1,0,0,0,0,0,-1,0,0,0))  # Brent - Laura >= 0
CI <- c(CI,0)

## second, find a consistent starting value, swapping Laura's strength for that of Brent:
swap  <- jj[3]  # Brent == jj[3]
jj[3] <- jj[9]  # Laura == jj[9]
jj[9] <- swap


jj2 <- # hot-start
c(0.109667788153559, 0.0785631722927578, 0.188668727599298, 0.029077960904546, 
0.124873146008887, 2.63812282535037e-09, 0.115674972706005, 0.0213623680862838, 
0.188668726563236, 0.0421144946907343, 0.0307905358378776, 1.18098133183577e-09
)

ans_constrained <- 
constrOptim(
    theta = jj2,
    f = function(p){-like_series(p,L)},  # 'L' created sequentially above
    grad=NULL,
    ui = UI, ci=CI,
    control=list(trace=100,maxit=10000000)
)

`ans_constrained_made_earlier` <- 
structure(list(par = c(0.109667788153559, 0.0785631722927578, 
0.188668727599298, 0.029077960904546, 0.124873146008887, 2.63812282535037e-09, 
0.115674972706005, 0.0213623680862838, 0.188668726563236, 0.0421144946907343, 
0.0307905358378776, 1.18098133183577e-09), value = 67.3764248634724, 
    counts = 0, convergence = 0L, message = NULL, outer.iterations = 1L, 
    barrier.value = 0.000219352800726824), .Names = c("par", 
"value", "counts", "convergence", "message", "outer.iterations", 
"barrier.value"))

