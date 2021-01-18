## ------------------------------------------------------------------------
library("hyper2",quietly=TRUE)

## ------------------------------------------------------------------------
H <- hyper2(pnames=c("S","a","b","c"))
H

## ------------------------------------------------------------------------
H[c("a","S")]     %<>% "+"(5  ) # a+S win 5 games
H["b"]            %<>% "+"(  1) # b wins 1 game
H[c("a","b","S")] %<>% "-"(5+1) # a+b+S play 1+5=6 games
H

## ------------------------------------------------------------------------
H[c("b","S")]     %<>% "+"(3  ) # b+S win 3 games
H["a"]            %<>% "+"(  1) # a wins 1 game 
H[c("a","b","S")] %<>% "-"(3+1) # a+b+S play 1+3=4 games

## ------------------------------------------------------------------------
H[c("a","S")]     %<>% "+"(4  ) # a+S win 4 games
H["c"]            %<>% "+"(  1) # c wins 1 game 
H[c("a","c","S")] %<>% "-"(4+1) # a+c+S play 1+4=5 games

## ------------------------------------------------------------------------
H["a"]            %<>% "+"(1  ) # a wins 1 game
H[c("c","S")]     %<>% "+"(  2) # c+S win 2 games
H[c("a","c","S")] %<>% "-"(1+2) # a+c+S play 1+2=3 games

## ------------------------------------------------------------------------
H

## ------------------------------------------------------------------------
p1 <- maxp(H)
p1

## ------------------------------------------------------------------------
maxlike_free <- loglik(H,indep(p1))
maxlike_free

## ------------------------------------------------------------------------
small <- 1e-4
maxlike_constrained <-
    maxp(
        H,
        startp=c(small/2, 1/3-small/2, 1/3-small/2),
        fcm=c(-1,0,0),  fcv=-small,give=TRUE
)

maxlike_constrained

## ------------------------------------------------------------------------
delta_support <- maxlike_free - maxlike_constrained$value
delta_support

## ------------------------------------------------------------------------
LR <- exp(delta_support)
LR

## ------------------------------------------------------------------------
Lambda <- 2*delta_support
pval <- pchisq(Lambda,df=1,lower.tail=FALSE)
pval 

## ------------------------------------------------------------------------
maxp(H)

## ------------------------------------------------------------------------
profile_likelihood <- function(S,give=FALSE){
    small <- 1e-2  # cannot make this too small, optimization routine needs wiggle room
    out <- maxp(
       H, startp=c(S+small/2 , small,small),
       give=TRUE,
       fcm=rbind(c(1,0,0),c(-1,0,0)),fcv=c(S,-S-small)  # S <= p[1] <= S+small
       )
    if(give){
       return(out)
      } else {
       return(out$value)
      }
}

## ------------------------------------------------------------------------
S_vals <- seq(from=0.02,to=0.85,len=30)  # possible values for S
prof_like <- sapply(S_vals,profile_likelihood)  
plot(S_vals, prof_like-max(prof_like),xlab="serve strength",ylab="profile likelihood")
abline(h= -2)

## ----dpi=72--------------------------------------------------------------
proflike <- function(bc,give=FALSE){
  B <- bc[1]
  C <- bc[2]
  if(B+C>=1){return(NA)}

  objective <- function(S){     #

    jj <- c(S,A=1-(S+B+C),B,C)  # B,C fixed
    loglik(H,indep(jj))         # no minus: we use maximum=T in optimize()
  }
  
  maxlike_constrained <-  # single DoF is S [B,C given, A is fillup]
    optimize(             # single DoF -> use optimize() not optim()
        f = objective,         
        interval = c(0,1-(B+C)),
        maximum = TRUE
    )

  if(give){
    return(maxlike_constrained)
  } else{
    return(maxlike_constrained$objective)
  }
} 

small <- 1e-5
n <- 50
p <- seq(from=small,to=1-small,length=n)
jj <- expand.grid(p,p)
support <- apply(jj,1,proflike)
support <- support-max(support,na.rm=TRUE)
supportx <- support
dim(supportx) <- c(n,n)

x <- jj[,1]+jj[,2]/2
y <- jj[,2]*sqrt(3)/2
plot(x,y,cex=(support+6)/12, 
   pch=16,asp=1,xlim=c(-0.1,1.1),ylim=c(-0.2,0.9),
   axes=FALSE,xlab="",ylab="",main="profile likelihood for B,C")
polygon(x=c(0,1/2,1),y=c(0,sqrt(3)/2,0))
text(0.07,-0.04,'B=0, C=0')
text(0.50,+0.90,'B=1, C=0')
text(0.93,-0.04,'B=0, C=1')

## ------------------------------------------------------------------------
logcontrast <- apply(jj,1,function(x){log(x[1]/x[2])})
plot(logcontrast,support,xlim=c(-5,5),ylim=c(-4,0))
abline(h=-2)

