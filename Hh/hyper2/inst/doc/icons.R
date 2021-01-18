## ------------------------------------------------------------------------
M <- matrix(c(
    5 , 3 , NA,  4, NA,   3,
    3 , NA,  5,  8, NA,   2,
    NA,  4,  9,  2, NA,   1,
    10,  3, NA,  3,  4,  NA,
    4 , NA,  5,  6,  3,  NA,
    NA,  4,  3,  1,  3,  NA,
    5 ,  1, NA, NA,  1,   2,
    5 , NA,  1, NA,  1,   1,
    NA,  9,  7, NA,  2,   0)
  , byrow=TRUE,ncol=6)
colnames(M) <- c("NB","L","PB","THC","OA","WAIS")
M

## ------------------------------------------------------------------------
library("hyper2",quietly=TRUE)
icons <- saffy(M)
icons

## ------------------------------------------------------------------------
mic <- maxp(icons)
mic
dotchart(mic,pch=16)

## ------------------------------------------------------------------------
L1 <- loglik(icons,indep(mic))
L1

## ------------------------------------------------------------------------
small <- 1e-6  #  ensure start at an interior point
maxlike_constrained <- 
    maxp(icons, startp=indep(equalp(icons))-c(small,0,0,0,0),
         give=TRUE, fcm=c(-1,0,0,0,0),fcv= -1/6)
maxlike_constrained	 

## ------------------------------------------------------------------------
ES <- L1-maxlike_constrained$value 
ES

## ------------------------------------------------------------------------
pchisq(2*ES,df=1,lower.tail=FALSE)

## ------------------------------------------------------------------------
o <- function(Ul,Cl,startp,give=FALSE){
    small <- 1e-6  #  ensure start at an interior point
    if(missing(startp)){startp <- small*(1:5)+rep(0.1,5)}			
    out <- maxp(icons, startp=small*(1:5)+rep(0.1,5), give=TRUE, fcm=Ul,fcv=Cl)
    if(give){
        return(out)
    }else{
        return(out$value)
    }
}

p2max <- o(c(-1, 1, 0, 0, 0), 0)
p3max <- o(c(-1, 0, 1, 0, 0), 0)
p4max <- o(c(-1, 0, 0, 1, 0), 0)
p5max <- o(c(-1, 0, 0, 0, 1), 0)
p6max <- o(c(-2,-1,-1,-1,-1),-1)

## ------------------------------------------------------------------------
likes <- c(p2max,p3max,p4max,p5max,p6max)
likes
ml <- max(likes) 
ml

## ------------------------------------------------------------------------
L1-ml

## ------------------------------------------------------------------------
o(c(-1, 1, 0, 0, 0), 0,give=TRUE)$par
o(c(-1, 0, 1, 0, 0), 0,give=TRUE)$par
o(c(-1, 0, 0, 1, 0), 0,give=TRUE)$par
o(c(-1, 0, 0, 0, 1), 0,give=TRUE)$par
o(c(-2,-1,-1,-1,-1),-1,give=TRUE)$par

## ------------------------------------------------------------------------
jj <- o(c(-1,-1,-1,-1,0) , -2/3, give=TRUE,start=indep((1:6)/21))$value
jj

## ------------------------------------------------------------------------
L1-jj

## ------------------------------------------------------------------------
start <- indep(c(small,small,small,small,0.5-2*small,0.5-2*small))
jj <- c(
   o(c(-1, 0, 0, 0, 1), 0,start=start),
   o(c( 0,-1, 0, 0, 1), 0,start=start),
   o(c( 0, 0,-1, 0, 1), 0,start=start),
   o(c( 0, 0, 0,-1, 1), 0,start=start),

   o(c(-2,-1,-1,-1,-1),-1,start=start),
   o(c(-1,-2,-1,-1,-1),-1,start=start),
   o(c(-1,-1,-2,-1,-1),-1,start=start),
   o(c(-1,-1,-1,-2,-1),-1,start=start)
   )
jj
max(jj)

## ------------------------------------------------------------------------
L1-max(jj)

## ------------------------------------------------------------------------
   o(c( 0, 0, 0,-1, 1), 0,give=TRUE,start=start)

