library(hyper2)


hotstart <- indep(maxp(collusion))

f <- function(LO,give=FALSE){
    small <- 1e-5

    eq <- hotstart
    m <- sum(eq[2:3])
    a <- exp(LO)
    eq[2:3] <- m/(1+a)*c(1,a)  # exact

    eq[2:3] <- eq[2:3] + c(-small,small)
    
    out <-
        maxp(collusion,
             startp=eq,
             fcm=c(0, -1, exp(LO),rep(0,22)),
             fcv=0,
             give=TRUE)

    if(give){return(out)}else{return(out$value)}
}

x <- seq(from= 1.5,to = 2.7,len=32)
jj <- sapply(x,f)

plot(x,jj-max(jj),type='b',xlab="log-odds",ylab="support")
abline(h=c(0,-2))
