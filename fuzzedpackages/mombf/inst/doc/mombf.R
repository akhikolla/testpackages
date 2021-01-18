### R code from vignette source 'mombf.Rnw'

###################################################
### code chunk number 1: quickstart
###################################################
library(mombf)
set.seed(1234)
x <- matrix(rnorm(100*3),nrow=100,ncol=3)
theta <- matrix(c(1,1,0),ncol=1)
y <- x %*% theta + rnorm(100)

#Default MOM prior on parameters
priorCoef <- momprior(tau=0.348)

#Beta-Binomial prior for model space
priorDelta <- modelbbprior(1,1)

#Model selection
fit1 <- modelSelection(y ~ x[,1]+x[,2]+x[,3], priorCoef=priorCoef, priorDelta=priorDelta)

#Posterior model probabilities
postProb(fit1)

#BMA estimates, 95% intervals, marginal post prob
coef(fit1)

#BMA predictions for y, 95% intervals
ypred= predict(fit1)
head(ypred)
cor(y, ypred[,1])


###################################################
### code chunk number 2: priorp2g
###################################################
priorp2g(.01, q=.2, prior='normalMom')
2*pmom(-.2, tau=0.3483356)


###################################################
### code chunk number 3: priorplota
###################################################
 thseq <- seq(-3,3,length=1000)
 plot(thseq,dnorm(thseq),type='l',ylab='Prior density')
 lines(thseq,dt(thseq,df=1),lty=2,col=2)
 legend('topright',c('Normal','Cauchy'),lty=1:2,col=1:2)


###################################################
### code chunk number 4: priorplotb
###################################################
 thseq <- seq(-3,3,length=1000)
 plot(thseq,dmom(thseq,tau=.348),type='l',ylab='Prior density',ylim=c(0,1.2))
 lines(thseq,demom(thseq,tau=.119),lty=2,col=2)
 lines(thseq,dimom(thseq,tau=.133),lty=3,col=4)
 legend('topright',c('MOM','eMOM','iMOM'),lty=1:3,col=c(1,2,4))


###################################################
### code chunk number 5: fpriorplota
###################################################
 thseq <- seq(-3,3,length=1000)
 plot(thseq,dnorm(thseq),type='l',ylab='Prior density')
 lines(thseq,dt(thseq,df=1),lty=2,col=2)
 legend('topright',c('Normal','Cauchy'),lty=1:2,col=1:2)


###################################################
### code chunk number 6: fpriorplotb
###################################################
 thseq <- seq(-3,3,length=1000)
 plot(thseq,dmom(thseq,tau=.348),type='l',ylab='Prior density',ylim=c(0,1.2))
 lines(thseq,demom(thseq,tau=.119),lty=2,col=2)
 lines(thseq,dimom(thseq,tau=.133),lty=3,col=4)
 legend('topright',c('MOM','eMOM','iMOM'),lty=1:3,col=c(1,2,4))


###################################################
### code chunk number 7: varsel1
###################################################
set.seed(2011*01*18)
x <- matrix(rnorm(100*3),nrow=100,ncol=3)
theta <- matrix(c(1,1,0),ncol=1)
y <- x %*% theta + rnorm(100)


###################################################
### code chunk number 8: varsel2
###################################################
priorCoef <- imomprior(tau=.133)
priorDelta <- modelbbprior(alpha.p=1,beta.p=1)
priorVar <- igprior(.01,.01)


###################################################
### code chunk number 9: varsel3
###################################################
fit1 <- modelSelection(y=y, x=x, center=FALSE, scale=FALSE,
priorCoef=priorCoef, priorDelta=priorDelta, priorVar=priorVar)
fit1$postMode
fit1$margpp
postProb(fit1)

fit2 <- modelSelection(y=y, x=x, center=FALSE, scale=FALSE,
priorCoef=priorCoef, priorDelta=priorDelta, priorVar=priorVar,
enumerate=FALSE, niter=1000)
fit2$postMode
fit2$margpp
postProb(fit2,method='norm')
postProb(fit2,method='exact')


###################################################
### code chunk number 10: varsel4
###################################################
priorCoef <- zellnerprior(tau=nrow(x))
fit3 <- modelSelection(y=y, x=x, center=FALSE, scale=FALSE, niter=10^2,
priorCoef=priorCoef, priorDelta=priorDelta, priorVar=priorVar,
method='Laplace')
postProb(fit3)


###################################################
### code chunk number 11: varsel5
###################################################
priorCoef <- imomprior(tau=.133)
fit4 <- modelSelection(y=y, x=x, center=FALSE, scale=FALSE,
priorCoef=priorCoef, priorDelta=priorDelta, priorVar=priorVar,
priorSkew=imomprior(tau=.133),family='auto')
head(postProb(fit4))


###################################################
### code chunk number 12: haldlm
###################################################
data(hald)
dim(hald)
lm1 <- lm(hald[,1] ~ hald[,2] + hald[,3] + hald[,4] + hald[,5])
summary(lm1)
V <- summary(lm1)$cov.unscaled
diag(V)


###################################################
### code chunk number 13: three
###################################################
mombf(lm1,coef=2,g=0.348)
imombf(lm1,coef=2,g=.133,method='adapt')
imombf(lm1,coef=2,g=.133,method='MC',B=10^5)


###################################################
### code chunk number 14: rnlp
###################################################
priorCoef <- momprior(tau=.348)
priorDelta <- modelbbprior(alpha.p=1,beta.p=1)
fit1 <- modelSelection(y=y, x=x, center=FALSE, scale=FALSE,
priorCoef=priorCoef, priorDelta=priorDelta)
b <- coef(fit1)
head(b)
th <- rnlp(msfit=fit1,priorCoef=priorCoef,niter=10000)
colMeans(th)
head(th)


###################################################
### code chunk number 15: rnlponemodel
###################################################
tau= 0.348
shrinkage= nrow(x)*tau/(1+nrow(x)*tau)
V= shrinkage * solve(t(x) %*% x)
m= as.vector(shrinkage * V %*% t(x) %*% y)
phi= mean((y - x%*%m)^2)
th= rnlp(m=m,V=phi * V,priorCoef=momprior(tau=tau))
colMeans(th)


###################################################
### code chunk number 16: bmsortho
###################################################
set.seed(1)
p <- 200; n <- 210
x <- scale(matrix(rnorm(n*p),nrow=n,ncol=p),center=TRUE,scale=TRUE)
S <- cov(x)
e <- eigen(cov(x))
x <- t(t(x %*% e$vectors)/sqrt(e$values))
th <- c(rep(0,p-3),c(.5,.75,1)); phi <- 1
y <- x %*% matrix(th,ncol=1) + rnorm(n,sd=sqrt(phi))
priorDelta=modelbinomprior(p=1/p)
priorVar=igprior(0.01,0.01)

priorCoef=zellnerprior(tau=n)
pm.zell <-
postModeOrtho(y,x=x,priorCoef=priorCoef,priorDelta=priorDelta,
priorVar=priorVar,bma=TRUE)
head(pm.zell$models)

priorCoef=momprior(tau=0.348)
pm.mom <- postModeOrtho(y,x=x,priorCoef=priorCoef,priorDelta=priorDelta,
priorVar=priorVar,bma=TRUE)
head(pm.mom$models)


###################################################
### code chunk number 17: bmsorthopp
###################################################
par(mar=c(5,5,1,1))
nvars <- sapply(strsplit(as.character(pm.zell$models$modelid),split=','),length)
plot(nvars,pm.zell$models$pp,ylab=expression(paste("p(",gamma,"|y)")),
xlab=expression(paste("|",gamma,"|")),cex.lab=1.5,ylim=0:1,xlim=c(0,50))
sel <- pm.zell$models$pp>.05
text(nvars[sel],pm.zell$models$pp[sel],pm.zell$models$modelid[sel],pos=4)
nvars <- sapply(strsplit(as.character(pm.mom$models$modelid),split=','),length)
points(nvars,pm.mom$models$pp,col='gray',pch=17)
sel <- pm.mom$models$pp>.05
text(nvars[sel],pm.mom$models$pp[sel],pm.mom$models$modelid[sel],pos=4,col='gray')
legend('topright',c('Zellner','MOM'),pch=c(1,17),col=c('black','gray'),cex=1.5)


###################################################
### code chunk number 18: bmaortho
###################################################
par(mar=c(5,5,1,1))
ols <- (t(x) %*% y) / colSums(x^2)
plot(ols,pm.zell$bma$coef,xlab='Least squares estimate',
ylab=expression(paste('E(',beta[j],'|y)')),cex.lab=1.5,cex.axis=1.2,col=1)
points(ols,pm.mom$bma$coef,pch=3,col='darkgray')
legend('topleft',c('Zellner','MOM'),pch=c(1,3),col=c('black','darkgray'))


###################################################
### code chunk number 19: bmsorthopp
###################################################
par(mar=c(5,5,1,1))
nvars <- sapply(strsplit(as.character(pm.zell$models$modelid),split=','),length)
plot(nvars,pm.zell$models$pp,ylab=expression(paste("p(",gamma,"|y)")),
xlab=expression(paste("|",gamma,"|")),cex.lab=1.5,ylim=0:1,xlim=c(0,50))
sel <- pm.zell$models$pp>.05
text(nvars[sel],pm.zell$models$pp[sel],pm.zell$models$modelid[sel],pos=4)
nvars <- sapply(strsplit(as.character(pm.mom$models$modelid),split=','),length)
points(nvars,pm.mom$models$pp,col='gray',pch=17)
sel <- pm.mom$models$pp>.05
text(nvars[sel],pm.mom$models$pp[sel],pm.mom$models$modelid[sel],pos=4,col='gray')
legend('topright',c('Zellner','MOM'),pch=c(1,17),col=c('black','gray'),cex=1.5)


###################################################
### code chunk number 20: bmaortho
###################################################
par(mar=c(5,5,1,1))
ols <- (t(x) %*% y) / colSums(x^2)
plot(ols,pm.zell$bma$coef,xlab='Least squares estimate',
ylab=expression(paste('E(',beta[j],'|y)')),cex.lab=1.5,cex.axis=1.2,col=1)
points(ols,pm.mom$bma$coef,pch=3,col='darkgray')
legend('topleft',c('Zellner','MOM'),pch=c(1,3),col=c('black','darkgray'))


###################################################
### code chunk number 21: bmsblockdiag1
###################################################
set.seed(1)
p <- 100; n <- 110
blocksize <- 10
blocks <- rep(1:(p/blocksize),each=blocksize)
x <- scale(matrix(rnorm(n*p),nrow=n,ncol=p),center=TRUE,scale=TRUE)
S <- cov(x)
e <- eigen(cov(x))
x <- t(t(x %*% e$vectors)/sqrt(e$values))
Sblock <- diag(blocksize)
Sblock[upper.tri(Sblock)] <- Sblock[lower.tri(Sblock)] <- 0.5
vv <- eigen(Sblock)$vectors
sqSblock <- vv %*% diag(sqrt(eigen(Sblock)$values)) %*% t(vv)
for (i in 1:(p/blocksize)) x[,blocks==i] <- x[,blocks==i] %*% sqSblock
th <- rep(0,ncol(x))
th[blocks==1] <- c(rep(0,blocksize-3),c(.5,.75,1))
th[blocks==2] <- c(rep(0,blocksize-2),c(.75,-1))
phi <- 1
y <- x %*% matrix(th,ncol=1) + rnorm(n,sd=sqrt(phi))


###################################################
### code chunk number 22: bmsblockdiag2
###################################################
priorCoef=zellnerprior(tau=n)
priorDelta=modelbinomprior(p=1/p)
priorVar=igprior(0.01,0.01)
pm <- postModeBlockDiag(y=y,x=x,blocks=blocks,priorCoef=priorCoef,
priorDelta=priorDelta,priorVar=priorVar,bma=TRUE)

head(pm$models)
head(pm$postmean.model)


###################################################
### code chunk number 23: coolblock
###################################################
maxvars=50
ylim=range(pm$postmean.model[,-1])
plot(NA,NA,xlab='Model size',
  ylab='Posterior mean given model',
  xlim=c(0,maxvars),ylim=ylim,cex.lab=1.5)
visited <- which(!is.na(pm$models$pp))
for (i in 2:ncol(pm$postmean.model)) {
  lines(pm$models$nvars[visited],pm$postmean.model[visited,i])
}
text(maxvars, pm$postmean.model[maxvars,which(th!=0)+1],
paste('X',which(th!=0),sep=''), pos=3)


###################################################
### code chunk number 24: coolblock
###################################################
maxvars=50
ylim=range(pm$postmean.model[,-1])
plot(NA,NA,xlab='Model size',
  ylab='Posterior mean given model',
  xlim=c(0,maxvars),ylim=ylim,cex.lab=1.5)
visited <- which(!is.na(pm$models$pp))
for (i in 2:ncol(pm$postmean.model)) {
  lines(pm$models$nvars[visited],pm$postmean.model[visited,i])
}
text(maxvars, pm$postmean.model[maxvars,which(th!=0)+1],
paste('X',which(th!=0),sep=''), pos=3)


###################################################
### code chunk number 25: six
###################################################
set.seed(4*2*2008)
n <- 50; theta <- c(log(2),0)
x <- matrix(NA,nrow=n,ncol=2)
x[,1] <- rnorm(n,0,1); x[,2] <- rnorm(n,.5*x[,1],1)
p <- pnorm(x %*% matrix(theta,ncol=1))
y <- rbinom(n,1,p)


###################################################
### code chunk number 26: pmomLM
###################################################
th <- pmomLM(y=y,x=x,xadj=rep(1,n),niter=10000)
head(th$postModel)
table(apply(th$postModel,1,paste,collapse=','))


###################################################
### code chunk number 27: seven
###################################################
glm1 <- glm(y~x[,1]+x[,2],family=binomial(link = "probit"))
thetahat <- coef(glm1)
V <- summary(glm1)$cov.scaled
g <- .5
bfmom.1 <- momknown(thetahat[2],V[2,2],n=n,g=g,sigma=1)
bfimom.1 <- imomknown(thetahat[2],V[2,2],n=n,nuisance.theta=2,g=g,sigma=1)
bfmom.1
bfimom.1


###################################################
### code chunk number 28: nine
###################################################
bfmom.2 <- momknown(thetahat[3],V[3,3],n=n,g=g,sigma=1)
bfimom.2 <- imomknown(thetahat[3],V[3,3],n=n,nuisance.theta=2,g=g,sigma=1)
bfmom.2
bfimom.2


###################################################
### code chunk number 29: ten
###################################################
bfmom.0 <- momknown(thetahat[2:3],V[2:3,2:3],n=n,g=g,sigma=1)
bfimom.0 <- imomknown(thetahat[2:3],V[2:3,2:3],n=n,nuisance.theta=2,g=g,sigma=1)
bfmom.0
bfimom.0


###################################################
### code chunk number 30: eleven
###################################################
prior.prob <- rep(1/4,4)
bf <- c(bfmom.0,bfmom.1,bfmom.2,1)
pos.prob <- prior.prob*bf/sum(prior.prob*bf)
pos.prob


###################################################
### code chunk number 31: mixtures1
###################################################
set.seed(1)
n=200; k=1:3
x= rnorm(n)


###################################################
### code chunk number 32: mixtures2
###################################################
fit= bfnormmix(x=x,k=k,q=3,q.niw=1,B=10^4)
postProb(fit)


###################################################
### code chunk number 33: mixturesNIW
###################################################
mcmcout= postSamples(fit)
names(mcmcout)
names(mcmcout[[2]])
colMeans(mcmcout[[2]]$eta)

parest= coef(fit)
names(parest)
parest[[2]]


###################################################
### code chunk number 34: mixturesMOMIW
###################################################
w= postSamples(fit)[[2]]$momiw.weight
apply(mcmcout[[2]]$eta, 2, weighted.mean, w=w)
apply(mcmcout[[2]]$mu, 2, weighted.mean, w=w)


###################################################
### code chunk number 35: mixture2comp
###################################################
set.seed(1)
n=200; k=1:3; probs= c(1/2,1/2)
mu1= -1.5; mu2=  1.5
x= rnorm(n) + ifelse(runif(n)<probs[1],mu1,mu2)
fit= bfnormmix(x=x,k=k,q=3,q.niw=1,B=10^4)
postProb(fit)


