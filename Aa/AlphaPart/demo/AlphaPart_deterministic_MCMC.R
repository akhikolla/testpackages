### AlphaPart_deterministic_MCMC.R
###-----------------------------------------------------------------------------
### $Id$
###-----------------------------------------------------------------------------

### DESCRIPTION OF A DEMONSTRATION
###-----------------------------------------------------------------------------

## A demonstration with a simple example to see in action the partitioning of
## additive genetic values by paths (Garcia-Cortes et al., 2008; Animal)
##       
## DETERMINISTIC SIMULATION (sort of)
##
## We have two locations (1 and 2). The first location has individualss with higher
## additive genetic value. Males from the first location are imported males to the
## second location from generation 2/3. This clearly leads to genetic gain in the
## second location. However, the second location can also perform their own selection
## so the question is how much genetic gain is due to the import of genes from the
## first location and due to their own selection.
##
## Above scenario will be tested with a simple example of a pedigree bellow. Two
## scenarios will be evaluated: without or with own selection in the second location.
## Selection will always be present in the first location.
##
## Additive genetic values are provided, but also infered from the data using MCMC
## method to show how to use AlphaPart with MCMC results!
##
## The idea of this example is not to do extensive simulations, but just to have
## a simple example to see how the partitioning of additive genetic values works.

### SETUP
###-----------------------------------------------------------------------------

options(width=200)

### EXAMPLE PEDIGREE & SETUP OF MENDELIAN SAMPLING - "DETERMINISTIC"
###-----------------------------------------------------------------------------

## Generation 0
  id0 <- c("01", "02", "03", "04")
 fid0 <- mid0 <- rep(NA, times=length(id0))
   h0 <- rep(c(1, 2), each=2)
   g0 <- rep(0, times=length(id0))
  w10 <- c( 2, 2, 0, 0)
  w20 <- c( 2, 2, 0, 0)

## Generation 1
  id1 <- c("11", "12", "13", "14")
 fid1 <- c("01", "01", "03", "03")
 mid1 <- c("02", "02", "04", "04")
   h1 <- h0
   g1 <- rep(1, times=length(id1))
  w11 <- c( 0.6,  0.2, -0.6,  0.2)
  w21 <- c( 0.6,  0.2,  0.6,  0.2)

## Generation 2
  id2 <- c("21", "22", "23", "24")
 fid2 <- c("12", "12", "12", "12")
 mid2 <- c("11", "11", "13", "14")
   h2 <- h0
   g2 <- rep(2, times=length(id2))
  w12 <- c( 0.6,  0.3, -0.2,  0.2)
  w22 <- c( 0.6,  0.3,  0.2,  0.2)

## Generation 3
  id3 <- c("31", "32", "33", "34")
 fid3 <- c("22", "22", "22", "22")
 mid3 <- c("21", "21", "23", "24")
   h3 <- h0
   g3 <- rep(3, times=length(id3))
  w13 <- c( 0.7,  0.1, -0.3,  0.3)
  w23 <- c( 0.7,  0.1,  0.3,  0.3)

## Generation 4
  id4 <- c("41", "42", "43", "44")
 fid4 <- c("32", "32", "32", "32")
 mid4 <- c("31", "31", "33", "34")
   h4 <- h0
   g4 <- rep(4, times=length(id4))
  w14 <- c( 0.8,  0.8, -0.1,  0.3)
  w24 <- c( 0.8,  0.8,  0.1,  0.3)

## Generation 5
  id5 <- c("51", "52", "53", "54")
 fid5 <- c("42", "42", "42", "42")
 mid5 <- c("41", "41", "43", "44")
   h5 <- h0
   g5 <- rep(5, times=length(id4))
  w15 <- c( 0.8,  1.0, -0.2,  0.3)
  w25 <- c( 0.8,  1.0,  0.2,  0.3)

ped <- data.frame( id=c( id0,  id1,  id2,  id3,  id4,  id5),
                  fid=c(fid0, fid1, fid2, fid3, fid4, fid5),
                  mid=c(mid0, mid1, mid2, mid3, mid4, mid5),
                  loc=c(  h0,   h1,   h2,   h3,   h4,  h5),
                  gen=c(  g0,   g1,   g2,   g3,   g4,  g5),
                   w1=c( w10,  w11,  w12,  w13,  w14,  w15),
                   w2=c( w20,  w21,  w22,  w23,  w24,  w25))
ped$sex <- 2
ped[ped$id %in% ped$fid, "sex"] <- 1
ped$loc.gen <- with(ped, paste(loc, gen, sep="-"))

### SIMULATE ADDITIVE GENETIC VALUES - SUM PARENT AVERAGE AND MENDELIAN SAMPLING
###-----------------------------------------------------------------------------

## Additive genetic mean in founders by location
mu1 <-  2
mu2 <-  0

## Additive genetic variance in population
sigma2 <- 1
sigma  <- sqrt(sigma2)

## Threshold value for Mendelian sampling for selection - only values above this
##  will be accepted in simulation
t <- 0

ped$agv1 <- ped$pa1 <- NA ## Scenario (trait) 1: No selection in the second location
ped$agv2 <- ped$pa2 <- NA ## Scenario (trait) 2:    Selection in the second location

## Generation 0  - founders (no parent average here - so setting it to zero)
ped[ped$gen == 0, c("pa1",  "pa2")] <- 0
ped[ped$gen == 0, c("agv1", "agv2")] <- ped[ped$gen == 0, c("w1", "w2")]

## Generation 1+ - non-founders (parent average + Mendelian sampling)
for(i in (length(g0)+1):nrow(ped)) {
  ped[i, "pa1"] <- 0.5 * (ped[ped$id %in% ped[i, "fid"], "agv1"] +
                          ped[ped$id %in% ped[i, "mid"], "agv1"])
  ped[i, "pa2"] <- 0.5 * (ped[ped$id %in% ped[i, "fid"], "agv2"] +
                          ped[ped$id %in% ped[i, "mid"], "agv2"])
  ped[i, c("agv1", "agv2")] <- ped[i, c("pa1", "pa2")] + ped[i, c("w1", "w2")]
}

## Recode identifications to 1:n
ped$ID  <- 1:nrow(ped)
ped$FID <- match(ped$fid, ped$id) 
ped$MID <- match(ped$mid, ped$id)

## pedigree object
ped2 <- with(ped, pedigree(sire=FID, dam=MID, label=ID))

### SIMULATE PHENOTYPES
###-----------------------------------------------------------------------------

mu     <- 100
h2     <- 0.9 ## set high heritability as this is a small example with little info
sigma2 <- (sigma2 / h2) - sigma2
sigma  <- sqrt(sigma2)

## Make repeated records to make pedigreemm() happy
dat <- rbind(ped, ped)
dat <- dat[order(dat$ID), ]

dat$y1 <- with(dat, round(rnorm(n=nrow(dat), mean=mu + agv1, sd=sigma), 1))
dat$y2 <- with(dat, round(rnorm(n=nrow(dat), mean=mu + agv2, sd=sigma), 1))

### PLOT INDIVIDUAL ADDITIVE GENETIC and PHENOTYPIC VALUES
###-----------------------------------------------------------------------------

par(mfrow=c(2, 2), bty="l", pty="m", mar=c(2, 2, 1, 1) + .1, mgp=c(0.7, 0.2, 0))

tmp <- ped$gen + seq(-1.5, 1.5, length=4) * 0.1
with(ped, plot(agv1 ~ tmp, pch=c(19, 21)[loc], ylab="Additive genetic value",
               xlab="Generation", main="Selection in location 1", axes=FALSE,
               ylim=range(c(agv1, agv2))))
axis(1, labels=FALSE, tick=FALSE); axis(2, labels=FALSE, tick=FALSE); box()
legend(x="topleft", legend=c(1, 2), title="Location", pch=c(19, 21), bty="n")

with(ped, plot(agv2 ~ tmp, pch=c(19, 21)[loc], ylab="Additive genetic value",
               xlab="Generation", main="Selection in locations 1 and 2", axes=FALSE,
               ylim=range(c(agv1, agv2))))
axis(1, labels=FALSE, tick=FALSE); axis(2, labels=FALSE, tick=FALSE); box()
legend(x="topleft", legend=c(1, 2), title="Location", pch=c(19, 21), bty="n")

tmp <- dat$gen + seq(-1.5, 1.5, length=8) * 0.1
with(dat, plot(y1 ~ tmp, pch=c(19, 21)[loc], ylab="Phenotypic value",
        xlab="Generation", main="Selection in location 1", axes=FALSE,
        ylim=range(c(y1, y2))))
axis(1, labels=FALSE, tick=FALSE); axis(2, labels=FALSE, tick=FALSE); box()
legend(x="topleft", legend=c(1, 2), title="Location", pch=c(19, 21), bty="n")

with(dat, plot(y2 ~ tmp, pch=c(19, 21)[loc], ylab="Phenotypic value",
        xlab="Generation", main="Selection in locations 1 and 2", axes=FALSE,
        ylim=range(c(y1, y2))))
axis(1, labels=FALSE, tick=FALSE); axis(2, labels=FALSE, tick=FALSE); box()
legend(x="topleft", legend=c(1, 2), title="Location", pch=c(19, 21), bty="n")

### INFER AGV VIA MCMC
###-----------------------------------------------------------------------------

## TODO: use MCMCglmm with fixed variances (make h2 high!)
fit1 <- pedigreemm(formula=y1 ~ (1 | ID), data=dat, pedigree=list(ID=ped2))
fit2 <- pedigreemm(formula=y2 ~ (1 | ID), data=dat, pedigree=list(ID=ped2))

par(mfrow=c(1, 2))
plot(ranef(fit1)$ID[, 1] ~ ped$agv1)
plot(ranef(fit2)$ID[, 1] ~ ped$agv2)
## TODO: is so bad correlation OK?

ranef(fit1, postVar=TRUE)
## Error in .local(object, ...) : 
##      code for applying pedigree and posterior variances not yet written
## TODO: check how postVar are obtained
## TODO: tell to Bates and Vasquez that this should not be needed as we need only diagonals

fit1MCMC <- mcmcsamp(fit1, n=1000, saveb=TRUE)
fit2MCMC <- mcmcsamp(fit2, n=1000, saveb=TRUE)

### PARTITION ADDITIVE GENETIC VALUES BY ...
###-----------------------------------------------------------------------------

## Compute partitions by location
(res <- AlphaPart(x=ped, colPath="loc", colBV=c("agv1", "agv2")))

## Summarize whole population
(ret <- summary(res))

## Summarize and plot by generation (=trend)
(ret <- summary(res, by="gen"))
plot(ret)

## Summarize and plot location specific trends
(ret <- summary(res, by="loc.gen"))
plot(ret)

## Summarize and plot location specific trends but only for location 1
(ret <- summary(res, by="loc.gen", subset=res[[1]]$loc == 1))
plot(ret)

## Summarize and plot location specific trends but only for location 2
(ret <- summary(res, by="loc.gen", subset=res[[1]]$loc == 2))
plot(ret)

## Compute partitions by location and sex
ped$loc.sex <- with(ped, paste(loc, sex, sep="-"))
(res <- AlphaPart(x=ped, colPath="loc.sex", colBV=c("agv1", "agv2")))

## Summarize and plot by generation (=trend)
(ret <- summary(res, by="gen"))
plot(ret)
plot(ret, lineTypeList=list("-1"=1, "-2"=2, def=3))

## Summarize and plot location specific trends
(ret <- summary(res, by="loc.gen"))
plot(ret, lineTypeList=list("-1"=1, "-2"=2, def=3))

## Summarize and plot location specific trends but only for location 1
(ret <- summary(res, by="loc.gen", subset=res[[1]]$loc == 1))
plot(ret)

## Summarize and plot location specific trends but only for location 2
(ret <- summary(res, by="loc.gen", subset=res[[1]]$loc == 2))
plot(ret)

###-----------------------------------------------------------------------------
### AlphaPart_deterministic_MCMC.R ends here
