### demo_stohastic.R
###-----------------------------------------------------------------------------
### $Id$
###-----------------------------------------------------------------------------

### DESCRIPTION OF A DEMONSTRATION
###-----------------------------------------------------------------------------

## A demonstration with a simple example to see in action the partitioning of
## additive genetic values by paths (Garcia-Cortes et al., 2008; Animal)
##       
## STOHASTIC SIMULATION (sort of)
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
## Additive genetic values are provided, i.e., no inference is being done!
##
## The idea of this example is not to do extensive simulations, but just to have
## a simple example to see how the partitioning of additive genetic values works.

### SETUP
###-----------------------------------------------------------------------------

options(width=200)

## install.packages(pkg=c("truncnorm"), dep=TRUE)
library(package="truncnorm")

### EXAMPLE PEDIGREE
###-----------------------------------------------------------------------------

## Generation 0
 id0 <- c("01", "02", "03", "04", "05", "06", "07", "08")
fid0 <- mid0 <- rep(NA, times=length(id0))
  h0 <- rep(c(1, 2), each=4)
  g0 <- rep(0, times=length(id0))

## Generation 1
 id1 <- c("11", "12", "13", "14", "15", "16", "17", "18")
fid1 <- c("02", "02", "02", "02", "06", "06", "06", "06")
mid1 <- c("01", "01", "03", "04", "05", "05", "07", "08")
  h1 <- h0
  g1 <- rep(1, times=length(id1))

## Generation 2
 id2 <- c("21", "22", "23", "24", "25", "26", "27", "28")
fid2 <- c("13", "13", "13", "13", "13", "13", "13", "13")
mid2 <- c("11", "12", "14", "14", "15", "16", "17", "18")
  h2 <- h0
  g2 <- rep(2, times=length(id2))

## Generation 3
 id3 <- c("31", "32", "33", "34", "35", "36", "37", "38")
fid3 <- c("24", "24", "24", "24", "24", "24", "24", "24")
mid3 <- c("21", "21", "22", "23", "25", "26", "27", "28")
  h3 <- h0
  g3 <- rep(3, times=length(id3))

## Generation 4
 id4 <- c("41", "42", "43", "44", "45", "46", "47", "48")
fid4 <- c("34", "34", "34", "34", "34", "34", "34", "34")
mid4 <- c("31", "32", "32", "33", "35", "36", "37", "38")
  h4 <- h0
  g4 <- rep(4, times=length(id4))

## Generation 5
 id5 <- c("51", "52", "53", "54", "55", "56", "57", "58")
fid5 <- c("44", "44", "44", "44", "44", "44", "44", "44")
mid5 <- c("41", "42", "43", "43", "45", "46", "47", "48")
  h5 <- h0
  g5 <- rep(5, times=length(id4))

ped <- data.frame( id=c( id0,  id1,  id2,  id3,  id4,  id5),
                  fid=c(fid0, fid1, fid2, fid3, fid4, fid5),
                  mid=c(mid0, mid1, mid2, mid3, mid4, mid5),
                  loc=c(  h0,   h1,   h2,   h3,   h4,   h5),
                  gen=c(  g0,   g1,   g2,   g3,   g4,   g5))
ped$sex <- 2
ped[ped$id %in% ped$fid, "sex"] <- 1
ped$loc.gen <- with(ped, paste(loc, gen, sep="-"))

### SIMULATE ADDITIVE GENETIC VALUES - STOHASTIC
###-----------------------------------------------------------------------------

## --- Parameters of simulation ---

## Additive genetic mean in founders by location
mu1 <- 2
mu2 <- 0

## Additive genetic variance in population
sigma2 <- 1
sigma  <- sqrt(sigma2)

## Threshold value for Mendelian sampling for selection - only values above this
##  will be accepted in simulation
t <- 0

## Set seed for simulation
set.seed(seed=19791123)

## --- Start of simulation ---

ped$agv1 <- NA ## Scenario (trait) 1: No selection in the second location
ped$agv2 <- NA ## Scenario (trait) 2:    Selection in the second location

## Generation 0  - founders (for simplicity set their values to the mean of location)
ped[ped$gen == 0 & ped$loc == 1, c("agv1", "agv2")] <- mu1
ped[ped$gen == 0 & ped$loc == 2, c("agv1", "agv2")] <- mu2

## Generation 1+ - non-founders
for(i in (length(g0)+1):nrow(ped)) {
  ## Scenario (trait) 1: selection only in the first location
  if(ped[i, "loc"] == 1) {
    w <- rtruncnorm(n=1, mean=0, sd=sqrt(sigma2/2), a=t)
  } else {
    w <-      rnorm(n=1, mean=0, sd=sqrt(sigma2/2))
  }
  ped[i, "agv1"] <- round(0.5 * ped[ped$id %in% ped[i, "fid"], "agv1"] +
                          0.5 * ped[ped$id %in% ped[i, "mid"], "agv1"] +
                          w, digits=1)
  ## Scenario (trait) 2: selection in both locations
  if(ped[i, "loc"] == 2) {
    w <- rtruncnorm(n=1, mean=0, sd=sqrt(sigma2/2), a=t)
  } ## for location 1 take the same values as above
  ped[i, "agv2"] <- round(0.5 * ped[ped$id %in% ped[i, "fid"], "agv2"] +
                          0.5 * ped[ped$id %in% ped[i, "mid"], "agv2"] +
                          w, digits=1)
}

### PLOT INDIVIDUAL ADDITIVE GENETIC VALUES
###-----------------------------------------------------------------------------

par(mfrow=c(2, 1), bty="l", pty="m", mar=c(2, 2, 1, 1) + .1, mgp=c(0.7, 0.2, 0))

tmp <- ped$gen + c(-1.5, -0.5, 0.5, 1.5) * 0.1
with(ped, plot(agv1 ~ tmp, pch=c(19, 21)[ped$loc], ylab="Additive genetic value",
               xlab="Generation", main="Selection in location 1", axes=FALSE,
               ylim=range(c(agv1, agv2))))
axis(1, labels=FALSE, tick=FALSE); axis(2, labels=FALSE, tick=FALSE); box()
legend(x="topleft", legend=c(1, 2), title="Location", pch=c(19, 21), bty="n")

with(ped, plot(agv2 ~ tmp, pch=c(19, 21)[ped$loc], ylab="Additive genetic value",
               xlab="Generation", main="Selection in locations 1 and 2", axes=FALSE,
               ylim=range(c(agv1, agv2))))
axis(1, labels=FALSE, tick=FALSE); axis(2, labels=FALSE, tick=FALSE); box()
legend(x="topleft", legend=c(1, 2), title="Location", pch=c(19, 21), bty="n")

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
### demo_stohastic.R ends here
