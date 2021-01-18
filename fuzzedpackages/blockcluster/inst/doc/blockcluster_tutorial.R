### R code from vignette source 'blockcluster_tutorial.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: prelim
###################################################
library(blockcluster)
bc.version <- packageDescription("blockcluster")$Version
bc.date <- packageDescription("blockcluster")$Date


###################################################
### code chunk number 2: blockcluster_tutorial.Rnw:414-416
###################################################
defaultstrategy <- coclusterStrategy()
summary(defaultstrategy)


###################################################
### code chunk number 3: blockcluster_tutorial.Rnw:422-423
###################################################
newstrategy <- coclusterStrategy(nbtry=5, nbxem=10, algo='BCEM')


###################################################
### code chunk number 4: blockcluster_tutorial.Rnw:599-603
###################################################
library(blockcluster)
data("binarydata")
out<-coclusterBinary(binarydata, nbcocluster=c(2,3))
summary(out)


###################################################
### code chunk number 5: blockcluster_tutorial.Rnw:615-616
###################################################
plot(out, asp = 0)


###################################################
### code chunk number 6: blockcluster_tutorial.Rnw:625-626
###################################################
plot(out, type = 'distribution')


###################################################
### code chunk number 7: blockcluster_tutorial.Rnw:791-795
###################################################
data(binarydata)
out<-coclusterBinary(binarydata,nbcocluster=c(3,2), model="pik_rhol_epsilon")
summary(out)
plot(out)


###################################################
### code chunk number 8: blockcluster_tutorial.Rnw:800-803
###################################################
data(categoricaldata)
out<-coclusterCategorical(categoricaldata,nbcocluster=c(3,2))
summary(out)


###################################################
### code chunk number 9: blockcluster_tutorial.Rnw:808-811
###################################################
data(contingencydataunknown)
out<-coclusterContingency( contingencydataunknown, nbcocluster=c(2,3))
summary(out)


###################################################
### code chunk number 10: blockcluster_tutorial.Rnw:815-821
###################################################
data(contingencydataunknown)
mui= rep(1,nrow(contingencydataunknown)) 
nuj= rep(1,ncol(contingencydataunknown)) 
out<-coclusterContingency( list(contingencydataunknown, mui, nuj)
                         , nbcocluster=c(2,3), model="pik_rhol_known")
summary(out)


###################################################
### code chunk number 11: blockcluster_tutorial.Rnw:826-829
###################################################
data(gaussiandata)
out<-coclusterContinuous(gaussiandata,nbcocluster=c(2,3))
summary(out)


###################################################
### code chunk number 12: blockcluster_tutorial.Rnw:833-836
###################################################
data(gaussiandata)
out<-coclusterContinuous(gaussiandata,nbcocluster=c(2,3), model="pik_rhol_sigma2")
summary(out)


