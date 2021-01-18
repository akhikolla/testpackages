### R code from vignette source 'Learning-Mixtures.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: prelim
###################################################
library(MixAll)
MixAll.version <- packageDescription("MixAll")$Version
MixAll.date <- packageDescription("MixAll")$Date
set.seed(39)


###################################################
### code chunk number 2: Learning-Mixtures.Rnw:225-233
###################################################
data(iris);
x <- as.matrix(iris[,1:4]); z <- as.vector(iris[,5]); n <- nrow(x); p <- ncol(x);
indexes <- matrix(c(round(runif(5,1,n)), round(runif(5,1,p))), ncol=2);
cbind(indexes, x[indexes]) # store true values
x[indexes] <- NA;          # and set them as missing
model <- learnDiagGaussian( data=x, labels = z
                          , models = clusterDiagGaussianNames(prop = "equal"))
summary(model)


###################################################
### code chunk number 3: Learning-Mixtures.Rnw:238-243
###################################################
# get estimated missing vallues
missingValues(model)
# compare predictions with true values
table(model@zi,model@ziFit)
plot(model)


###################################################
### code chunk number 4: Learning-Mixtures.Rnw:251-261
###################################################
data(birds)
## add 10 missing values
x <- as.matrix(birds[,2:5]); z <- as.vector(birds[,1]); n <- nrow(x); p <- ncol(x);
indexes <- matrix(c(round(runif(5,1,n)), round(runif(5,1,p))), ncol=2)
cbind(indexes, x[indexes]) # print true values
x[indexes] <- NA;          # set them as missing
model <- learnCategorical( data=x, labels=z
                         , models = clusterCategoricalNames(prop = "equal")
                         , algo="simul", nbIter = 2)
summary(model)


###################################################
### code chunk number 5: Learning-Mixtures.Rnw:266-271
###################################################
# get estimated missing vallues
missingValues(model)
# compare predictions
table(model@zi,model@ziFit)
plot(model)


###################################################
### code chunk number 6: Learning-Mixtures.Rnw:280-292
###################################################
data(iris)
x <- as.matrix(iris[,1:4]); z <- as.vector(iris[,5]); n <- nrow(x); p <- ncol(x);
indexes <- matrix(c(round(runif(5,1,n)), round(runif(5,1,p))), ncol=2);
cbind(indexes, x[indexes]) # print true values
x[indexes] <- NA;          # set them as missing
model <- learnGamma( data=x, labels= z
                   , models = clusterGammaNames(prop = "equal")
                   , algo = "simul", nbIter = 2, epsilon = 1e-08
                )
summary(model)
# get estimated missing values
missingValues(model)


###################################################
### code chunk number 7: Learning-Mixtures.Rnw:297-300
###################################################
  # compare predictions
table(model@zi,model@ziFit)
plot(model)


###################################################
### code chunk number 8: Learning-Mixtures.Rnw:310-322
###################################################
data(DebTrivedi)
x <- DebTrivedi[, c(1, 6, 8, 15)]; z <- DebTrivedi$medicaid; n <- nrow(x); p <- ncol(x);
indexes <- matrix(c(round(runif(5,1,n)), round(runif(5,1,p))), ncol=2);
cbind(indexes, x[indexes]) # print true values
x[indexes] <- NA;          # set them as missing
model <- learnPoisson( data=x, labels=z
                     , models = clusterPoissonNames(prop = "equal")
                     , algo="simul", nbIter = 2, epsilon =  1e-08
)
summary(model)
# get estimated missing vallues
missingValues(model)


###################################################
### code chunk number 9: Learning-Mixtures.Rnw:327-330
###################################################
# compare predictions
table(model@zi,model@ziFit)
plot(model)


###################################################
### code chunk number 10: Learning-Mixtures.Rnw:362-370
###################################################
data(HeartDisease.cat)
data(HeartDisease.cont)
data(HeartDisease.target)
ldata = list(HeartDisease.cat, HeartDisease.cont);
models = c("categorical_pk_pjk","gaussian_pk_sjk")
z<-HeartDisease.target[[1]];
model <- learnMixedData(ldata, models, z, algo="simul", nbIter=2)
summary(model)


###################################################
### code chunk number 11: Learning-Mixtures.Rnw:375-380
###################################################
# get estimated missing values
missingValues(model)
# compare predictions
table(model@zi,model@ziFit)
plot(model)


