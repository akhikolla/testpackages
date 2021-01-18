library(MixAll)
## get data and target from iris data set
data(iris)
x <- as.matrix(iris[,1:4]); z <- as.vector(iris[,5]); n <- nrow(x); p <- ncol(x)
## add missing values at random
indexes <- matrix(c(round(runif(5,1,n)), round(runif(5,1,p))), ncol=2)
cbind(indexes, x[indexes])
x[indexes] <- NA
## learn continuous model
model <- learnDiagGaussian( data=x, labels= z, prop = c(1/3,1/3,1/3)
                          , models = clusterDiagGaussianNames(prop = "equal")
                          , algo = "simul", nbIter = 2, epsilon = 1e-08
                          )
missingValues(model)
print(model)
model <- learnDiagGaussian( data=x, labels= z,
                          , models = clusterDiagGaussianNames(prop = "equal")
                          , algo = "impute", nbIter = 2, epsilon = 1e-08)
missingValues(model)
print(model)
model <- learnGamma( data=x, labels= z,
                   , models = clusterGammaNames(prop = "equal")
                   , algo = "simul", nbIter = 2, epsilon = 1e-08
                   )
missingValues(model)
print(model)

## get data and target from DebTrivedi data set
data(DebTrivedi)
x <- DebTrivedi[, c(1, 6,8, 15)]
z <- DebTrivedi$medicaid;
n <- nrow(x); p <- ncol(x);
model <- learnPoisson( data=x, labels=z
                     , models = clusterPoissonNames(prop = "equal")
                     , algo="simul", nbIter = 2, epsilon =  1e-08
                     )
print(model)

## get data and target from bird data set
data(birds)

## add 10 missing values
x <- birds[,2:5]; x = as.matrix(x); z <- birds[,1]; n <- nrow(x); p <- ncol(x);
indexes <- matrix(c(round(runif(5,1,n)), round(runif(5,1,p))), ncol=2);
cbind(indexes, x[indexes])
x[indexes] <- NA;
model <- learnCategorical( data=x, labels=z
                         , models = clusterCategoricalNames(prop = "equal")
                         , algo="simul", nbIter = 2, epsilon =  1e-08
                         )
missingValues(model)
print(model)

## A quantitative example with the heart disease data set
data(HeartDisease.cat)
data(HeartDisease.cont)
data(HeartDisease.target)
## with default values
lcomponent = list(HeartDisease.cat, HeartDisease.cont);
models = c("categorical_pk_pjk","gaussian_pk_sjk")
z<-HeartDisease.target[[1]];
model <- learnMixedData(lcomponent, models, z, algo="simul", nbIter=2)
missingValues(model)
print(model)

