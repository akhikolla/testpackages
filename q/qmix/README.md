# qmix: An R Package for Finite Quantile Mixture Models


The package `qmix` contains functions to estimate finite quantile mixture models using MCMC methods. Both fixed- and random quantile specifications are allowed as pre-specified inputs to the estiamtion fucntion.  

Caution: The package is still under initial development and there is absolutely NO guarantee that the package is functional. You are using this package on your own risk!

## Installation

```r
# Make sure that the following packages have been installed in your local R environment
if(!require(rstan)) install.packages("rstan")

# Install cirque from github
if(!require(devtools)) install.packages("devtools")
devtools::install_github("xiao-lu-research/qmix")
```


## Usage

```r

# Load the package
library(qmix)

# Get help
?qmix

# simulate a mixture of 2 ALDs
k <- 2
N <- 50
beta1 <- -10
beta2 <- 10
set.seed(34324)
x1 <- rnorm(N,0,1)
x2 <- rnorm(N,0,1)
xb1 <- x1*beta1
xb2 <- x2*beta2
y1 <- y2 <- NA
p1 <- 0.1
p2 <- 0.9
for (i in 1:N){
y1[i] <- rald(1,mu = xb1[i],p = p1,sigma = 1)
y2[i] <- rald(1,mu = xb2[i],p = p2,sigma = 1)
}
y <- c(y1,y2)
x <- c(x1,x2)
dat <- as.data.frame(cbind(y,x))
dat$z = rnorm(N)

# Estimate the models using both the fixed- and random-quantile specification
model1 <- qmix(y ~ x+z, data = dat, nmix = 2, design = "fixed", q = c(0.1, 0.9))
model2 <- qmix(y ~ x+z, data = dat, nmix = 2, design = "random")

# Summarize the results
coef(model1)
coef(model2)

print(model1)
print(model2)

# check traceplots
plot(model1)
plot(model2)

```

## References

Lu, Xiao (2019). Beyond the Average: Conditional Hypothesis Testing with Quantile Mixture. Working Paper. 
