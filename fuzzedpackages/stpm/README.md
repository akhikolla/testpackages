# Stochastic Process Model (SPM)
## Features
### Data simulation
* Continuous (one- and multiple-dimensions)
* Discrete (one- and multiple-dimensions)

### Optimisation
* Continuous (one- and multiple-dimensions)
* Discrete (one- and multiple-dimensions)
* Time-dependant coefficients (one-dimensional optimisation)

### Data imputation (censored time-to-event data imputation)

### How to install

**Note:** for Windows, please install Rtools: https://cran.r-project.org/bin/windows/Rtools/
**Note:** compilation under Windows may fail if you run Windows on a virtual machine!
Then:
```
install.packages("devtools")
library(devtools)
install_github("izhbannikov/spm")
```

## Examples

### Test discrete simulation
```
library(stpm)
## Test discrete
data <- simdata_discr(N=1000)
pars <- spm_discrete(data)
pars
```

### Test projection

```
library(stpm)
# Setting up the model
model.par <- list()
model.par$a <- matrix(c(-0.05, 1e-3, 2e-3, -0.05), nrow=2, ncol=2, byrow=TRUE)
model.par$f1 <- matrix(c(90, 35), nrow=1, ncol=2)
model.par$Q <- matrix(c(1e-8, 1e-9, 1e-9, 1e-8), nrow=2, ncol=2, byrow=TRUE)
model.par$f <- matrix(c(80, 27), nrow=1, ncol=2)
model.par$b <- matrix(c(6, 2), nrow=2, ncol=2)
model.par$mu0 <- 1e-5
model.par$theta <- 0.11
# Projection
# Discrete-time model
data.proj.discrete <- spm_projection(model.par, N=100, ystart=c(80, 27), tstart=c(30, 60))
plot(data.proj.discrete$stat$srv.prob, xlim = c(30, 115))
# Continuous-time model
data.proj.continuous <- spm_projection(model.par, N=100, ystart=c(80, 27), model="continuous", gomp=TRUE)
plot(data.proj.continuous$stat$srv.prob, xlim = c(30, 115))
```

### Time-dependent model
```
library(stpm)
model.par <- list(at = "-0.05", f1t = "80", Qt = "2e-8", ft= "80", bt = "5", mu0t = "1e-5*exp(0.11*t)")
data.proj.time_dependent <- spm_projection(model.par, N=100, ystart=80, model="time-dependent")
plot(data.proj.time_dependent$stat$srv.prob, xlim = c(30,105))
```

### Test prepare_data()
```
library(stpm)
data <- prepare_data(x=system.file("data","longdat.csv",package="stpm"))
head(data[[1]])
head(data[[2]])
```

### Test sim_pobs()

```
library(stpm)
dat <- sim_pobs(N=500)
head(dat)
```

### Test spm_pobs()
```
library(stpm)
#Reading the data:
data <- sim_pobs(N=100)
head(data)
#Parameters estimation:
pars <- spm_pobs(x=data)
pars
```

### Test simdata_cont()
```
library(stpm)
dat <- simdata_cont(N=50)
head(dat)
```

### Test simdata_discr()
```
library(stpm)
data <- simdata_discr(N=100)
head(data)
```

### Test simdata_time-dep()
```
library(stpm)
dat <- simdata_time_dep(N=100)
head(dat)
```

### Test spm_continuous()
```
library(stpm)
data <- simdata_cont(N=50)
pars <- spm_continuous(dat=data,a=-0.05, f1=80, Q=2e-8, f=80, b=5, mu0=2e-5)
pars
```

### Test spm_discrete()
```
library(stpm)
data <- simdata_discr(N=10)
pars <- spm_discrete(data)
pars
```

### Test spm()
```
library(stpm)
data.continuous <- simdata_cont(N=1000)
data.discrete <- simdata_discr(N=1000)
data <- list(data.continuous, data.discrete)
p.discr.model <- spm(data)
p.discr.model
p.cont.model <- spm(data, model="continuous")
p.cont.model
p.td.model <- spm(data, model="time-dependent",formulas=list(at="aa*t+bb", f1t="f1", Qt="Q", ft="f", bt="b", mu0t="mu0"), start=list(a=-0.001, bb=0.05, f1=80, Q=2e-8, f=80, b=5, mu0=1e-3))
p.td.model
```

### Multiple imputation with spm.impute(...)

The SPM offers longitudinal data imputation with results that are better than from other imputation tools since it preserves data structure, i.e. relation between 
Y(t) and mu(Y(t),t). Below there are two examples of multiple data imputation with function spm.impute(...).

```
library(stpm)

#######################################################
############## One dimensional case ###################
#######################################################

# Data preparation (short format)#
data <- simdata_discr(N=1000, dt = 2, format="short")

miss.id <- sample(x=dim(data)[1], size=round(dim(data)[1]/4)) # ~25% missing data
incomplete.data <- data
incomplete.data[miss.id,4] <- NA
# End of data preparation #

##### Multiple imputation with SPM #####
imp.data <- spm.impute(x=incomplete.data, id=1, case="xi", t1=3, covariates="y1", minp=1, theta_range=seq(0.075, 0.09, by=0.001))$imputed

##### Look at the incomplete data with missings #####
head(incomplete.data)

##### Look at the imputed data #####
head(imp.data)


#########################################################
################ Two-dimensional case ###################
#########################################################

# Parameters for data simulation #
a <- matrix(c(-0.05, 0.01, 0.01, -0.05), nrow=2)
f1 <- matrix(c(90, 30), nrow=1, byrow=FALSE)
Q <- matrix(c(1e-7, 1e-8, 1e-8, 1e-7), nrow=2)
f0 <- matrix(c(80, 25), nrow=1, byrow=FALSE)
b <- matrix(c(5, 3), nrow=2, byrow=TRUE)
mu0 <- 1e-04
theta <- 0.07
ystart <- matrix(c(80, 25), nrow=2, byrow=TRUE)

# Data preparation #
data <- simdata_discr(N=1000, a=a, f1=f1, Q=Q, f=f0, b=b, ystart=ystart, mu0 = mu0, theta=theta, dt=2, format="short")

# Delete some observations in order to have approx. 25% missing data
incomplete.data <- data
miss.id <- sample(x=dim(data)[1], size=round(dim(data)[1]/4)) 
incomplete.data <- data
incomplete.data[miss.id,4] <- NA
miss.id <- sample(x=dim(data)[1], size=round(dim(data)[1]/4)) 
incomplete.data[miss.id,5] <- NA
# End of data preparation #

##### Multiple imputation with SPM #####
imp.data <- spm.impute(x=incomplete.data, id=1, case="xi", t1=3, covariates=c("y1", "y2"), minp=1, theta_range=seq(0.060, 0.07, by=0.001))$imputed

##### Look at the incomplete data with missings #####
head(incomplete.data)

##### Look at the imputed data #####
head(imp.data)


```